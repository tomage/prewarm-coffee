use crate::physics::{HeatSnapshot, K_AIR};
use crate::storage::c_to_f;
use leptos::*;
use wasm_bindgen::JsCast;

// Profile chart coordinate constants (within SVG viewBox 0 0 500 240)
pub(crate) const PROFILE_X: f64 = 50.0;
pub(crate) const PROFILE_W: f64 = 430.0;
pub(crate) const PROFILE_Y_BOTTOM: f64 = 200.0;
pub(crate) const PROFILE_H: f64 = 180.0;

/// Look up the temperature at a given radius from a snapshot, applying
/// boundary-layer exponential decay for points beyond the outer surface.
fn temperature_at_radius(snap: &HeatSnapshot, r_m: f64, r_outer_m: f64, h_conv: f64, t_room: f64) -> f64 {
    if r_m <= r_outer_m {
        // Inside the mug: find nearest grid point
        snap.radii.iter().zip(snap.temperatures.iter())
            .min_by(|(r1, _), (r2, _)| {
                (*r1 - r_m).abs().total_cmp(&(*r2 - r_m).abs())
            })
            .map(|(_, &t)| t)
            .unwrap_or(t_room)
    } else {
        // Air region: interpolate surface temp then apply exponential decay
        let mut t_s = t_room;
        for i in 1..snap.radii.len() {
            if snap.radii[i - 1] <= r_outer_m && snap.radii[i] >= r_outer_m {
                let frac = (r_outer_m - snap.radii[i - 1]) / (snap.radii[i] - snap.radii[i - 1]);
                t_s = snap.temperatures[i - 1] + frac * (snap.temperatures[i] - snap.temperatures[i - 1]);
                break;
            }
        }
        let k_air = K_AIR;
        let delta = k_air * r_outer_m / (h_conv * r_outer_m + k_air);
        t_room + (t_s - t_room) * (-(r_m - r_outer_m) / delta).exp()
    }
}

/// Build SVG polyline points string for a temperature profile snapshot.
///
/// Plots actual grid data up to `r_outer`, then a smooth exponential boundary-layer
/// decay through the air region beyond it.
fn profile_polyline_points(
    snap: &HeatSnapshot,
    t_max: f64,
    t_min: f64,
    r_max_m: f64,
    r_outer_m: f64,
    h_conv: f64,
    t_room: f64,
) -> String {
    let t_range = t_max - t_min;

    // Find surface temperature at r_outer via interpolation
    let t_surface = {
        let mut t_s = t_room;
        for w in snap.radii.windows(2).zip(snap.temperatures.windows(2)) {
            let (rs, ts) = w;
            if rs[0] <= r_outer_m && rs[1] >= r_outer_m {
                let frac = (r_outer_m - rs[0]) / (rs[1] - rs[0]);
                t_s = ts[0] + frac * (ts[1] - ts[0]);
                break;
            }
        }
        t_s
    };
    let k_air = K_AIR;
    let delta = k_air * r_outer_m / (h_conv * r_outer_m + k_air);

    let mut pts: Vec<String> = Vec::new();
    for i in 0..snap.radii.len() {
        let r = snap.radii[i];
        if r > r_outer_m {
            break;
        }
        let t = snap.temperatures[i];
        let x = PROFILE_X + (r / r_max_m) * PROFILE_W;
        let y = PROFILE_Y_BOTTOM - ((t - t_min) / t_range * PROFILE_H).clamp(0.0, PROFILE_H);
        pts.push(format!("{:.1},{:.1}", x, y));
    }
    // Smooth boundary layer curve in air region
    let n_air_pts = 20;
    for j in 0..=n_air_pts {
        let frac = j as f64 / n_air_pts as f64;
        let r = r_outer_m + frac * (r_max_m - r_outer_m);
        let t_bl = t_room + (t_surface - t_room) * (-(r - r_outer_m) / delta).exp();
        let x = PROFILE_X + (r / r_max_m) * PROFILE_W;
        let y = PROFILE_Y_BOTTOM - ((t_bl - t_min) / t_range * PROFILE_H).clamp(0.0, PROFILE_H);
        pts.push(format!("{:.1},{:.1}", x, y));
    }
    pts.join(" ")
}

#[component]
pub fn ProfileChart(
    current_snapshot: Memo<HeatSnapshot>,
    initial_snapshot: Memo<HeatSnapshot>,
    inner_radius_cm: Memo<f64>,
    outer_radius_cm: Memo<f64>,
    chart_max_r_cm: Memo<f64>,
    temp_max: Memo<f64>,
    temp_min: Memo<f64>,
    convection_coeff: RwSignal<f64>,
    room_temp_c: RwSignal<f64>,
    use_fahrenheit: RwSignal<bool>,
    simulation_data_blended: Memo<Vec<HeatSnapshot>>,
    selected_time_idx: ReadSignal<usize>,
    set_selected_time_idx: WriteSignal<usize>,
) -> impl IntoView {
    // Hover state created internally
    let (profile_hover, set_profile_hover) = create_signal(None::<(f64, f64)>);

    view! {
        <div class="sim-viz">
            <h3>"Temperature Profile T(r)"</h3>
            <div class="time-slider">
                <label>"Time: " {move || {
                    let snap = current_snapshot.get();
                    let secs = snap.time_s;
                    if secs >= 60.0 {
                        format!("{:.0}m {:02.0}s", (secs / 60.0).floor(), secs % 60.0)
                    } else {
                        format!("{:.0}s", secs)
                    }
                }}</label>
                <input
                    type="range"
                    min="0"
                    max=move || simulation_data_blended.get().len().saturating_sub(1)
                    step="1"
                    prop:value=move || selected_time_idx.get()
                    on:input=move |ev| {
                        if let Ok(v) = event_target_value(&ev).parse::<usize>() {
                            set_selected_time_idx.set(v);
                        }
                    }
                />
            </div>

            <svg viewBox="0 0 500 240" class="profile-chart" role="img" aria-labelledby="profile-chart-title profile-chart-desc">
                <title id="profile-chart-title">"Temperature Profile T(r)"</title>
                <desc id="profile-chart-desc">"Radial temperature profile showing liquid, shell, and air regions at the selected time"</desc>
                // Background
                <rect x="50" y="20" width="430" height="180" fill="var(--color-bg)" rx="4"/>

                // Y-axis labels
                <text x="45" y="30" class="axis-label" text-anchor="end">
                    {move || { let v = temp_max.get(); if use_fahrenheit.get() { format!("{:.0}\u{00B0}F", c_to_f(v)) } else { format!("{:.0}\u{00B0}C", v) } }}
                </text>
                <text x="45" y="110" class="axis-label" text-anchor="end">
                    {move || { let v = (temp_max.get() + temp_min.get()) / 2.0; if use_fahrenheit.get() { format!("{:.0}\u{00B0}F", c_to_f(v)) } else { format!("{:.0}\u{00B0}C", v) } }}
                </text>
                <text x="45" y="195" class="axis-label" text-anchor="end">
                    {move || { let v = temp_min.get(); if use_fahrenheit.get() { format!("{:.0}\u{00B0}F", c_to_f(v)) } else { format!("{:.0}\u{00B0}C", v) } }}
                </text>

                // X-axis labels
                <text x="50" y="220" class="axis-label">"0"</text>
                <text x="240" y="220" class="axis-label" text-anchor="middle">"r"</text>
                <text x="475" y="220" class="axis-label" text-anchor="end">
                    {move || format!("{:.0}cm", chart_max_r_cm.get())}
                </text>

                // Liquid region background
                <rect
                    x="50"
                    y="20"
                    width=move || {
                        let inner_r = inner_radius_cm.get();
                        (inner_r / chart_max_r_cm.get() * PROFILE_W) as f64
                    }
                    height="180"
                    fill="var(--color-cold)"
                    opacity="0.2"
                />

                // Shell region background
                <rect
                    x=move || {
                        let inner_r = inner_radius_cm.get();
                        PROFILE_X + (inner_r / chart_max_r_cm.get() * PROFILE_W)
                    }
                    y="20"
                    width=move || {
                        let inner_r = inner_radius_cm.get();
                        let outer_r = outer_radius_cm.get();
                        ((outer_r - inner_r) / chart_max_r_cm.get() * PROFILE_W) as f64
                    }
                    height="180"
                    fill="var(--color-secondary)"
                    opacity="0.4"
                />

                // Inner boundary line
                <line
                    x1=move || {
                        let inner_r = inner_radius_cm.get();
                        PROFILE_X + (inner_r / chart_max_r_cm.get() * PROFILE_W)
                    }
                    y1="20"
                    x2=move || {
                        let inner_r = inner_radius_cm.get();
                        PROFILE_X + (inner_r / chart_max_r_cm.get() * PROFILE_W)
                    }
                    y2="200"
                    stroke="var(--color-cold)"
                    stroke-width="1"
                    stroke-dasharray="4 4"
                />
                <text
                    x=move || {
                        let inner_r = inner_radius_cm.get();
                        PROFILE_X + (inner_r / chart_max_r_cm.get() * PROFILE_W)
                    }
                    y="15"
                    class="axis-label"
                    text-anchor="middle"
                >"r"<tspan baseline-shift="sub" font-size="8">"i"</tspan></text>

                // Outer boundary line
                <line
                    x1=move || {
                        let outer_r = outer_radius_cm.get();
                        PROFILE_X + (outer_r / chart_max_r_cm.get() * PROFILE_W)
                    }
                    y1="20"
                    x2=move || {
                        let outer_r = outer_radius_cm.get();
                        PROFILE_X + (outer_r / chart_max_r_cm.get() * PROFILE_W)
                    }
                    y2="200"
                    stroke="var(--color-primary)"
                    stroke-width="1"
                    stroke-dasharray="4 4"
                />
                <text
                    x=move || {
                        let outer_r = outer_radius_cm.get();
                        PROFILE_X + (outer_r / chart_max_r_cm.get() * PROFILE_W)
                    }
                    y="15"
                    class="axis-label"
                    text-anchor="middle"
                >"r"<tspan baseline-shift="sub" font-size="8">"o"</tspan></text>

                // Initial profile (faded)
                <polyline
                    class="temp-profile-line initial"
                    points=move || {
                        profile_polyline_points(
                            &initial_snapshot.get(),
                            temp_max.get(), temp_min.get(),
                            chart_max_r_cm.get() / 100.0,
                            outer_radius_cm.get() / 100.0,
                            convection_coeff.get(), room_temp_c.get(),
                        )
                    }
                />

                // Current temperature profile
                <polyline
                    class="temp-profile-line"
                    points=move || {
                        profile_polyline_points(
                            &current_snapshot.get(),
                            temp_max.get(), temp_min.get(),
                            chart_max_r_cm.get() / 100.0,
                            outer_radius_cm.get() / 100.0,
                            convection_coeff.get(), room_temp_c.get(),
                        )
                    }
                />

                // Hover indicator
                {move || profile_hover.get().map(|(r_cm, temp)| {
                    let t_max = temp_max.get();
                    let t_min = temp_min.get();
                    let t_range = t_max - t_min;
                    let x = PROFILE_X + (r_cm / chart_max_r_cm.get()) * PROFILE_W;
                    let y = PROFILE_Y_BOTTOM - ((temp - t_min) / t_range * PROFILE_H).clamp(0.0, PROFILE_H);

                    view! {
                        <line x1=x y1="20" x2=x y2="200" class="hover-line"/>
                        <circle cx=x cy=y r="5" class="hover-dot"/>
                        <text x=x+10.0 y=y-10.0 class="hover-text">
                            {if use_fahrenheit.get() { format!("r={:.1}cm, T={:.1}\u{00B0}F", r_cm, c_to_f(temp)) } else { format!("r={:.1}cm, T={:.1}\u{00B0}C", r_cm, temp) }}
                        </text>
                    }
                })}

                // Legend
                <g class="chart-legend" transform="translate(405, 28)">
                    <rect x="-6" y="-6" width="72" height="60" fill="var(--color-surface)" opacity="0.85" rx="3"/>
                    <rect x="0" y="0" width="12" height="8" fill="var(--color-cold)" opacity="0.4"/>
                    <text x="16" y="7" class="legend-label">"Liquid"</text>
                    <rect x="0" y="14" width="12" height="8" fill="var(--color-secondary)" opacity="0.6"/>
                    <text x="16" y="21" class="legend-label">"Shell"</text>
                    <text x="0" y="35" fill="var(--color-text-muted)" font-size="9">"Air (outside)"</text>
                    <line x1="0" y1="44" x2="16" y2="44" stroke="var(--color-primary)" stroke-width="1" opacity="0.5"/>
                    <text x="20" y="47" class="legend-label">"t=0"</text>
                </g>

                // Interactive overlay
                <rect
                    x="0" y="0" width="500" height="240"
                    fill="transparent"
                    pointer-events="all"
                    style="cursor: crosshair"
                    on:mousemove=move |ev| {
                        let Some(target) = ev.target() else { return };
                        let rect_el = target.unchecked_into::<web_sys::Element>();
                        let Some(svg) = rect_el.parent_element() else { return };
                        let rendered_width = svg.client_width() as f64;
                        if rendered_width <= 0.0 { return; }

                        let vb_x = ev.offset_x() as f64 / rendered_width * 500.0;
                        let x_normalized = ((vb_x - PROFILE_X) / PROFILE_W).clamp(0.0, 1.0);
                        let r_cm = x_normalized * chart_max_r_cm.get();
                        let snap = current_snapshot.get();
                        let temp = temperature_at_radius(
                            &snap, r_cm / 100.0, outer_radius_cm.get() / 100.0,
                            convection_coeff.get(), room_temp_c.get(),
                        );
                        set_profile_hover.set(Some((r_cm, temp)));
                    }
                    on:mouseleave=move |_| {
                        set_profile_hover.set(None);
                    }
                />
            </svg>
        </div>

        <div class="sim-stats sim-stats--4col">
            <div class="stat liquid">
                <span class="stat-value">{move || { let v = current_snapshot.get().temperatures[0]; if use_fahrenheit.get() { format!("{:.1}\u{00B0}F", c_to_f(v)) } else { format!("{:.1}\u{00B0}C", v) } }}</span>
                <span class="stat-label">"Liquid center"</span>
            </div>
            <div class="stat liquid">
                <span class="stat-value">{move || {
                    let snap = current_snapshot.get();
                    let v = snap.temperatures[snap.inner_interface_idx];
                    if use_fahrenheit.get() { format!("{:.1}\u{00B0}F", c_to_f(v)) } else { format!("{:.1}\u{00B0}C", v) }
                }}</span>
                <span class="stat-label">"Inner surface"</span>
            </div>
            <div class="stat shell">
                <span class="stat-value">{move || {
                    let snap = current_snapshot.get();
                    let inner = snap.inner_interface_idx;
                    let outer = snap.outer_interface_idx;
                    let mid = (inner + outer) / 2;
                    let v = snap.temperatures[mid];
                    if use_fahrenheit.get() { format!("{:.1}\u{00B0}F", c_to_f(v)) } else { format!("{:.1}\u{00B0}C", v) }
                }}</span>
                <span class="stat-label">"Shell mid"</span>
            </div>
            <div class="stat shell">
                <span class="stat-value">{move || {
                    let snap = current_snapshot.get();
                    let v = snap.temperatures[snap.outer_interface_idx];
                    if use_fahrenheit.get() { format!("{:.1}\u{00B0}F", c_to_f(v)) } else { format!("{:.1}\u{00B0}C", v) }
                }}</span>
                <span class="stat-label">"Outer surface"</span>
            </div>
        </div>
    }
}

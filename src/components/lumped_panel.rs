use crate::components::{MobileSummaryBar, MugControls, SimulationResults};
use crate::components::profile_chart::{PROFILE_X, PROFILE_W, PROFILE_Y_BOTTOM, PROFILE_H};
use crate::components::solver_toggle::create_solver_signals;
use crate::physics::{
    compute_comparison, simulate_lumped_analytical, simulate_lumped_numerical,
    CoolingDataPoint, MugParameters, MAX_SIMULATION_MINUTES,
};
use crate::storage::{c_to_f, SharedSettings};
use leptos::*;
use std::rc::Rc;

/// Timed simulation result: data points + how long computation took
#[derive(Clone, PartialEq)]
struct TimedResult {
    data: Vec<CoolingDataPoint>,
    elapsed_ms: f64,
}

fn timed_numerical(params: &MugParameters, duration: usize) -> TimedResult {
    let start = js_sys::Date::now();
    let data = simulate_lumped_numerical(params, duration);
    TimedResult { data, elapsed_ms: js_sys::Date::now() - start }
}

fn timed_analytical(params: &MugParameters, duration: usize) -> TimedResult {
    let start = js_sys::Date::now();
    let data = simulate_lumped_analytical(params, duration);
    TimedResult { data, elapsed_ms: js_sys::Date::now() - start }
}

#[component]
pub fn LumpedPanel() -> impl IntoView {
    // Use shared settings singleton — same RwSignals across all routes
    let settings = SharedSettings::instance();
    let volume_ml = settings.volume_ml;
    let wall_thickness_mm = settings.wall_thickness_mm;
    let coffee_temp_c = settings.coffee_temp_c;
    let room_temp_c = settings.room_temp_c;
    let min_drinkable_c = settings.min_drinkable_c;
    let use_fahrenheit = settings.use_fahrenheit;
    let duration_minutes = settings.duration_minutes;
    let preheat_temp_c = settings.preheat_temp_c;

    let (use_analytical, set_use_analytical, sim_ready, set_sim_ready) = create_solver_signals(settings.use_analytical_lumped);

    let make_params = move |preheated: bool| settings.make_mug_params(preheated);

    // Timed simulation memos — compute data + elapsed time together
    let timed_cold = create_memo(move |_| {
        let params = make_params(false);
        if use_analytical.get() {
            timed_analytical(&params, MAX_SIMULATION_MINUTES)
        } else {
            timed_numerical(&params, MAX_SIMULATION_MINUTES)
        }
    });

    let timed_hot = create_memo(move |_| {
        let params = make_params(true);
        if use_analytical.get() {
            timed_analytical(&params, MAX_SIMULATION_MINUTES)
        } else {
            timed_numerical(&params, MAX_SIMULATION_MINUTES)
        }
    });

    // Total computation time (both cold + hot runs)
    let compute_time_ms = create_memo(move |_| {
        let total = timed_cold.get().elapsed_ms + timed_hot.get().elapsed_ms;
        set_sim_ready.set(true);
        total
    });

    // Convert CoolingDataPoint data to (time_minutes, temp_c) tuples
    let cold_coffee = create_memo(move |_| CoolingDataPoint::coffee_series(&timed_cold.get().data));
    let hot_coffee = create_memo(move |_| CoolingDataPoint::coffee_series(&timed_hot.get().data));
    let cold_mug = create_memo(move |_| CoolingDataPoint::mug_series(&timed_cold.get().data));
    let hot_mug = create_memo(move |_| CoolingDataPoint::mug_series(&timed_hot.get().data));
    let temp_top = create_memo(move |_| 100.0_f64);

    // Bundled comparison metrics
    let metrics = create_memo(move |_| {
        compute_comparison(&cold_coffee.get(), &hot_coffee.get(), min_drinkable_c.get())
    });

    // Crossing times for MobileSummaryBar
    let crossing_time_cold = create_memo(move |_| metrics.get().crossing_time_cold);
    let crossing_time_hot = create_memo(move |_| metrics.get().crossing_time_hot);

    // --- Profile chart signals ---
    let (selected_time_idx, set_selected_time_idx) = create_signal(0_usize);

    // Derived: inner radius from volume (V = 4/3 π r³ → r = (3V/4π)^(1/3))
    let inner_radius_cm = create_memo(move |_| {
        let v = volume_ml.get();
        (3.0 * v / (4.0 * std::f64::consts::PI)).cbrt()
    });

    // Derived: outer radius = inner + shell thickness
    let outer_radius_cm = create_memo(move |_| {
        inner_radius_cm.get() + wall_thickness_mm.get() / 10.0
    });

    let chart_max_r_cm = 10.0_f64;

    // Y-axis range
    let temp_max = create_memo(move |_| {
        let hi = coffee_temp_c.get().max(preheat_temp_c.get()).max(room_temp_c.get());
        (hi * 1.1 / 10.0).ceil() * 10.0
    });
    let temp_min = create_memo(move |_| {
        let lo = coffee_temp_c.get().min(preheat_temp_c.get()).min(room_temp_c.get());
        (lo * 0.9 / 10.0).floor() * 10.0
    });

    view! {
        <div class="model-container">
            <MugControls settings=settings />
            <SimulationResults
                cold_coffee=cold_coffee
                hot_coffee=hot_coffee
                cold_mug=cold_mug
                hot_mug=hot_mug
                metrics=metrics
                temp_top=temp_top
                room_temp=room_temp_c.read_only()
                min_drinkable=min_drinkable_c.read_only()
                duration_minutes=duration_minutes.read_only()
                use_fahrenheit=use_fahrenheit.read_only()
                use_analytical=use_analytical
                set_use_analytical=set_use_analytical
                compute_time_ms=compute_time_ms.into()
                sim_ready=Signal::from(sim_ready)
                max_sim_minutes=MAX_SIMULATION_MINUTES as f64
                footer=Rc::new(move || {
                    // Current time index clamped to data length
                    let time_idx = move || selected_time_idx.get().min(hot_coffee.get().len().saturating_sub(1));

                    // Temperatures at selected time
                    let coffee_t = move || hot_coffee.get().get(time_idx()).map(|&(_, t)| t).unwrap_or(0.0);
                    let mug_t = move || hot_mug.get().get(time_idx()).map(|&(_, t)| t).unwrap_or(0.0);
                    let room_t = move || room_temp_c.get();

                    // Initial temperatures (t=0)
                    let coffee_t0 = move || hot_coffee.get().first().map(|&(_, t)| t).unwrap_or(0.0);
                    let mug_t0 = move || hot_mug.get().first().map(|&(_, t)| t).unwrap_or(0.0);

                    // Time display
                    let time_display = move || {
                        let data = hot_coffee.get();
                        let idx = time_idx();
                        let minutes = data.get(idx).map(|&(t, _)| t).unwrap_or(0.0);
                        let secs = minutes * 60.0;
                        if secs >= 60.0 {
                            format!("{:.0}m {:02.0}s", (secs / 60.0).floor(), secs % 60.0)
                        } else {
                            format!("{:.0}s", secs)
                        }
                    };

                    // Helper: map temperature to SVG y coordinate
                    let t_to_y = move |t: f64| {
                        let t_range = temp_max.get() - temp_min.get();
                        if t_range <= 0.0 { return PROFILE_Y_BOTTOM; }
                        PROFILE_Y_BOTTOM - ((t - temp_min.get()) / t_range * PROFILE_H).clamp(0.0, PROFILE_H)
                    };

                    // Helper: map radius (cm) to SVG x coordinate
                    let r_to_x = move |r_cm: f64| {
                        PROFILE_X + (r_cm / chart_max_r_cm) * PROFILE_W
                    };

                    // Build flat polyline points for a region
                    let flat_line = move |r_start_cm: f64, r_end_cm: f64, temp: f64| {
                        let x1 = r_to_x(r_start_cm);
                        let x2 = r_to_x(r_end_cm);
                        let y = t_to_y(temp);
                        format!("{:.1},{:.1} {:.1},{:.1}", x1, y, x2, y)
                    };

                    view! {
                        <div class="sim-viz">
                            <h3>"Temperature Profile T(r)"</h3>
                            <div class="time-slider">
                                <label>"Time: " {time_display}</label>
                                <input
                                    type="range"
                                    min="0"
                                    max=move || {
                                        let dur = duration_minutes.get() as f64;
                                        let data = hot_coffee.get();
                                        data.iter()
                                            .rposition(|&(t, _)| t <= dur)
                                            .unwrap_or(data.len().saturating_sub(1))
                                    }
                                    step="1"
                                    prop:value=move || selected_time_idx.get()
                                    on:input=move |ev| {
                                        if let Ok(v) = event_target_value(&ev).parse::<usize>() {
                                            set_selected_time_idx.set(v);
                                        }
                                    }
                                />
                            </div>

                            <svg viewBox="0 0 500 240" class="profile-chart" role="img" aria-labelledby="lumped-profile-title lumped-profile-desc">
                                <title id="lumped-profile-title">"Lumped Temperature Profile"</title>
                                <desc id="lumped-profile-desc">"Flat-line temperature profile showing uniform coffee, mug, and air temperatures at the selected time"</desc>
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
                                    {move || format!("{:.0}cm", chart_max_r_cm)}
                                </text>

                                // Liquid region background
                                <rect
                                    x="50"
                                    y="20"
                                    width=move || (inner_radius_cm.get() / chart_max_r_cm * PROFILE_W) as f64
                                    height="180"
                                    fill="var(--color-cold)"
                                    opacity="0.2"
                                />

                                // Shell region background
                                <rect
                                    x=move || PROFILE_X + (inner_radius_cm.get() / chart_max_r_cm * PROFILE_W)
                                    y="20"
                                    width=move || ((outer_radius_cm.get() - inner_radius_cm.get()) / chart_max_r_cm * PROFILE_W) as f64
                                    height="180"
                                    fill="var(--color-secondary)"
                                    opacity="0.4"
                                />

                                // Inner boundary line
                                <line
                                    x1=move || PROFILE_X + (inner_radius_cm.get() / chart_max_r_cm * PROFILE_W)
                                    y1="20"
                                    x2=move || PROFILE_X + (inner_radius_cm.get() / chart_max_r_cm * PROFILE_W)
                                    y2="200"
                                    stroke="var(--color-cold)"
                                    stroke-width="1"
                                    stroke-dasharray="4 4"
                                />
                                <text
                                    x=move || PROFILE_X + (inner_radius_cm.get() / chart_max_r_cm * PROFILE_W)
                                    y="15"
                                    class="axis-label"
                                    text-anchor="middle"
                                >"r"<tspan baseline-shift="sub" font-size="8">"i"</tspan></text>

                                // Outer boundary line
                                <line
                                    x1=move || PROFILE_X + (outer_radius_cm.get() / chart_max_r_cm * PROFILE_W)
                                    y1="20"
                                    x2=move || PROFILE_X + (outer_radius_cm.get() / chart_max_r_cm * PROFILE_W)
                                    y2="200"
                                    stroke="var(--color-primary)"
                                    stroke-width="1"
                                    stroke-dasharray="4 4"
                                />
                                <text
                                    x=move || PROFILE_X + (outer_radius_cm.get() / chart_max_r_cm * PROFILE_W)
                                    y="15"
                                    class="axis-label"
                                    text-anchor="middle"
                                >"r"<tspan baseline-shift="sub" font-size="8">"o"</tspan></text>

                                // Initial profile (t=0, faded) — flat lines per region
                                <polyline
                                    class="temp-profile-line initial"
                                    points=move || flat_line(0.0, inner_radius_cm.get(), coffee_t0())
                                />
                                <polyline
                                    class="temp-profile-line initial"
                                    points=move || flat_line(inner_radius_cm.get(), outer_radius_cm.get(), mug_t0())
                                />
                                <polyline
                                    class="temp-profile-line initial"
                                    points=move || flat_line(outer_radius_cm.get(), chart_max_r_cm, room_t())
                                />

                                // Current temperature profile — flat lines per region
                                <polyline
                                    class="temp-profile-line"
                                    stroke="var(--color-cold)"
                                    points=move || flat_line(0.0, inner_radius_cm.get(), coffee_t())
                                />
                                <polyline
                                    class="temp-profile-line"
                                    stroke="var(--color-primary)"
                                    points=move || flat_line(inner_radius_cm.get(), outer_radius_cm.get(), mug_t())
                                />
                                <polyline
                                    class="temp-profile-line"
                                    stroke="var(--color-text-muted)"
                                    points=move || flat_line(outer_radius_cm.get(), chart_max_r_cm, room_t())
                                />

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
                            </svg>
                        </div>

                        <div class="sim-stats">
                            <div class="stat liquid">
                                <span class="stat-value">{move || { let v = coffee_t(); if use_fahrenheit.get() { format!("{:.1}\u{00B0}F", c_to_f(v)) } else { format!("{:.1}\u{00B0}C", v) } }}</span>
                                <span class="stat-label">"Coffee temp"</span>
                            </div>
                            <div class="stat shell">
                                <span class="stat-value">{move || { let v = mug_t(); if use_fahrenheit.get() { format!("{:.1}\u{00B0}F", c_to_f(v)) } else { format!("{:.1}\u{00B0}C", v) } }}</span>
                                <span class="stat-label">"Mug temp"</span>
                            </div>
                            <div class="stat ambient">
                                <span class="stat-value">{move || { let v = room_t(); if use_fahrenheit.get() { format!("{:.1}\u{00B0}F", c_to_f(v)) } else { format!("{:.1}\u{00B0}C", v) } }}</span>
                                <span class="stat-label">"Room temp"</span>
                            </div>
                        </div>
                    }.into()
                })
            />
        </div>
        <MobileSummaryBar
            time_cold=crossing_time_cold
            time_hot=crossing_time_hot
        />
    }
}

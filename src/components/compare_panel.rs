use crate::components::cooling_chart::{
    lowest_plotted, time_axis_bottom, to_chart_points, CHART_BOTTOM, CHART_HEIGHT, CHART_LEFT,
    CHART_RIGHT,
};
use crate::components::throttled::throttled;
use crate::components::{ControlsLayout, MugControls, PanelSettings, SimTab};
use crate::physics::{
    convert_snapshots, liquid_average_temperature, seconds_to_minutes, simulate_lumped_analytical,
    simulate_lumped_numerical, simulate_sphere_heat_numerical, simulate_three_region_relinearised,
    CoolingDataPoint, HeatMaterial, MugParameters, OuterBoundary, SphereHeatParams, ThreeRegionMug,
    ANALYTICAL_SEGMENTS, MAX_SIMULATION_MINUTES,
};
use crate::storage::{display_temp as fmt_temp, SharedSettings};
use crate::timing::SOLVE_THROTTLE_MS;
use leptos::*;

/// Everything the four solves depend on, built once from the settings and
/// grouped by solve so each solve can follow only its own inputs.
#[derive(Clone, PartialEq)]
struct CompareInputs {
    /// Both lumped solves.
    lumped: MugParameters,
    analytical: HeatEqInput<ThreeRegionMug>,
    numerical: HeatEqInput<SphereHeatParams>,
}

/// A heat-equation solve's inputs: the solver's parameters plus run length
/// and output frames.
#[derive(Clone, PartialEq)]
struct HeatEqInput<P> {
    params: P,
    duration_s: f64,
    num_snapshots: usize,
}

/// Timed result: data points in (time_minutes, temp_c) format + compute time
#[derive(Clone, PartialEq)]
struct TimedSeries {
    data: Vec<(f64, f64)>,
    elapsed_ms: f64,
}

/// The Compare panel's own controls: snapshots, series terms and domain
/// size. They render inside `MugControls` (its `layout` decides the heading
/// and the ids).
#[component]
pub fn CompareExtraControls(#[prop(optional)] layout: ControlsLayout) -> impl IntoView {
    let panel = PanelSettings::instance();
    let (num_snapshots, set_num_snapshots) = panel.num_snapshots.split();
    let (num_modes, set_num_modes) = panel.num_modes.split();
    let (r_max_mult, set_r_max_mult) = panel.r_max_mult.split();
    let prefix = layout.id_prefix();
    let id = move |name: &str| format!("{prefix}{name}");

    view! {
        {(layout == ControlsLayout::Flat).then(|| view! { <h3>"Compare Settings"</h3> })}

        <div class="control-group">
            <label for=id("cmp-snapshots")>"Snapshots: " {move || format!("{}", num_snapshots.get())}</label>
            <input
                type="range"
                id=id("cmp-snapshots")
                min="20"
                max="200"
                step="10"
                prop:value=move || num_snapshots.get()
                on:input=move |ev| {
                    if let Ok(v) = event_target_value(&ev).parse::<usize>() {
                        set_num_snapshots.set(v);
                    }
                }
            />
            <span class="control-hint">"Time resolution (all methods)"</span>
        </div>

        <div class="control-group">
            <label for=id("cmp-modes")>"Series terms: " {move || num_modes.get()}</label>
            <input
                type="range"
                id=id("cmp-modes")
                min="1"
                max="100"
                step="1"
                prop:value=move || num_modes.get()
                on:input=move |ev| {
                    if let Ok(v) = event_target_value(&ev).parse::<usize>() {
                        set_num_modes.set(v);
                    }
                }
            />
            <span class="control-hint">"Heat Eq Analytical only"</span>
        </div>

        <div class="control-group">
            <label for=id("cmp-rmax")>"r_max: " {move || format!("{:.1}x", r_max_mult.get())}</label>
            <input
                type="range"
                id=id("cmp-rmax")
                min="1.5"
                max="10"
                step="0.5"
                prop:value=move || r_max_mult.get()
                on:input=move |ev| {
                    if let Ok(v) = event_target_value(&ev).parse::<f64>() {
                        set_r_max_mult.set(v);
                    }
                }
            />
            <span class="control-hint">"Heat Equation only"</span>
        </div>
    }
}

/// `with_controls = false` leaves the controls to the caller (the overlay
/// renders them outside the panel).
#[component]
pub fn ComparePanel(#[prop(default = true)] with_controls: bool) -> impl IntoView {
    let settings = SharedSettings::instance();
    let material = settings.material;
    let volume_ml = settings.volume_ml;
    let wall_thickness_mm = settings.wall_thickness_mm;
    let coffee_temp_c = settings.coffee_temp_c;
    let room_temp_c = settings.room_temp_c;
    let min_drinkable_c = settings.min_drinkable_c;
    let use_fahrenheit = settings.use_fahrenheit;
    let duration_minutes = settings.duration_minutes;
    let preheat_temp_c = settings.preheat_temp_c;
    let h_conv = settings.h_conv;
    let emissivity = settings.emissivity;
    let exposure_pct = settings.exposure_pct;

    let panel = PanelSettings::instance();
    let num_snapshots = panel.num_snapshots;
    let num_modes = panel.num_modes;
    let r_max_mult = panel.r_max_mult;
    let (use_preheated, set_use_preheated) = panel.compare_use_preheated.split();

    // Derived geometry
    let inner_radius_cm = create_memo(move |_| {
        let v = volume_ml.get();
        (3.0 * v / (4.0 * std::f64::consts::PI)).cbrt()
    });
    let outer_radius_cm =
        create_memo(move |_| inner_radius_cm.get() + wall_thickness_mm.get() / 10.0);
    let r_max_cm = create_memo(move |_| outer_radius_cm.get() * r_max_mult.get());

    // Shell material from mug signal
    let shell_material = move || {
        let mat = material.get();
        HeatMaterial::from_mug(mat, emissivity.get())
    };

    // Shell temperature depends on preheated toggle
    let shell_temp = move || {
        if use_preheated.get() {
            preheat_temp_c.get()
        } else {
            room_temp_c.get()
        }
    };

    let compare_inputs = create_memo(move |_| {
        let duration_s = duration_minutes.get() as f64 * 60.0;
        let num_snapshots = num_snapshots.get();
        CompareInputs {
            lumped: settings.make_mug_params(use_preheated.get()),
            analytical: HeatEqInput {
                params: ThreeRegionMug {
                    r_inner: inner_radius_cm.get() / 100.0,
                    r_outer: outer_radius_cm.get() / 100.0,
                    r_max: r_max_cm.get() / 100.0,
                    shell: shell_material(),
                    liquid_temp_c: coffee_temp_c.get(),
                    shell_temp_c: shell_temp(),
                    room_temp_c: room_temp_c.get(),
                    h_conv: h_conv.get(),
                    exposure: exposure_pct.get() / 100.0,
                    num_modes: num_modes.get(),
                },
                duration_s,
                num_snapshots,
            },
            numerical: HeatEqInput {
                params: SphereHeatParams {
                    outer_radius_m: outer_radius_cm.get() / 100.0,
                    inner_radius_m: inner_radius_cm.get() / 100.0,
                    shell_temp_c: shell_temp(),
                    liquid_temp_c: coffee_temp_c.get(),
                    room_temp_c: room_temp_c.get(),
                    shell_material: shell_material(),
                    liquid_material: HeatMaterial::water(),
                    medium_material: HeatMaterial::air(),
                    exposure: exposure_pct.get() / 100.0,
                    convection_coeff: h_conv.get(),
                    r_max_multiplier: r_max_mult.get(),
                    outer_boundary: OuterBoundary::Ambient,
                },
                duration_s,
                num_snapshots,
            },
        }
    });
    // The heat-equation solves are expensive: during a drag all four follow
    // the inputs at most once per window (so the curves stay comparable), and
    // always finish on the latest inputs.
    let solve_inputs = throttled(compare_inputs, SOLVE_THROTTLE_MS);
    // One slice per solve: memo equality stops a control that only one solve
    // reads (e.g. series terms) from re-running the others.
    let lumped_input = create_memo(move |_| solve_inputs.with(|inputs| inputs.lumped.clone()));
    let analytical_input =
        create_memo(move |_| solve_inputs.with(|inputs| inputs.analytical.clone()));
    let numerical_input =
        create_memo(move |_| solve_inputs.with(|inputs| inputs.numerical.clone()));

    // 1. Lumped Analytical
    let lumped_analytical = create_memo(move |_| {
        let start = js_sys::Date::now();
        let data = lumped_input.with(|mug| simulate_lumped_analytical(mug, MAX_SIMULATION_MINUTES));
        let series = CoolingDataPoint::coffee_series(&data);
        TimedSeries {
            data: series,
            elapsed_ms: js_sys::Date::now() - start,
        }
    });

    // 2. Lumped Numerical
    let lumped_numerical = create_memo(move |_| {
        let start = js_sys::Date::now();
        let data = lumped_input.with(|mug| simulate_lumped_numerical(mug, MAX_SIMULATION_MINUTES));
        let series = CoolingDataPoint::coffee_series(&data);
        TimedSeries {
            data: series,
            elapsed_ms: js_sys::Date::now() - start,
        }
    });

    // 3. Heat Equation Analytical
    let heateq_analytical = create_memo(move |_| {
        let start = js_sys::Date::now();
        let heat_snaps = analytical_input.with(|input| {
            let display_r_max = 10.0_f64 / 100.0; // 10cm
            convert_snapshots(simulate_three_region_relinearised(
                &input.params,
                ANALYTICAL_SEGMENTS,
                input.duration_s,
                input.num_snapshots,
                80,
                display_r_max,
            ))
        });

        let avg = liquid_average_temperature(&heat_snaps);
        let series = seconds_to_minutes(&avg);
        TimedSeries {
            data: series,
            elapsed_ms: js_sys::Date::now() - start,
        }
    });

    // 4. Heat Equation Numerical
    let heateq_numerical = create_memo(move |_| {
        let start = js_sys::Date::now();
        let snaps = numerical_input.with(|input| {
            simulate_sphere_heat_numerical(&input.params, input.duration_s, input.num_snapshots)
        });
        let avg = liquid_average_temperature(&snaps);
        let series = seconds_to_minutes(&avg);
        TimedSeries {
            data: series,
            elapsed_ms: js_sys::Date::now() - start,
        }
    });

    // Temperature range for Y-axis
    let temp_top = create_memo(move |_| 100.0_f64);
    let temp_bot = create_memo(move |_| {
        let lowest = lowest_plotted(
            &[
                &lumped_analytical.get().data,
                &lumped_numerical.get().data,
                &heateq_analytical.get().data,
                &heateq_numerical.get().data,
            ],
            duration_minutes.get() as f64,
        );
        time_axis_bottom(room_temp_c.get(), lowest)
    });

    let display_temp = move |c: f64| fmt_temp(c, use_fahrenheit.get());

    // SVG y of a temperature on the chart's axis
    let temp_to_y = move |temp: f64| {
        let t_bot = temp_bot.get();
        let t_range = temp_top.get() - t_bot;
        if t_range.abs() < 1e-6 {
            100.0
        } else {
            CHART_BOTTOM - ((temp - t_bot) / t_range * CHART_HEIGHT)
        }
    };
    let threshold_y = create_memo(move |_| temp_to_y(min_drinkable_c.get()));
    let room_y = create_memo(move |_| temp_to_y(room_temp_c.get()));

    // Temperature at midpoint for stats table
    let midpoint_temp = move |data: &[(f64, f64)]| -> f64 {
        let mid_min = duration_minutes.get() as f64 / 2.0;
        data.iter()
            .min_by(|(t1, _), (t2, _)| {
                (t1 - mid_min)
                    .abs()
                    .partial_cmp(&(t2 - mid_min).abs())
                    .unwrap()
            })
            .map(|&(_, temp)| temp)
            .unwrap_or(0.0)
    };

    view! {
        <div class="model-container compare-panel">
            {with_controls.then(|| view! {
                <MugControls
                    settings=settings
                    on_reset=Callback::new(|_| SimTab::Compare.reset_own_settings())
                >
                    <CompareExtraControls/>
                </MugControls>
            })}

            <div class="temperature-display">
                // Scenario toggle
                <div class="solver-toggle">
                    <button
                        class:active=move || use_preheated.get()
                        on:click=move |_| set_use_preheated.set(true)
                    >"Preheated"</button>
                    <button
                        class:active=move || !use_preheated.get()
                        on:click=move |_| set_use_preheated.set(false)
                    >"Cold Mug"</button>
                    <span class="compute-time">
                        {move || {
                            let total = lumped_analytical.get().elapsed_ms
                                + lumped_numerical.get().elapsed_ms
                                + heateq_analytical.get().elapsed_ms
                                + heateq_numerical.get().elapsed_ms;
                            format!("{:.1} ms total", total)
                        }}
                    </span>
                </div>

                // Chart
                <div class="chart-placeholder">
                    <h4>"Coffee Temperature \u{2014} 4 Methods"</h4>
                    <p class="wip-note">
                        "Lumped (a cylinder with a base) and heat equation (a sphere) are different "
                        "idealisations of the mug, so their curves are not expected to coincide; "
                        "the theory pages list each model\u{2019}s assumptions."
                    </p>
                    <svg viewBox="0 0 400 200" class="temp-chart" role="img" aria-labelledby="compare-chart-title compare-chart-desc">
                        <title id="compare-chart-title">"Coffee Temperature Comparison"</title>
                        <desc id="compare-chart-desc">"Line chart comparing coffee temperatures across four solver methods: lumped analytical, lumped numerical, heat equation analytical, and heat equation numerical"</desc>
                        // Y-axis labels
                        <text x="5" y="20" class="axis-label">
                            {move || display_temp(temp_top.get())}
                        </text>
                        <text x="5" y=move || threshold_y.get() + 4.0 class="axis-label">
                            {move || display_temp(min_drinkable_c.get())}
                        </text>
                        <text x="5" y="168" class="axis-label">
                            {move || display_temp(temp_bot.get())}
                        </text>

                        // X-axis labels
                        <text x="40" y="195" class="axis-label">"0"</text>
                        <text x="210" y="195" class="axis-label">
                            {move || format!("{}", duration_minutes.get() / 2)}
                        </text>
                        <text x="380" y="195" class="axis-label">
                            {move || format!("{} min", duration_minutes.get())}
                        </text>

                        // Room temperature line
                        <line
                            x1=CHART_LEFT
                            y1=move || room_y.get()
                            x2=CHART_RIGHT
                            y2=move || room_y.get()
                            class="grid-line room-temp"
                        />

                        // Drinkable threshold line
                        <line
                            x1=CHART_LEFT
                            y1=move || threshold_y.get()
                            x2=CHART_RIGHT
                            y2=move || threshold_y.get()
                            class="grid-line drinkable"
                        />

                        // 1. Lumped Analytical — blue, solid
                        <polyline
                            class="temp-line lumped-analytical"
                            points=move || to_chart_points(
                                &lumped_analytical.get().data,
                                duration_minutes.get() as f64,
                                temp_top.get(),
                                temp_bot.get(),
                            )
                        />

                        // 2. Lumped Numerical — blue, dashed
                        <polyline
                            class="temp-line lumped-numerical"
                            points=move || to_chart_points(
                                &lumped_numerical.get().data,
                                duration_minutes.get() as f64,
                                temp_top.get(),
                                temp_bot.get(),
                            )
                        />

                        // 3. Heat Eq Analytical — coral, solid
                        <polyline
                            class="temp-line heateq-analytical"
                            points=move || to_chart_points(
                                &heateq_analytical.get().data,
                                duration_minutes.get() as f64,
                                temp_top.get(),
                                temp_bot.get(),
                            )
                        />

                        // 4. Heat Eq Numerical — coral, dashed
                        <polyline
                            class="temp-line heateq-numerical"
                            points=move || to_chart_points(
                                &heateq_numerical.get().data,
                                duration_minutes.get() as f64,
                                temp_top.get(),
                                temp_bot.get(),
                            )
                        />
                    </svg>

                    // Legend
                    <div class="compare-legend">
                        <div class="compare-legend-item">
                            <span class="compare-swatch lumped-analytical"></span>
                            "Lumped Analytical"
                        </div>
                        <div class="compare-legend-item">
                            <span class="compare-swatch lumped-numerical"></span>
                            "Lumped Numerical"
                        </div>
                        <div class="compare-legend-item">
                            <span class="compare-swatch heateq-analytical"></span>
                            "Heat Eq Analytical"
                        </div>
                        <div class="compare-legend-item">
                            <span class="compare-swatch heateq-numerical"></span>
                            "Heat Eq Numerical"
                        </div>
                    </div>
                </div>

                <div class="compare-stats">
                    <table>
                        <thead>
                            <tr>
                                <th>"Method"</th>
                                <th>{move || format!("Temp @ {} min", duration_minutes.get() / 2)}</th>
                                <th>"Compute"</th>
                            </tr>
                        </thead>
                        <tbody>
                            <tr class="lumped-a-row">
                                <td>"Lumped Analytical"</td>
                                <td>{move || display_temp(midpoint_temp(&lumped_analytical.get().data))}</td>
                                <td>{move || format!("{:.1} ms", lumped_analytical.get().elapsed_ms)}</td>
                            </tr>
                            <tr class="lumped-n-row">
                                <td>"Lumped Numerical"</td>
                                <td>{move || display_temp(midpoint_temp(&lumped_numerical.get().data))}</td>
                                <td>{move || format!("{:.1} ms", lumped_numerical.get().elapsed_ms)}</td>
                            </tr>
                            <tr class="heateq-a-row">
                                <td>"Heat Eq Analytical"</td>
                                <td>{move || display_temp(midpoint_temp(&heateq_analytical.get().data))}</td>
                                <td>{move || format!("{:.1} ms", heateq_analytical.get().elapsed_ms)}</td>
                            </tr>
                            <tr class="heateq-n-row">
                                <td>"Heat Eq Numerical"</td>
                                <td>{move || display_temp(midpoint_temp(&heateq_numerical.get().data))}</td>
                                <td>{move || format!("{:.1} ms", heateq_numerical.get().elapsed_ms)}</td>
                            </tr>
                        </tbody>
                    </table>
                </div>
            </div>
        </div>
    }
}

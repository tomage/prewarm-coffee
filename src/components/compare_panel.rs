use crate::components::cooling_chart::{
    to_chart_points, CHART_BOTTOM, CHART_HEIGHT, CHART_LEFT, CHART_RIGHT,
};
use crate::components::MugControls;
use crate::physics::{
    blend_exposure, compute_eigenmodes, convert_snapshots,
    liquid_average_temperature, seconds_to_minutes,
    simulate_lumped_analytical, simulate_lumped_numerical, simulate_sphere_heat_numerical,
    simulate_three_region_analytical, CoolingDataPoint, HeatMaterial,
    SphereHeatParams, ThreeRegionParams, MAX_SIMULATION_MINUTES,
};
use crate::storage::{display_temp as fmt_temp, SharedSettings};
use leptos::*;

/// Timed result: data points in (time_minutes, temp_c) format + compute time
#[derive(Clone, PartialEq)]
struct TimedSeries {
    data: Vec<(f64, f64)>,
    elapsed_ms: f64,
}

#[component]
pub fn ComparePanel() -> impl IntoView {
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

    // Panel-local signals
    let (num_snapshots, set_num_snapshots) = create_signal(100_usize);
    let (num_modes, set_num_modes) = create_signal(50_usize);
    let (r_max_mult, set_r_max_mult) = create_signal(10.0_f64);
    let (use_preheated, set_use_preheated) = create_signal(true);

    // Derived geometry
    let inner_radius_cm = create_memo(move |_| {
        let v = volume_ml.get();
        (3.0 * v / (4.0 * std::f64::consts::PI)).cbrt()
    });
    let outer_radius_cm = create_memo(move |_| {
        inner_radius_cm.get() + wall_thickness_mm.get() / 10.0
    });
    let r_max_cm = create_memo(move |_| outer_radius_cm.get() * r_max_mult.get());

    // Shell material from mug signal
    let shell_material = move || {
        let mat = material.get();
        HeatMaterial {
            conductivity: mat.thermal_conductivity(),
            density: mat.density(),
            specific_heat: mat.specific_heat(),
            emissivity: emissivity.get(),
        }
    };

    // Shell temperature depends on preheated toggle
    let shell_temp = move || {
        if use_preheated.get() {
            preheat_temp_c.get()
        } else {
            room_temp_c.get()
        }
    };

    let make_mug_params = move || settings.make_mug_params(use_preheated.get());

    // 1. Lumped Analytical
    let lumped_analytical = create_memo(move |_| {
        let start = js_sys::Date::now();
        let params = make_mug_params();
        let data = simulate_lumped_analytical(&params, MAX_SIMULATION_MINUTES);
        let series = CoolingDataPoint::coffee_series(&data);
        TimedSeries {
            data: series,
            elapsed_ms: js_sys::Date::now() - start,
        }
    });

    // 2. Lumped Numerical
    let lumped_numerical = create_memo(move |_| {
        let start = js_sys::Date::now();
        let params = make_mug_params();
        let data = simulate_lumped_numerical(&params, MAX_SIMULATION_MINUTES);
        let series = CoolingDataPoint::coffee_series(&data);
        TimedSeries {
            data: series,
            elapsed_ms: js_sys::Date::now() - start,
        }
    });

    // 3. Heat Equation Analytical
    let heateq_analytical = create_memo(move |_| {
        let start = js_sys::Date::now();
        let water = HeatMaterial::water();
        let air = HeatMaterial::air();
        let shell_mat = shell_material();
        let inner_r = inner_radius_cm.get() / 100.0;
        let outer_r = outer_radius_cm.get() / 100.0;
        let r_max = r_max_cm.get() / 100.0;
        let liq_temp = coffee_temp_c.get();
        let s_temp = shell_temp();
        let room = room_temp_c.get();
        let h = h_conv.get();
        let dur_s = duration_minutes.get() as f64 * 60.0;

        // Enhance h_conv with linearized radiation from mug outer surface
        let h_conv_eff = {
            let sigma: f64 = 5.670_374e-8;
            let t_avg_k = 0.5 * (liq_temp + room) + 273.15;
            let h_rad_outer = 4.0 * shell_mat.emissivity * sigma * t_avg_k.powi(3);
            h + h_rad_outer
        };

        let params = ThreeRegionParams {
            r_inner: inner_r,
            r_outer: outer_r,
            r_max,
            alpha_1: water.diffusivity(),
            alpha_2: shell_mat.diffusivity(),
            alpha_3: air.diffusivity(),
            k_1: water.conductivity,
            k_2: shell_mat.conductivity,
            k_3: air.conductivity,
            t_liquid_init: liq_temp,
            t_shell_init: s_temp,
            t_ambient: room,
            num_modes: num_modes.get(),
            h_conv: h_conv_eff,
        };
        let modes = compute_eigenmodes(&params);
        let display_r_max = 10.0_f64 / 100.0; // 10cm
        let snaps = simulate_three_region_analytical(
            &params,
            &modes,
            dur_s,
            num_snapshots.get(),
            80,
            display_r_max,
        );
        let heat_snaps = convert_snapshots(snaps);

        // Apply exposure blending if needed
        let exposure = exposure_pct.get() / 100.0;
        let blended = blend_exposure(
            heat_snaps, exposure, liq_temp, room, h, inner_r, 50, dur_s,
        );

        let avg = liquid_average_temperature(&blended);
        let series = seconds_to_minutes(&avg);
        TimedSeries {
            data: series,
            elapsed_ms: js_sys::Date::now() - start,
        }
    });

    // 4. Heat Equation Numerical
    let heateq_numerical = create_memo(move |_| {
        let start = js_sys::Date::now();
        let shell_mat = shell_material();
        let inner_r = inner_radius_cm.get() / 100.0;
        let outer_r = outer_radius_cm.get() / 100.0;
        let liq_temp = coffee_temp_c.get();
        let s_temp = shell_temp();
        let room = room_temp_c.get();
        let h = h_conv.get();
        let dur_s = duration_minutes.get() as f64 * 60.0;

        let params = SphereHeatParams {
            outer_radius_m: outer_r,
            inner_radius_m: inner_r,
            shell_temp_c: s_temp,
            liquid_temp_c: liq_temp,
            room_temp_c: room,
            shell_material: shell_mat,
            liquid_material: HeatMaterial::water(),
            medium_material: HeatMaterial::air(),
            evaporation_factor: exposure_pct.get() / 100.0,
            convection_coeff: h,
            r_max_multiplier: r_max_mult.get(),
        };
        let snaps =
            simulate_sphere_heat_numerical(&params, dur_s, num_snapshots.get());
        let avg = liquid_average_temperature(&snaps);
        let series = seconds_to_minutes(&avg);
        TimedSeries {
            data: series,
            elapsed_ms: js_sys::Date::now() - start,
        }
    });

    // Temperature range for Y-axis
    let temp_top = create_memo(move |_| 100.0_f64);
    let temp_bot = create_memo(move |_| room_temp_c.get());

    let display_temp = move |c: f64| fmt_temp(c, use_fahrenheit.get());

    let threshold_y = create_memo(move |_| {
        let t_top = temp_top.get();
        let t_bot = temp_bot.get();
        let t_range = t_top - t_bot;
        if t_range.abs() < 1e-6 {
            100.0
        } else {
            CHART_BOTTOM - ((min_drinkable_c.get() - t_bot) / t_range * CHART_HEIGHT)
        }
    });

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
            <MugControls
                settings=settings
                on_reset=Callback::new(move |_| {
                    set_num_snapshots.set(100);
                    set_num_modes.set(50);
                    set_r_max_mult.set(10.0);
                    set_use_preheated.set(true);
                })
            >
                <h3>"Compare Settings"</h3>

                <div class="control-group">
                    <label for="cmp-snapshots">"Snapshots: " {move || format!("{}", num_snapshots.get())}</label>
                    <input
                        type="range"
                        id="cmp-snapshots"
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
                    <label for="cmp-modes">"Series terms: " {move || num_modes.get()}</label>
                    <input
                        type="range"
                        id="cmp-modes"
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
                    <label for="cmp-rmax">"r_max: " {move || format!("{:.1}x", r_max_mult.get())}</label>
                    <input
                        type="range"
                        id="cmp-rmax"
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
            </MugControls>

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

                <p class="wip-note">
                    "Temperatures may differ substantially between models due to differing "
                    "assumptions (effective conductivity, boundary conditions, evaporation handling) "
                    "and possible errors. See the theory pages for details."
                </p>

                // Chart
                <div class="chart-placeholder">
                    <h4>"Coffee Temperature \u{2014} 4 Methods"</h4>
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
                            {move || display_temp(room_temp_c.get())}
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
                        <line x1=CHART_LEFT y1=CHART_BOTTOM x2=CHART_RIGHT y2=CHART_BOTTOM class="grid-line room-temp"/>

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

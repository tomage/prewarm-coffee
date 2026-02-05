use crate::components::{MugControls, ProfileChart, SimulationResults};
use crate::components::solver_toggle::create_solver_signals;
use crate::physics::{
    blend_exposure, compute_comparison, convert_snapshots,
    seconds_to_minutes, compute_eigenmodes,
    liquid_average_temperature, liquid_near_wall_temperature, shell_average_temperature,
    simulate_sphere_heat_numerical, simulate_three_region_analytical,
    HeatMaterial, HeatSnapshot,
    SphereHeatParams, ThreeRegionParams,
    MAX_SIMULATION_MINUTES,
};
use crate::storage::{SharedSettings, StoredSettings};
use leptos::*;
use std::cell::Cell;
use std::rc::Rc;

/// Run a heat simulation (analytical or numerical) with the given parameters.
///
/// Centralises the duplicated simulation-construction logic that differs only
/// in shell temperature, shell material, convection multiplier, and
/// evaporation factor.
fn run_heat_sim(
    use_analytical: bool,
    inner_r_cm: f64,
    outer_r_cm: f64,
    r_max_cm: f64,
    chart_max_r_cm: f64,
    liquid_temp: f64,
    shell_temp: f64,
    room_temp: f64,
    shell_mat: HeatMaterial,
    duration_s: f64,
    num_snapshots: usize,
    num_modes: usize,
    h_conv: f64,
    r_max_mult: f64,
    evaporation_factor: f64,
) -> Vec<HeatSnapshot> {
    if use_analytical {
        let water = HeatMaterial::water();
        let air = HeatMaterial::air();

        // Enhance h_conv with linearized radiation from the mug outer surface.
        //
        // The analytical eigenfunction expansion can only model heat loss at
        // the outer boundary (r = r_outer). Inner-surface effects (evaporation,
        // lid-opening convection/radiation) act at r = r_inner and are only
        // captured by the numerical solver. Remapping them to the outer boundary
        // destroys the eigenvalue structure (τ_1 collapses from ~1000s to ~150s).
        let h_conv_eff = {
            let sigma: f64 = 5.670_374e-8;
            let t_avg_k = 0.5 * (liquid_temp + room_temp) + 273.15;
            let h_rad_outer = 4.0 * shell_mat.emissivity * sigma * t_avg_k.powi(3);
            h_conv + h_rad_outer
        };

        let params = ThreeRegionParams {
            r_inner: inner_r_cm / 100.0,
            r_outer: outer_r_cm / 100.0,
            r_max: r_max_cm / 100.0,
            alpha_1: water.diffusivity(),
            alpha_2: shell_mat.diffusivity(),
            alpha_3: air.diffusivity(),
            k_1: water.conductivity,
            k_2: shell_mat.conductivity,
            k_3: air.conductivity,
            t_liquid_init: liquid_temp,
            t_shell_init: shell_temp,
            t_ambient: room_temp,
            num_modes,
            h_conv: h_conv_eff,
        };
        let modes = compute_eigenmodes(&params);
        let display_r_max = chart_max_r_cm / 100.0;
        convert_snapshots(simulate_three_region_analytical(
            &params, &modes, duration_s, num_snapshots, 80, display_r_max,
        ))
    } else {
        let params = SphereHeatParams {
            outer_radius_m: outer_r_cm / 100.0,
            inner_radius_m: inner_r_cm / 100.0,
            shell_temp_c: shell_temp,
            liquid_temp_c: liquid_temp,
            room_temp_c: room_temp,
            shell_material: shell_mat,
            liquid_material: HeatMaterial::water(),
            medium_material: HeatMaterial::air(),
            evaporation_factor,
            convection_coeff: h_conv,
            r_max_multiplier: r_max_mult,
        };
        simulate_sphere_heat_numerical(&params, duration_s, num_snapshots)
    }
}

#[component]
pub fn HeatEquationPanel() -> impl IntoView {
    // Use shared settings singleton — same RwSignals across all routes
    let settings = SharedSettings::instance();
    let material = settings.material;
    let volume_ml = settings.volume_ml;
    let shell_thickness_mm = settings.wall_thickness_mm;
    let shell_temp_c = settings.preheat_temp_c;
    let liquid_temp_c = settings.coffee_temp_c;
    let use_fahrenheit = settings.use_fahrenheit;
    let room_temp_c = settings.room_temp_c;
    let min_drinkable_c = settings.min_drinkable_c;
    let duration_minutes = settings.duration_minutes;
    let convection_coeff = settings.h_conv;
    let (exposure_pct, set_exposure_pct) = settings.exposure_pct.split();
    let (emissivity, set_emissivity) = settings.emissivity.split();

    // Panel-local signals
    let (selected_time_idx, set_selected_time_idx) = create_signal(2_usize);
    let (num_snapshots, set_num_snapshots) = create_signal(100_usize);
    let (num_modes, set_num_modes) = create_signal(50_usize);
    let (r_max_mult, set_r_max_mult) = create_signal(10.0_f64);
    let (use_near_wall, set_use_near_wall) = settings.use_near_wall.split();

    let (use_analytical, set_use_analytical, sim_ready, set_sim_ready) = create_solver_signals(settings.use_analytical_heat_eq);

    // Track compute time for display — Cells written inside memos, read by effect
    let (compute_time_ms, set_compute_time_ms) = create_signal(0.0_f64);
    let sim_time = Rc::new(Cell::new(0.0_f64));
    let sim_cold_time = Rc::new(Cell::new(0.0_f64));

    // Derived: inner radius from volume (V = 4/3 π r³, so r = (3V/4π)^(1/3))
    // volume_ml = cm³
    let inner_radius_cm = create_memo(move |_| {
        let v = volume_ml.get();
        (3.0 * v / (4.0 * std::f64::consts::PI)).cbrt()
    });

    // Derived: outer radius = inner + shell thickness
    let outer_radius_cm = create_memo(move |_| {
        inner_radius_cm.get() + shell_thickness_mm.get() / 10.0
    });

    // Chart radial extent — fixed at 10 cm for consistent display
    let chart_max_r_cm = create_memo(move |_| 10.0_f64);

    // Simulation domain extent — can be larger than chart display
    let r_max_cm = create_memo(move |_| outer_radius_cm.get() * r_max_mult.get());

    // Helper: build HeatMaterial from the current mug material signal
    let shell_material_from_mug = move || {
        let mat = material.get();
        HeatMaterial {
            conductivity: mat.thermal_conductivity(),
            density: mat.density(),
            specific_heat: mat.specific_heat(),
            emissivity: emissivity.get(),
        }
    };

    // Hot mug simulation: shell starts at preheat temperature
    // Exposure fraction drives evaporation + convection + radiation at the coffee
    // surface directly, rather than relying on post-hoc blending.
    let simulation_data = {
        let timing = sim_time.clone();
        create_memo(move |_| {
            let start = js_sys::Date::now();
            let result = run_heat_sim(
                use_analytical.get(),
                inner_radius_cm.get(), outer_radius_cm.get(),
                r_max_cm.get(), chart_max_r_cm.get(),
                liquid_temp_c.get(), shell_temp_c.get(), room_temp_c.get(),
                shell_material_from_mug(),
                duration_minutes.get() as f64 * 60.0, num_snapshots.get(), num_modes.get(),
                convection_coeff.get(), r_max_mult.get(), exposure_pct.get() / 100.0,
            );
            timing.set(js_sys::Date::now() - start);
            result
        })
    };

    // Cold mug simulation: shell starts at room temperature
    let simulation_data_cold = {
        let timing = sim_cold_time.clone();
        create_memo(move |_| {
            let start = js_sys::Date::now();
            let result = run_heat_sim(
                use_analytical.get(),
                inner_radius_cm.get(), outer_radius_cm.get(),
                r_max_cm.get(), chart_max_r_cm.get(),
                liquid_temp_c.get(), room_temp_c.get(), room_temp_c.get(),
                shell_material_from_mug(),
                duration_minutes.get() as f64 * 60.0, num_snapshots.get(), num_modes.get(),
                convection_coeff.get(), r_max_mult.get(), exposure_pct.get() / 100.0,
            );
            timing.set(js_sys::Date::now() - start);
            result
        })
    };

    // Blend liquid-region temperatures with the exposed-sphere model.
    // For analytical mode with exposure > 0%, each liquid grid point's temperature
    // is blended: T = (1 - exp) * T_3region(r,t) + exp * T_exposed_avg(t)
    let simulation_data_blended = create_memo(move |_| -> Vec<HeatSnapshot> {
        let data = simulation_data.get();
        let exposure = exposure_pct.get() / 100.0;
        if !use_analytical.get() || exposure <= 0.0 {
            return data;
        }
        blend_exposure(
            data, exposure, liquid_temp_c.get(), room_temp_c.get(),
            convection_coeff.get(), inner_radius_cm.get() / 100.0,
            50, duration_minutes.get() as f64 * 60.0,
        )
    });

    let simulation_data_cold_blended = create_memo(move |_| -> Vec<HeatSnapshot> {
        let data = simulation_data_cold.get();
        let exposure = exposure_pct.get() / 100.0;
        if !use_analytical.get() || exposure <= 0.0 {
            return data;
        }
        blend_exposure(
            data, exposure, liquid_temp_c.get(), room_temp_c.get(),
            convection_coeff.get(), inner_radius_cm.get() / 100.0,
            50, duration_minutes.get() as f64 * 60.0,
        )
    });

    // Read timing from Cells written by simulation memos (memos run before effects in Leptos 0.6)
    create_effect(move |_| {
        let _ = simulation_data_blended.get();
        let _ = simulation_data_cold_blended.get();
        let total = sim_time.get() + sim_cold_time.get();
        set_compute_time_ms.set(total);
        set_sim_ready.set(true);
    });

    let liquid_avg_temps = create_memo(move |_| {
        liquid_average_temperature(&simulation_data_blended.get())
    });

    let shell_avg_temps = create_memo(move |_| {
        shell_average_temperature(&simulation_data_blended.get())
    });

    let liquid_avg_temps_cold = create_memo(move |_| {
        liquid_average_temperature(&simulation_data_cold_blended.get())
    });

    let shell_avg_temps_cold = create_memo(move |_| {
        shell_average_temperature(&simulation_data_cold_blended.get())
    });

    // Near-wall liquid temperatures (outer 10% of liquid radius)
    let wall_fraction = 0.10;
    let liquid_wall_temps = create_memo(move |_| {
        liquid_near_wall_temperature(&simulation_data_blended.get(), wall_fraction)
    });
    let liquid_wall_temps_cold = create_memo(move |_| {
        liquid_near_wall_temperature(&simulation_data_cold_blended.get(), wall_fraction)
    });

    // Active liquid temps — switches between avg and near-wall based on toggle
    let active_liquid_temps = create_memo(move |_| {
        if use_near_wall.get() { liquid_wall_temps.get() } else { liquid_avg_temps.get() }
    });
    let active_liquid_temps_cold = create_memo(move |_| {
        if use_near_wall.get() { liquid_wall_temps_cold.get() } else { liquid_avg_temps_cold.get() }
    });

    // Convert PDE data (time_seconds, temp) to (time_minutes, temp) for charts
    let cold_coffee_chart = create_memo(move |_| seconds_to_minutes(&active_liquid_temps_cold.get()));
    let hot_coffee_chart = create_memo(move |_| seconds_to_minutes(&active_liquid_temps.get()));
    let cold_mug_chart = create_memo(move |_| seconds_to_minutes(&shell_avg_temps_cold.get()));
    let hot_mug_chart = create_memo(move |_| seconds_to_minutes(&shell_avg_temps.get()));
    let temp_top_chart = create_memo(move |_| liquid_temp_c.get());

    // Bundled comparison metrics (times in minutes for display)
    let metrics = create_memo(move |_| {
        compute_comparison(&cold_coffee_chart.get(), &hot_coffee_chart.get(), min_drinkable_c.get())
    });

    let current_snapshot = create_memo(move |_| {
        let data = simulation_data_blended.get();
        let idx = selected_time_idx.get().min(data.len() - 1);
        data[idx].clone()
    });

    let initial_snapshot = create_memo(move |_| {
        let data = simulation_data_blended.get();
        data[0].clone()
    });

    let temp_max = create_memo(move |_| {
        let hi = liquid_temp_c.get().max(shell_temp_c.get()).max(room_temp_c.get());
        (hi * 1.1 / 10.0).ceil() * 10.0
    });
    let temp_min = create_memo(move |_| {
        let lo = liquid_temp_c.get().min(shell_temp_c.get()).min(room_temp_c.get());
        (lo * 0.9 / 10.0).floor() * 10.0
    });

    view! {
        <div class="model-container">
            <MugControls
                settings=settings
                on_reset=Callback::new(move |_| {
                    let d = StoredSettings::default();
                    set_selected_time_idx.set(2);
                    set_num_snapshots.set(100);
                    set_num_modes.set(50);
                    set_r_max_mult.set(10.0);
                    set_use_near_wall.set(false);
                    set_exposure_pct.set(d.exposure_pct);
                    set_emissivity.set(d.emissivity);
                })
            >
                <h3>"Heat Equation"</h3>

                <div class="control-group">
                    <label for="pde-snapshots">"Snapshots: " {move || format!("{}", num_snapshots.get())}</label>
                    <input
                        type="range"
                        id="pde-snapshots"
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
                    <span class="control-hint">"Time resolution (output frames)"</span>
                </div>

                <div class="control-group" class:disabled=move || !use_analytical.get()>
                    <label for="pde-modes">"Series terms: " {move || num_modes.get()}</label>
                    <input
                        type="range"
                        id="pde-modes"
                        min="1"
                        max="100"
                        step="1"
                        prop:value=move || num_modes.get()
                        prop:disabled=move || !use_analytical.get()
                        on:input=move |ev| {
                            if let Ok(v) = event_target_value(&ev).parse::<usize>() {
                                set_num_modes.set(v);
                            }
                        }
                    />
                    <span class="control-hint">"Eigenfunction terms (analytical only)"</span>
                </div>

                <div class="control-group">
                    <label for="pde-rmax">"r_max: " {move || format!("{:.1}x r_outer ({:.1} cm)", r_max_mult.get(), r_max_cm.get())}</label>
                    <input
                        type="range"
                        id="pde-rmax"
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
                    <span class="control-hint">"Simulation domain extent"</span>
                </div>

            </MugControls>

            // Right side — charts + stats
            <SimulationResults
                cold_coffee=cold_coffee_chart
                hot_coffee=hot_coffee_chart
                cold_mug=cold_mug_chart
                hot_mug=hot_mug_chart
                metrics=metrics
                temp_top=temp_top_chart
                room_temp=room_temp_c.read_only()
                min_drinkable=min_drinkable_c.read_only()
                duration_minutes=duration_minutes.read_only()
                use_fahrenheit=use_fahrenheit.read_only()
                use_analytical=use_analytical
                set_use_analytical=set_use_analytical
                compute_time_ms=compute_time_ms.into()
                sim_ready=Signal::from(sim_ready)
                max_sim_minutes=MAX_SIMULATION_MINUTES as f64
                footer=Rc::new(move || view! {
                    <ProfileChart
                        current_snapshot=current_snapshot
                        initial_snapshot=initial_snapshot
                        inner_radius_cm=inner_radius_cm
                        outer_radius_cm=outer_radius_cm
                        chart_max_r_cm=chart_max_r_cm
                        temp_max=temp_max
                        temp_min=temp_min
                        convection_coeff=convection_coeff
                        room_temp_c=room_temp_c
                        use_fahrenheit=use_fahrenheit
                        simulation_data_blended=simulation_data_blended
                        selected_time_idx=selected_time_idx
                        set_selected_time_idx=set_selected_time_idx
                    />
                }.into())
            >
                <div class="solver-toggle">
                    <button
                        class:active=move || !use_near_wall.get()
                        on:click=move |_| set_use_near_wall.set(false)
                    >
                        "Liquid avg"
                    </button>
                    <button
                        class:active=move || use_near_wall.get()
                        on:click=move |_| set_use_near_wall.set(true)
                    >
                        "Near-wall"
                    </button>
                </div>
            </SimulationResults>
        </div>
    }
}

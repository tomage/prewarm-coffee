use crate::components::solver_toggle::create_solver_signals;
use crate::components::throttled::throttled;
use crate::components::{
    ControlsLayout, MugControls, PanelSettings, ProfileChart, SimTab, SimulationResults,
};
use crate::physics::{
    compute_comparison, convert_snapshots, liquid_average_temperature,
    liquid_near_wall_temperature, seconds_to_minutes, shell_average_temperature,
    simulate_sphere_heat_numerical, simulate_three_region_relinearised, HeatMaterial, HeatSnapshot,
    OuterBoundary, SphereHeatParams, ThreeRegionMug, ANALYTICAL_SEGMENTS, MAX_SIMULATION_MINUTES,
};
use crate::storage::SharedSettings;
use crate::timing::SOLVE_THROTTLE_MS;
use leptos::*;
use std::cell::Cell;
use std::rc::Rc;

/// Everything one heat simulation run depends on. The hot and cold runs
/// differ only in `shell_temp`.
#[derive(Clone, PartialEq)]
struct HeatSimInput {
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
    /// Lid openness, 0–1.
    exposure: f64,
}

/// The inputs of the panel's two runs: preheated and cold mug.
#[derive(Clone, PartialEq)]
struct HeatRunInputs {
    hot: HeatSimInput,
    cold: HeatSimInput,
}

/// Run a heat simulation (analytical or numerical) for the given input.
fn run_heat_sim(input: HeatSimInput) -> Vec<HeatSnapshot> {
    let HeatSimInput {
        use_analytical,
        inner_r_cm,
        outer_r_cm,
        r_max_cm,
        chart_max_r_cm,
        liquid_temp,
        shell_temp,
        room_temp,
        shell_mat,
        duration_s,
        num_snapshots,
        num_modes,
        h_conv,
        r_max_mult,
        exposure,
    } = input;
    if use_analytical {
        let mug = ThreeRegionMug {
            r_inner: inner_r_cm / 100.0,
            r_outer: outer_r_cm / 100.0,
            r_max: r_max_cm / 100.0,
            shell: shell_mat,
            liquid_temp_c: liquid_temp,
            shell_temp_c: shell_temp,
            room_temp_c: room_temp,
            h_conv,
            exposure,
            num_modes,
        };
        convert_snapshots(simulate_three_region_relinearised(
            &mug,
            ANALYTICAL_SEGMENTS,
            duration_s,
            num_snapshots,
            80,
            chart_max_r_cm / 100.0,
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
            exposure,
            convection_coeff: h_conv,
            r_max_multiplier: r_max_mult,
            outer_boundary: OuterBoundary::Ambient,
        };
        simulate_sphere_heat_numerical(&params, duration_s, num_snapshots)
    }
}

/// Inner radius (cm) of the spherical mug holding `volume_ml`
/// (V = 4/3 π r³, so r = (3V/4π)^(1/3); ml = cm³).
fn inner_radius_for_volume(volume_ml: f64) -> f64 {
    (3.0 * volume_ml / (4.0 * std::f64::consts::PI)).cbrt()
}

/// The heat-equation panel's own controls: snapshots, series terms and
/// domain size. They render inside `MugControls` (its `layout` decides the
/// heading and the ids).
#[component]
pub fn HeatEqExtraControls(#[prop(optional)] layout: ControlsLayout) -> impl IntoView {
    let settings = SharedSettings::instance();
    let panel = PanelSettings::instance();
    let (num_snapshots, set_num_snapshots) = panel.num_snapshots.split();
    let (num_modes, set_num_modes) = panel.num_modes.split();
    let (r_max_mult, set_r_max_mult) = panel.r_max_mult.split();
    let use_analytical = settings.use_analytical_heat_eq.read_only();
    let r_max_cm = move || {
        (inner_radius_for_volume(settings.volume_ml.get())
            + settings.wall_thickness_mm.get() / 10.0)
            * r_max_mult.get()
    };
    let prefix = layout.id_prefix();
    let id = move |name: &str| format!("{prefix}{name}");

    view! {
        {(layout == ControlsLayout::Flat).then(|| view! { <h3>"Heat Equation"</h3> })}

        <div class="control-group">
            <label for=id("pde-snapshots")>"Snapshots: " {move || format!("{}", num_snapshots.get())}</label>
            <input
                type="range"
                id=id("pde-snapshots")
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
            <label for=id("pde-modes")>"Series terms: " {move || num_modes.get()}</label>
            <input
                type="range"
                id=id("pde-modes")
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
            <label for=id("pde-rmax")>"r_max: " {move || format!("{:.1}x r_outer ({:.1} cm)", r_max_mult.get(), r_max_cm())}</label>
            <input
                type="range"
                id=id("pde-rmax")
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
    }
}

/// `with_controls = false` leaves the controls to the caller (the overlay
/// renders them outside the panel).
#[component]
pub fn HeatEquationPanel(#[prop(default = true)] with_controls: bool) -> impl IntoView {
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
    let exposure_pct = settings.exposure_pct.read_only();
    let emissivity = settings.emissivity.read_only();

    let panel = PanelSettings::instance();
    let num_snapshots = panel.num_snapshots;
    let num_modes = panel.num_modes;
    let r_max_mult = panel.r_max_mult;
    let (selected_time_idx, set_selected_time_idx) = panel.selected_time_idx.split();
    let (use_near_wall, set_use_near_wall) = settings.use_near_wall.split();

    let (use_analytical, set_use_analytical, sim_ready, set_sim_ready) =
        create_solver_signals(settings.use_analytical_heat_eq);

    // Track compute time for display — Cells written inside memos, read by effect
    let (compute_time_ms, set_compute_time_ms) = create_signal(0.0_f64);
    let sim_time = Rc::new(Cell::new(0.0_f64));
    let sim_cold_time = Rc::new(Cell::new(0.0_f64));

    let inner_radius_cm = create_memo(move |_| inner_radius_for_volume(volume_ml.get()));

    // Derived: outer radius = inner + shell thickness
    let outer_radius_cm =
        create_memo(move |_| inner_radius_cm.get() + shell_thickness_mm.get() / 10.0);

    // Chart radial extent — fixed at 10 cm for consistent display
    let chart_max_r_cm = create_memo(move |_| 10.0_f64);

    // Simulation domain extent — can be larger than chart display
    let r_max_cm = create_memo(move |_| outer_radius_cm.get() * r_max_mult.get());

    // Helper: build HeatMaterial from the current mug material signal
    let shell_material_from_mug = move || {
        let mat = material.get();
        HeatMaterial::from_mug(mat, emissivity.get())
    };

    // Simulation input from the current settings; only the shell's starting
    // temperature differs between the hot and cold runs.
    let sim_input = move |shell_temp: f64| HeatSimInput {
        use_analytical: use_analytical.get(),
        inner_r_cm: inner_radius_cm.get(),
        outer_r_cm: outer_radius_cm.get(),
        r_max_cm: r_max_cm.get(),
        chart_max_r_cm: chart_max_r_cm.get(),
        liquid_temp: liquid_temp_c.get(),
        shell_temp,
        room_temp: room_temp_c.get(),
        shell_mat: shell_material_from_mug(),
        duration_s: duration_minutes.get() as f64 * 60.0,
        num_snapshots: num_snapshots.get(),
        num_modes: num_modes.get(),
        h_conv: convection_coeff.get(),
        r_max_mult: r_max_mult.get(),
        exposure: exposure_pct.get() / 100.0,
    };

    // Hot mug starts at the preheat temperature, cold mug at room temperature.
    let run_inputs = create_memo(move |_| HeatRunInputs {
        hot: sim_input(shell_temp_c.get()),
        cold: sim_input(room_temp_c.get()),
    });
    // Both solves are expensive: during a drag they follow the inputs at most
    // once per window, and always finish on the latest inputs.
    let solve_inputs = throttled(run_inputs, SOLVE_THROTTLE_MS);
    // One slice per solve: memo equality stops a change that only affects one
    // run (e.g. the preheat temperature) from re-running the other.
    let hot_input = create_memo(move |_| solve_inputs.with(|inputs| inputs.hot.clone()));
    let cold_input = create_memo(move |_| solve_inputs.with(|inputs| inputs.cold.clone()));

    // Hot mug simulation
    let simulation_data = {
        let timing = sim_time.clone();
        create_memo(move |_| {
            let start = js_sys::Date::now();
            let result = run_heat_sim(hot_input.get());
            timing.set(js_sys::Date::now() - start);
            result
        })
    };

    // Cold mug simulation
    let simulation_data_cold = {
        let timing = sim_cold_time.clone();
        create_memo(move |_| {
            let start = js_sys::Date::now();
            let result = run_heat_sim(cold_input.get());
            timing.set(js_sys::Date::now() - start);
            result
        })
    };

    // Read timing from Cells written by simulation memos (memos run before effects in Leptos 0.6)
    create_effect(move |_| {
        let _ = simulation_data.get();
        let _ = simulation_data_cold.get();
        let total = sim_time.get() + sim_cold_time.get();
        set_compute_time_ms.set(total);
        set_sim_ready.set(true);
    });

    let liquid_avg_temps = create_memo(move |_| liquid_average_temperature(&simulation_data.get()));

    let shell_avg_temps = create_memo(move |_| shell_average_temperature(&simulation_data.get()));

    let liquid_avg_temps_cold =
        create_memo(move |_| liquid_average_temperature(&simulation_data_cold.get()));

    let shell_avg_temps_cold =
        create_memo(move |_| shell_average_temperature(&simulation_data_cold.get()));

    // Near-wall liquid temperatures (outer 10% of liquid radius)
    let wall_fraction = 0.10;
    let liquid_wall_temps =
        create_memo(move |_| liquid_near_wall_temperature(&simulation_data.get(), wall_fraction));
    let liquid_wall_temps_cold = create_memo(move |_| {
        liquid_near_wall_temperature(&simulation_data_cold.get(), wall_fraction)
    });

    // Active liquid temps — switches between avg and near-wall based on toggle
    let active_liquid_temps = create_memo(move |_| {
        if use_near_wall.get() {
            liquid_wall_temps.get()
        } else {
            liquid_avg_temps.get()
        }
    });
    let active_liquid_temps_cold = create_memo(move |_| {
        if use_near_wall.get() {
            liquid_wall_temps_cold.get()
        } else {
            liquid_avg_temps_cold.get()
        }
    });

    // Convert PDE data (time_seconds, temp) to (time_minutes, temp) for charts
    let cold_coffee_chart =
        create_memo(move |_| seconds_to_minutes(&active_liquid_temps_cold.get()));
    let hot_coffee_chart = create_memo(move |_| seconds_to_minutes(&active_liquid_temps.get()));
    let cold_mug_chart = create_memo(move |_| seconds_to_minutes(&shell_avg_temps_cold.get()));
    let hot_mug_chart = create_memo(move |_| seconds_to_minutes(&shell_avg_temps.get()));
    let temp_top_chart = create_memo(move |_| liquid_temp_c.get());

    // Bundled comparison metrics (times in minutes for display)
    let metrics = create_memo(move |_| {
        compute_comparison(
            &cold_coffee_chart.get(),
            &hot_coffee_chart.get(),
            min_drinkable_c.get(),
        )
    });

    let current_snapshot = create_memo(move |_| {
        let data = simulation_data.get();
        let idx = selected_time_idx.get().min(data.len() - 1);
        data[idx].clone()
    });

    let initial_snapshot = create_memo(move |_| {
        let data = simulation_data.get();
        data[0].clone()
    });

    let temp_max = create_memo(move |_| {
        let hi = liquid_temp_c
            .get()
            .max(shell_temp_c.get())
            .max(room_temp_c.get());
        (hi * 1.1 / 10.0).ceil() * 10.0
    });
    // Lower bound covers the lowest plotted value, which can be below room
    // once evaporation cools the coffee below room temperature.
    let temp_min = create_memo(move |_| {
        let lowest_snapshot = simulation_data.with(|snaps| {
            snaps
                .iter()
                .flat_map(|s| s.temperatures.iter().copied())
                .fold(f64::INFINITY, f64::min)
        });
        let lo = liquid_temp_c
            .get()
            .min(shell_temp_c.get())
            .min(room_temp_c.get())
            .min(lowest_snapshot);
        ((lo - 0.1 * lo.abs()) / 10.0).floor() * 10.0
    });

    view! {
        <div class="model-container">
            {with_controls.then(|| view! {
                <MugControls
                    settings=settings
                    on_reset=Callback::new(|_| SimTab::HeatEquation.reset_own_settings())
                >
                    <HeatEqExtraControls/>
                </MugControls>
            })}

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
                        snapshots=simulation_data
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
                        "Coffee avg"
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

use leptos::*;
use leptos_router::use_query_map;

/// Create panel-local solver signals backed by a shared `RwSignal<bool>`.
///
/// The `?solver=numerical` query param overrides the stored value on first load.
/// Returns `(use_analytical, set_use_analytical, sim_ready, set_sim_ready)`.
pub fn create_solver_signals(shared: RwSignal<bool>) -> (ReadSignal<bool>, WriteSignal<bool>, ReadSignal<bool>, WriteSignal<bool>) {
    let query = use_query_map();
    let has_numerical_override = query.with_untracked(|q| {
        q.get("solver").map(|v| v == "numerical").unwrap_or(false)
    });
    if has_numerical_override {
        shared.set(false);
    }
    let (use_analytical, set_use_analytical) = shared.split();
    let (sim_ready, set_sim_ready) = create_signal(false);
    (use_analytical, set_use_analytical, sim_ready, set_sim_ready)
}

#[component]
pub fn SolverToggle(
    use_analytical: ReadSignal<bool>,
    set_use_analytical: WriteSignal<bool>,
    compute_time_ms: Signal<f64>,
    /// Whether the simulation has completed at least once (for e2e test hooks)
    sim_ready: Signal<bool>,
) -> impl IntoView {
    view! {
        <div class="solver-toggle" attr:data-sim-ready=move || sim_ready.get()>
            <button
                class:active=move || use_analytical.get()
                on:click=move |_| set_use_analytical.set(true)
            >
                "Analytical"
            </button>
            <button
                class:active=move || !use_analytical.get()
                on:click=move |_| set_use_analytical.set(false)
            >
                "Numerical"
            </button>
            <span class="compute-time">
                {move || format!("{:.1} ms", compute_time_ms.get())}
            </span>
        </div>
    }
}

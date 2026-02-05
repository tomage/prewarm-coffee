use crate::components::{ComparisonMetrics, CoolingChart, SolverToggle};
use crate::physics::ComparisonData;
use leptos::*;

/// Reusable wrapper around SolverToggle + ComparisonMetrics + CoolingChart.
///
/// Both LumpedPanel and HeatEquationPanel share this same layout; only the
/// data computation differs.  An optional `children` slot renders between the
/// SolverToggle and ComparisonMetrics (used by HeatEquationPanel for the
/// near-wall toggle).  An optional `footer` slot renders after the CoolingChart
/// (used by HeatEquationPanel for the temperature profile + sim stats).
#[component]
pub fn SimulationResults(
    // Chart data (all in time_minutes)
    cold_coffee: Memo<Vec<(f64, f64)>>,
    hot_coffee: Memo<Vec<(f64, f64)>>,
    cold_mug: Memo<Vec<(f64, f64)>>,
    hot_mug: Memo<Vec<(f64, f64)>>,
    // Metrics (pre-computed, times in minutes)
    metrics: Memo<ComparisonData>,
    // Display settings
    temp_top: Memo<f64>,
    room_temp: ReadSignal<f64>,
    min_drinkable: ReadSignal<f64>,
    duration_minutes: ReadSignal<usize>,
    use_fahrenheit: ReadSignal<bool>,
    // Solver toggle
    use_analytical: ReadSignal<bool>,
    set_use_analytical: WriteSignal<bool>,
    compute_time_ms: Signal<f64>,
    sim_ready: Signal<bool>,
    // Optional
    #[prop(default = f64::MAX)]
    max_sim_minutes: f64,
    #[prop(optional)]
    children: Option<Children>,
    #[prop(optional)]
    footer: Option<ChildrenFn>,
) -> impl IntoView {
    // Derive crossing times from metrics for CoolingChart
    let crossing_cold = create_memo(move |_| metrics.get().crossing_time_cold);
    let crossing_hot = create_memo(move |_| metrics.get().crossing_time_hot);

    view! {
        <div class="temperature-display">
            <SolverToggle
                use_analytical=use_analytical
                set_use_analytical=set_use_analytical
                compute_time_ms=compute_time_ms
                sim_ready=sim_ready
            />
            {children.map(|c| c())}
            <ComparisonMetrics
                metrics=metrics
                min_drinkable=min_drinkable
                duration_minutes=duration_minutes
                use_fahrenheit=use_fahrenheit
                max_sim_minutes=max_sim_minutes
            />
            <CoolingChart
                cold_coffee=cold_coffee
                hot_coffee=hot_coffee
                cold_mug=cold_mug
                hot_mug=hot_mug
                temp_top=temp_top
                room_temp=room_temp
                min_drinkable=min_drinkable
                duration_minutes=duration_minutes
                use_fahrenheit=use_fahrenheit
                crossing_cold=crossing_cold
                crossing_hot=crossing_hot
            />
            {footer.as_ref().map(|f| f())}
        </div>
    }
}

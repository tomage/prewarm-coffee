use crate::physics::ComparisonData;
use crate::storage::{display_temp as fmt_temp, format_minutes};
use leptos::*;

#[component]
pub fn ComparisonMetrics(
    metrics: Memo<ComparisonData>,
    min_drinkable: ReadSignal<f64>,
    duration_minutes: ReadSignal<usize>,
    use_fahrenheit: ReadSignal<bool>,
    /// Maximum simulation time in minutes. Shows "X+" when time is at this cap.
    #[prop(default = f64::MAX)]
    max_sim_minutes: f64,
) -> impl IntoView {
    let display_temp = move |c: f64| fmt_temp(c, use_fahrenheit.get());

    let display_time = move |minutes: f64| -> String {
        if minutes >= max_sim_minutes - 0.5 {
            format!("{:.0}+ min", max_sim_minutes)
        } else {
            format_minutes(minutes)
        }
    };

    view! {
        <div class="metrics">
            <div class="metric">
                <span class="metric-value cold">{move || display_temp(metrics.get().avg_temp_cold)}</span>
                <span class="metric-label">"Cold mug avg"</span>
            </div>
            <div class="metric">
                <span class="metric-value hot">{move || display_temp(metrics.get().avg_temp_hot)}</span>
                <span class="metric-label">"Preheated avg"</span>
            </div>
            <div class="metric highlight">
                <span class="metric-value">{move || {
                    let diff = metrics.get().temp_diff;
                    if use_fahrenheit.get() { format!("+{:.0}\u{00B0}F", diff * 9.0 / 5.0) } else { format!("+{:.0}\u{00B0}C", diff) }
                }}</span>
                <span class="metric-label">"Warmer"</span>
            </div>
        </div>
        <p class="metrics-note">
            {move || format!("(Average for first {}, until cold mug drops below {})",
                display_time(metrics.get().crossing_time_cold),
                display_temp(min_drinkable.get())
            )}
        </p>

        <div class="drinkable-time">
            <h4>{move || format!("Time above {} (drinkable)", display_temp(min_drinkable.get()))}</h4>
            <div class="time-comparison">
                <div class="time-bar cold">
                    <div
                        class="time-fill"
                        style=move || format!("width: {}%", (metrics.get().crossing_time_cold / duration_minutes.get() as f64 * 100.0).min(100.0))
                    ></div>
                    <span>{move || display_time(metrics.get().crossing_time_cold)}</span>
                </div>
                <div class="time-bar hot">
                    <div
                        class="time-fill"
                        style=move || format!("width: {}%", (metrics.get().crossing_time_hot / duration_minutes.get() as f64 * 100.0).min(100.0))
                    ></div>
                    <span>{move || display_time(metrics.get().crossing_time_hot)}</span>
                </div>
            </div>
            <div class="time-labels">
                <span class="cold-label">"Cold mug"</span>
                <span class="hot-label">"Preheated"</span>
                <span class="time-diff">{move || {
                    let m = metrics.get();
                    if m.crossing_time_hot >= max_sim_minutes - 0.5 {
                        return String::new();
                    }
                    let diff = m.crossing_time_hot - m.crossing_time_cold;
                    if diff > 0.5 { format!("+{}", format_minutes(diff)) } else { String::new() }
                }}</span>
            </div>
        </div>
    }
}

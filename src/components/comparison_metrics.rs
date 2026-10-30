use crate::physics::{extra_time, ComparisonData, ExtraTime};
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

    // A missing crossing is a real answer: the coffee never got undrinkable.
    let display_crossing = move |crossing: Option<f64>| -> String {
        match crossing {
            Some(minutes) => display_time(minutes),
            None => format!("stays above for the whole {} min", duration_minutes.get()),
        }
    };

    // Bars show time spent drinkable; no crossing means the whole run.
    let bar_width = move |crossing: Option<f64>| -> String {
        let fraction = match crossing {
            Some(minutes) => minutes / duration_minutes.get() as f64,
            None => 1.0,
        };
        format!("width: {}%", (fraction * 100.0).min(100.0))
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
            {move || match metrics.get().crossing_time_cold {
                Some(minutes) => format!("(Average for first {}, until cold mug drops below {})",
                    display_time(minutes),
                    display_temp(min_drinkable.get())
                ),
                None => format!("(Average for the whole {} min; cold mug never drops below {})",
                    duration_minutes.get(),
                    display_temp(min_drinkable.get())
                ),
            }}
        </p>

        <div class="drinkable-time">
            <h4>{move || format!("Time above {} (drinkable)", display_temp(min_drinkable.get()))}</h4>
            <div class="time-comparison">
                <div class="time-bar cold">
                    <div
                        class="time-fill"
                        style=move || bar_width(metrics.get().crossing_time_cold)
                    ></div>
                    <span>{move || display_crossing(metrics.get().crossing_time_cold)}</span>
                </div>
                <div class="time-bar hot">
                    <div
                        class="time-fill"
                        style=move || bar_width(metrics.get().crossing_time_hot)
                    ></div>
                    <span>{move || display_crossing(metrics.get().crossing_time_hot)}</span>
                </div>
            </div>
            <div class="time-labels">
                <span class="cold-label">"Cold mug"</span>
                <span class="hot-label">"Preheated"</span>
                <span class="time-diff">{move || {
                    let m = metrics.get();
                    let duration = duration_minutes.get() as f64;
                    match extra_time(m.crossing_time_cold, m.crossing_time_hot, duration) {
                        Some(ExtraTime::Exact(diff)) if diff > 0.5 => {
                            let pct = match m.crossing_time_cold {
                                Some(cold) if cold > 0.5 => format!(" (+{:.0}%)", diff / cold * 100.0),
                                _ => String::new(),
                            };
                            format!("+{}{}", format_minutes(diff), pct)
                        }
                        Some(ExtraTime::AtLeast(diff)) if diff > 0.5 => {
                            format!("at least +{}", format_minutes(diff))
                        }
                        _ => String::new(),
                    }
                }}</span>
            </div>
        </div>
    }
}

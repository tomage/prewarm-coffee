use crate::physics::{extra_time, ExtraTime};
use crate::storage::format_minutes;
use leptos::*;

#[component]
pub fn MobileSummaryBar(
    time_cold: Memo<Option<f64>>,
    time_hot: Memo<Option<f64>>,
    duration_minutes: ReadSignal<usize>,
) -> impl IntoView {
    // No crossing means the coffee stayed drinkable for the whole run.
    let display = move |crossing: Option<f64>| match crossing {
        Some(minutes) => format_minutes(minutes),
        None => format!("all {} min", duration_minutes.get()),
    };

    view! {
        <div class="mobile-summary">
            <div class="summary-content">
                <div class="summary-item">
                    <span class="summary-value cold">{move || display(time_cold.get())}</span>
                    <span class="summary-label">"Cold"</span>
                </div>
                <div class="summary-item">
                    <span class="summary-value hot">{move || display(time_hot.get())}</span>
                    <span class="summary-label">"Preheated"</span>
                </div>
                <div class="summary-item">
                    <span class="summary-value diff">{move || {
                        let duration = duration_minutes.get() as f64;
                        match extra_time(time_cold.get(), time_hot.get(), duration) {
                            Some(ExtraTime::Exact(d)) if d > 0.5 => format!("+{}", format_minutes(d)),
                            Some(ExtraTime::AtLeast(d)) if d > 0.5 => {
                                format!("at least +{}", format_minutes(d))
                            }
                            _ => "—".to_string(),
                        }
                    }}</span>
                    <span class="summary-label">"Extra time"</span>
                </div>
            </div>
        </div>
    }
}

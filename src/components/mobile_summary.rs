use crate::storage::format_minutes;
use leptos::*;

#[component]
pub fn MobileSummaryBar(
    time_cold: Memo<f64>,
    time_hot: Memo<f64>,
) -> impl IntoView {
    let time_diff = create_memo(move |_| time_hot.get() - time_cold.get());

    view! {
        <div class="mobile-summary">
            <div class="summary-content">
                <div class="summary-item">
                    <span class="summary-value cold">{move || format_minutes(time_cold.get())}</span>
                    <span class="summary-label">"Cold"</span>
                </div>
                <div class="summary-item">
                    <span class="summary-value hot">{move || format_minutes(time_hot.get())}</span>
                    <span class="summary-label">"Preheated"</span>
                </div>
                <div class="summary-item">
                    <span class="summary-value diff">{move || {
                        let diff = time_diff.get();
                        if diff > 0.5 { format!("+{}", format_minutes(diff)) } else { "—".to_string() }
                    }}</span>
                    <span class="summary-label">"Extra time"</span>
                </div>
            </div>
        </div>
    }
}

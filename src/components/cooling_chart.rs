use crate::storage::{display_temp as fmt_temp, format_minutes};
use leptos::*;

// Chart layout constants (SVG coordinate space within viewBox="0 0 400 200")
pub(crate) const CHART_LEFT: f64 = 40.0;
pub(crate) const CHART_WIDTH: f64 = 350.0;
pub(crate) const CHART_RIGHT: f64 = CHART_LEFT + CHART_WIDTH; // 390
pub(crate) const CHART_BOTTOM: f64 = 165.0;
pub(crate) const CHART_HEIGHT: f64 = 145.0;

/// Convert time/temp data points to SVG polyline coordinates.
pub(crate) fn to_chart_points(
    data: &[(f64, f64)],
    duration: f64,
    t_top: f64,
    t_bot: f64,
) -> String {
    let t_range = t_top - t_bot;
    if t_range.abs() < 1e-6 {
        return String::new();
    }
    data.iter()
        .filter(|(t, _)| *t <= duration)
        .map(|(time, temp)| {
            let x = CHART_LEFT + (time / duration) * CHART_WIDTH;
            let y = CHART_BOTTOM - ((temp - t_bot) / t_range * CHART_HEIGHT);
            format!("{:.1},{:.1}", x, y)
        })
        .collect::<Vec<_>>()
        .join(" ")
}

/// Lowest temperature among the points plotted up to `duration` minutes,
/// or `+∞` when nothing is plotted.
pub(crate) fn lowest_plotted(series: &[&[(f64, f64)]], duration: f64) -> f64 {
    series
        .iter()
        .flat_map(|s| s.iter())
        .filter(|(t, _)| *t <= duration)
        .map(|&(_, temp)| temp)
        .fold(f64::INFINITY, f64::min)
}

/// Bottom of a temperature-over-time y-axis: room temperature, or the whole
/// degree below the lowest plotted value when coffee settles below room.
pub(crate) fn time_axis_bottom(room: f64, lowest: f64) -> f64 {
    room.min(lowest.floor())
}

/// Shared Temperature Over Time chart used by both lumped and PDE simulations.
///
/// All time-series data must be in (time_minutes, temp_celsius) format.
#[component]
pub fn CoolingChart(
    /// Cold mug coffee temperature: (time_minutes, temp_c)
    cold_coffee: Memo<Vec<(f64, f64)>>,
    /// Preheated mug coffee temperature: (time_minutes, temp_c)
    hot_coffee: Memo<Vec<(f64, f64)>>,
    /// Cold mug shell temperature: (time_minutes, temp_c)
    cold_mug: Memo<Vec<(f64, f64)>>,
    /// Preheated mug shell temperature: (time_minutes, temp_c)
    hot_mug: Memo<Vec<(f64, f64)>>,
    /// Temperature at top of Y-axis (e.g. 100.0 for lumped, coffee_temp for PDE)
    temp_top: Memo<f64>,
    /// Room temperature (drawn as a reference line)
    room_temp: ReadSignal<f64>,
    /// Minimum drinkable temperature (threshold line)
    min_drinkable: ReadSignal<f64>,
    /// Chart X-axis duration in minutes
    duration_minutes: ReadSignal<usize>,
    /// Display temperatures in Fahrenheit
    use_fahrenheit: ReadSignal<bool>,
    /// Time in minutes when cold mug crosses threshold
    crossing_cold: Memo<Option<f64>>,
    /// Time in minutes when preheated mug crosses threshold
    crossing_hot: Memo<Option<f64>>,
) -> impl IntoView {
    let display_temp = move |c: f64| fmt_temp(c, use_fahrenheit.get());

    let temp_bot = create_memo(move |_| {
        let lowest = lowest_plotted(
            &[
                &cold_coffee.get(),
                &hot_coffee.get(),
                &cold_mug.get(),
                &hot_mug.get(),
            ],
            duration_minutes.get() as f64,
        );
        time_axis_bottom(room_temp.get(), lowest)
    });

    // SVG y of a temperature on this chart's axis
    let temp_to_y = move |temp: f64| {
        let t_bot = temp_bot.get();
        let t_range = temp_top.get() - t_bot;
        if t_range.abs() < 1e-6 {
            100.0
        } else {
            CHART_BOTTOM - ((temp - t_bot) / t_range * CHART_HEIGHT)
        }
    };
    let threshold_y = create_memo(move |_| temp_to_y(min_drinkable.get()));
    let room_y = create_memo(move |_| temp_to_y(room_temp.get()));

    view! {
        <div class="chart-placeholder">
            <h4>"Temperature Over Time"</h4>
            <svg viewBox="0 0 400 200" class="temp-chart" role="img" aria-labelledby="cooling-chart-title cooling-chart-desc">
                <title id="cooling-chart-title">"Temperature Over Time"</title>
                <desc id="cooling-chart-desc">"Line chart showing coffee and mug temperatures over time for preheated and cold mugs"</desc>
                // Y-axis labels
                <text x="5" y="20" class="axis-label">
                    {move || display_temp(temp_top.get())}
                </text>
                <text x="5" y=move || threshold_y.get() + 4.0 class="axis-label">
                    {move || display_temp(min_drinkable.get())}
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

                // Cold mug — coffee
                <polyline
                    class="temp-line cold"
                    points=move || to_chart_points(
                        &cold_coffee.get(),
                        duration_minutes.get() as f64,
                        temp_top.get(),
                        temp_bot.get(),
                    )
                />

                // Preheated — coffee
                <polyline
                    class="temp-line hot"
                    points=move || to_chart_points(
                        &hot_coffee.get(),
                        duration_minutes.get() as f64,
                        temp_top.get(),
                        temp_bot.get(),
                    )
                />

                // Cold mug — shell (mug temp)
                <polyline
                    class="temp-line cold mug-temp"
                    points=move || to_chart_points(
                        &cold_mug.get(),
                        duration_minutes.get() as f64,
                        temp_top.get(),
                        temp_bot.get(),
                    )
                />

                // Preheated — shell (mug temp)
                <polyline
                    class="temp-line hot mug-temp"
                    points=move || to_chart_points(
                        &hot_mug.get(),
                        duration_minutes.get() as f64,
                        temp_top.get(),
                        temp_bot.get(),
                    )
                />

                // Crossing marker — cold mug
                {move || {
                    let duration = duration_minutes.get() as f64;
                    // No crossing: nothing to mark, the line stays above the threshold.
                    let mins = crossing_cold.get().filter(|m| *m < duration);
                    if let Some(mins) = mins {
                        let x = CHART_LEFT + (mins / duration) * CHART_WIDTH;
                        let y = threshold_y.get();
                        Some(view! {
                            <circle cx=x cy=y r="5" class="crossing-marker cold"/>
                            <text x=x y=y-10.0 class="crossing-label cold">{format_minutes(mins)}</text>
                        })
                    } else {
                        None
                    }
                }}

                // Crossing marker — preheated
                {move || {
                    let duration = duration_minutes.get() as f64;
                    // No crossing: nothing to mark, the line stays above the threshold.
                    let mins = crossing_hot.get().filter(|m| *m < duration);
                    if let Some(mins) = mins {
                        let x = CHART_LEFT + (mins / duration) * CHART_WIDTH;
                        let y = threshold_y.get();
                        Some(view! {
                            <circle cx=x cy=y r="5" class="crossing-marker hot"/>
                            <text x=x y=y-10.0 class="crossing-label hot">{format_minutes(mins)}</text>
                        })
                    } else {
                        None
                    }
                }}
            </svg>
            <div class="chart-legend">
                <div class="legend-row">
                    <span class="legend-label">"Coffee"</span>
                    <span class="legend-item cold">"Cold mug"</span>
                    <span class="legend-item hot">"Preheated"</span>
                </div>
                <div class="legend-row">
                    <span class="legend-label">"Mug"</span>
                    <span class="legend-item cold dotted">"Cold mug"</span>
                    <span class="legend-item hot dotted">"Preheated"</span>
                </div>
            </div>
        </div>
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn time_axis_bottom_covers_values_below_room() {
        let cold = [(0.0, 30.0), (60.0, 19.4), (130.0, 15.0)];
        let hot = [(0.0, 30.0), (60.0, 20.2)];
        let lowest = lowest_plotted(&[&cold, &hot], 120.0);
        assert_eq!(lowest, 19.4);
        assert_eq!(time_axis_bottom(22.0, lowest), 19.0);
        // Coffee that stays above room keeps room as the bottom.
        assert_eq!(time_axis_bottom(22.0, 40.0), 22.0);
        assert_eq!(lowest_plotted(&[], 120.0), f64::INFINITY);
    }
}

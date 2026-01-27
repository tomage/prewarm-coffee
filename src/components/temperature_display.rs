use crate::physics::CoolingDataPoint;
use crate::storage::c_to_f;
use leptos::*;

// Must match MAX_SIMULATION_MINUTES in app.rs
const MAX_SIMULATION_MINUTES: f64 = 120.0;

#[component]
pub fn TemperatureDisplay(
    cold_mug_data: Memo<Vec<CoolingDataPoint>>,
    hot_mug_data: Memo<Vec<CoolingDataPoint>>,
    room_temp: ReadSignal<f64>,
    min_drinkable: ReadSignal<f64>,
    use_fahrenheit: ReadSignal<bool>,
    duration_minutes: ReadSignal<usize>,
    set_duration_minutes: WriteSignal<usize>,
) -> impl IntoView {
    // Helper to display temperature in current unit
    let display_temp = move |c: f64| -> String {
        if use_fahrenheit.get() {
            format!("{:.0}°F", c_to_f(c))
        } else {
            format!("{:.0}°C", c)
        }
    };

    // Helper to display time, with "+" if at max simulation duration
    let display_time = |minutes: f64| -> String {
        if minutes >= MAX_SIMULATION_MINUTES - 0.5 {
            format!("{:.0}+", MAX_SIMULATION_MINUTES)
        } else {
            format!("{:.0}", minutes)
        }
    };
    // Find when coffee crosses threshold (interpolate between data points)
    let crossing_time_cold = create_memo(move |_| {
        let threshold = min_drinkable.get();
        let data = cold_mug_data.get();
        for i in 1..data.len() {
            let prev = &data[i - 1];
            let curr = &data[i];
            if prev.coffee_temp_c >= threshold && curr.coffee_temp_c < threshold {
                // Linear interpolation
                let t =
                    (threshold - prev.coffee_temp_c) / (curr.coffee_temp_c - prev.coffee_temp_c);
                return Some(prev.time_minutes + t);
            }
        }
        // Never crossed - stays above threshold for entire duration
        data.last().map(|p| p.time_minutes)
    });

    let crossing_time_hot = create_memo(move |_| {
        let threshold = min_drinkable.get();
        let data = hot_mug_data.get();
        for i in 1..data.len() {
            let prev = &data[i - 1];
            let curr = &data[i];
            if prev.coffee_temp_c >= threshold && curr.coffee_temp_c < threshold {
                // Linear interpolation
                let t =
                    (threshold - prev.coffee_temp_c) / (curr.coffee_temp_c - prev.coffee_temp_c);
                return Some(prev.time_minutes + t);
            }
        }
        // Never crossed - stays above threshold for entire duration
        data.last().map(|p| p.time_minutes)
    });

    // Time above threshold in minutes (using crossing time)
    let time_above_threshold_cold = create_memo(move |_| crossing_time_cold.get().unwrap_or(0.0));

    let time_above_threshold_hot = create_memo(move |_| crossing_time_hot.get().unwrap_or(0.0));

    // Calculate average temperature from 0 until cold mug hits threshold
    let avg_temp_cold = create_memo(move |_| {
        let threshold = min_drinkable.get();
        let data = cold_mug_data.get();
        let drinkable_points: Vec<_> = data
            .iter()
            .take_while(|p| p.coffee_temp_c >= threshold)
            .collect();
        if drinkable_points.is_empty() {
            return 0.0;
        }
        let sum: f64 = drinkable_points.iter().map(|p| p.coffee_temp_c).sum();
        sum / drinkable_points.len() as f64
    });

    let avg_temp_hot = create_memo(move |_| {
        let threshold = min_drinkable.get();
        let cold_data = cold_mug_data.get();
        let hot_data = hot_mug_data.get();
        // Use same time window as cold mug (until cold mug hits threshold)
        let cold_drinkable_count = cold_data
            .iter()
            .take_while(|p| p.coffee_temp_c >= threshold)
            .count();
        if cold_drinkable_count == 0 {
            return 0.0;
        }
        let sum: f64 = hot_data
            .iter()
            .take(cold_drinkable_count)
            .map(|p| p.coffee_temp_c)
            .sum();
        sum / cold_drinkable_count as f64
    });

    let temp_diff = create_memo(move |_| avg_temp_hot.get() - avg_temp_cold.get());

    view! {
        <div class="temperature-display">
            <div class="metrics">
                <div class="metric">
                    <span class="metric-value cold">{move || {
                        let c = avg_temp_cold.get();
                        if use_fahrenheit.get() { format!("{:.0}°F", c_to_f(c)) } else { format!("{:.0}°C", c) }
                    }}</span>
                    <span class="metric-label">"Cold mug avg"</span>
                </div>
                <div class="metric">
                    <span class="metric-value hot">{move || {
                        let c = avg_temp_hot.get();
                        if use_fahrenheit.get() { format!("{:.0}°F", c_to_f(c)) } else { format!("{:.0}°C", c) }
                    }}</span>
                    <span class="metric-label">"Preheated avg"</span>
                </div>
                <div class="metric highlight">
                    <span class="metric-value">{move || {
                        let diff = temp_diff.get();
                        if use_fahrenheit.get() { format!("+{:.0}°F", diff * 9.0 / 5.0) } else { format!("+{:.0}°C", diff) }
                    }}</span>
                    <span class="metric-label">"Warmer"</span>
                </div>
            </div>
            <p class="metrics-note">
                {move || format!("(Average for first {} min, until cold mug drops below {})",
                    display_time(time_above_threshold_cold.get()),
                    display_temp(min_drinkable.get())
                )}
            </p>

            <div class="drinkable-time">
                <h4>{move || format!("Time above {} (drinkable)", display_temp(min_drinkable.get()))}</h4>
                <div class="time-comparison">
                    <div class="time-bar cold">
                        <div
                            class="time-fill"
                            style=move || format!("width: {}%", (time_above_threshold_cold.get() / duration_minutes.get() as f64 * 100.0).min(100.0))
                        ></div>
                        <span>{move || format!("{} min", display_time(time_above_threshold_cold.get()))}</span>
                    </div>
                    <div class="time-bar hot">
                        <div
                            class="time-fill"
                            style=move || format!("width: {}%", (time_above_threshold_hot.get() / duration_minutes.get() as f64 * 100.0).min(100.0))
                        ></div>
                        <span>{move || format!("{} min", display_time(time_above_threshold_hot.get()))}</span>
                    </div>
                </div>
                <div class="time-labels">
                    <span class="cold-label">"Cold mug"</span>
                    <span class="hot-label">"Preheated"</span>
                    <span class="time-diff">{move || {
                        let hot_time = time_above_threshold_hot.get();
                        // Don't show diff if hot mug is at max (we don't know actual time)
                        if hot_time >= MAX_SIMULATION_MINUTES - 0.5 {
                            return String::new();
                        }
                        let diff = hot_time - time_above_threshold_cold.get();
                        if diff > 0.5 { format!("+{:.0} min", diff) } else { String::new() }
                    }}</span>
                </div>
            </div>

            <div class="chart-placeholder">
                <h4>"Temperature Over Time"</h4>
                <svg viewBox="0 0 400 200" class="temp-chart">
                    // Y-axis labels
                    <text x="5" y="20" class="axis-label">{move || if use_fahrenheit.get() { "212°F" } else { "100°C" }}</text>
                    <text x="5" y=move || 165.0 - ((min_drinkable.get() - room_temp.get()) / (100.0 - room_temp.get()) * 145.0) + 4.0 class="axis-label">{move || display_temp(min_drinkable.get())}</text>
                    <text x="5" y="168" class="axis-label">{move || display_temp(room_temp.get())}</text>

                    // X-axis labels
                    <text x="40" y="195" class="axis-label">"0"</text>
                    <text x="210" y="195" class="axis-label">{move || format!("{}", duration_minutes.get() / 2)}</text>
                    <text x="380" y="195" class="axis-label">{move || format!("{} min", duration_minutes.get())}</text>

                    // Room temperature line (baseline)
                    <line
                        x1="40"
                        y1="165"
                        x2="390"
                        y2="165"
                        class="grid-line room-temp"
                    />

                    // Drinkable threshold line (calculated based on scale)
                    <line
                        x1="40"
                        y1=move || 165.0 - ((min_drinkable.get() - room_temp.get()) / (100.0 - room_temp.get()) * 145.0)
                        x2="390"
                        y2=move || 165.0 - ((min_drinkable.get() - room_temp.get()) / (100.0 - room_temp.get()) * 145.0)
                        class="grid-line drinkable"
                    />

                    // Cold mug line
                    <polyline
                        class="temp-line cold"
                        points=move || {
                            let duration = duration_minutes.get() as f64;
                            cold_mug_data.get()
                                .iter()
                                .filter(|p| p.time_minutes <= duration)
                                .map(|p| {
                                    let x = 40.0 + (p.time_minutes / duration) * 350.0;
                                    let y = 165.0 - ((p.coffee_temp_c - room_temp.get()) / (100.0 - room_temp.get()) * 145.0);
                                    format!("{:.1},{:.1}", x, y)
                                })
                                .collect::<Vec<_>>()
                                .join(" ")
                        }
                    />

                    // Hot mug line (coffee)
                    <polyline
                        class="temp-line hot"
                        points=move || {
                            let duration = duration_minutes.get() as f64;
                            hot_mug_data.get()
                                .iter()
                                .filter(|p| p.time_minutes <= duration)
                                .map(|p| {
                                    let x = 40.0 + (p.time_minutes / duration) * 350.0;
                                    let y = 165.0 - ((p.coffee_temp_c - room_temp.get()) / (100.0 - room_temp.get()) * 145.0);
                                    format!("{:.1},{:.1}", x, y)
                                })
                                .collect::<Vec<_>>()
                                .join(" ")
                        }
                    />

                    // Cold mug temperature line (mug itself)
                    <polyline
                        class="temp-line cold mug-temp"
                        points=move || {
                            let duration = duration_minutes.get() as f64;
                            cold_mug_data.get()
                                .iter()
                                .filter(|p| p.time_minutes <= duration)
                                .map(|p| {
                                    let x = 40.0 + (p.time_minutes / duration) * 350.0;
                                    let y = 165.0 - ((p.mug_temp_c - room_temp.get()) / (100.0 - room_temp.get()) * 145.0);
                                    format!("{:.1},{:.1}", x, y)
                                })
                                .collect::<Vec<_>>()
                                .join(" ")
                        }
                    />

                    // Hot mug temperature line (mug itself)
                    <polyline
                        class="temp-line hot mug-temp"
                        points=move || {
                            let duration = duration_minutes.get() as f64;
                            hot_mug_data.get()
                                .iter()
                                .filter(|p| p.time_minutes <= duration)
                                .map(|p| {
                                    let x = 40.0 + (p.time_minutes / duration) * 350.0;
                                    let y = 165.0 - ((p.mug_temp_c - room_temp.get()) / (100.0 - room_temp.get()) * 145.0);
                                    format!("{:.1},{:.1}", x, y)
                                })
                                .collect::<Vec<_>>()
                                .join(" ")
                        }
                    />

                    // Marker for cold mug crossing threshold
                    {move || crossing_time_cold.get().map(|t| {
                        let x = 40.0 + (t / duration_minutes.get() as f64) * 350.0;
                        let y = 165.0 - ((min_drinkable.get() - room_temp.get()) / (100.0 - room_temp.get()) * 145.0);
                        view! {
                            <circle cx=x cy=y r="5" class="crossing-marker cold"/>
                            <text x=x y=y-10.0 class="crossing-label cold">{format!("{:.0}m", t)}</text>
                        }
                    })}

                    // Marker for hot mug crossing threshold
                    {move || crossing_time_hot.get().map(|t| {
                        let x = 40.0 + (t / duration_minutes.get() as f64) * 350.0;
                        let y = 165.0 - ((min_drinkable.get() - room_temp.get()) / (100.0 - room_temp.get()) * 145.0);
                        view! {
                            <circle cx=x cy=y r="5" class="crossing-marker hot"/>
                            <text x=x y=y-10.0 class="crossing-label hot">{format!("{:.0}m", t)}</text>
                        }
                    })}
                </svg>
                <div class="chart-legend">
                    <div class="legend-row">
                        <span class="legend-label">"Coffee"</span>
                        <span class="legend-item cold">"When roomtemp"</span>
                        <span class="legend-item hot">"When preheated"</span>
                    </div>
                    <div class="legend-row">
                        <span class="legend-label">"Mug"</span>
                        <span class="legend-item cold dotted">"When roomtemp"</span>
                        <span class="legend-item hot dotted">"When preheated"</span>
                    </div>
                </div>
                <div class="chart-duration">
                    <label for="duration">"Duration: " {move || format!("{} min", duration_minutes.get())}</label>
                    <input
                        type="range"
                        id="duration"
                        min="10"
                        max="120"
                        step="5"
                        prop:value=move || duration_minutes.get()
                        on:input=move |ev| {
                            if let Ok(v) = event_target_value(&ev).parse::<usize>() {
                                set_duration_minutes.set(v);
                            }
                        }
                    />
                </div>
            </div>
        </div>
    }
}

use crate::components::{MobileSummaryBar, MugControls, TemperatureDisplay};
use crate::physics::{simulate_cooling, CoolingDataPoint, MugParameters};
use crate::storage::{load_settings, material_to_string, save_settings, StoredSettings};
use leptos::*;

/// Find when coffee crosses below a threshold temperature (interpolated)
fn find_crossing_time(data: &[CoolingDataPoint], threshold: f64) -> Option<f64> {
    for i in 1..data.len() {
        let prev = &data[i - 1];
        let curr = &data[i];
        if prev.coffee_temp_c >= threshold && curr.coffee_temp_c < threshold {
            // Linear interpolation
            let t = (threshold - prev.coffee_temp_c) / (curr.coffee_temp_c - prev.coffee_temp_c);
            return Some(prev.time_minutes + t);
        }
    }
    // Never crossed - stays above threshold for entire duration
    data.last().map(|p| p.time_minutes)
}

#[component]
pub fn App() -> impl IntoView {
    // Load saved settings
    let saved = load_settings();

    // Mug parameters as signals (initialized from storage)
    let (material, set_material) = create_signal(saved.material_enum());
    let (volume_ml, set_volume_ml) = create_signal(saved.volume_ml);
    let (wall_thickness_mm, set_wall_thickness_mm) = create_signal(saved.wall_thickness_mm);
    let (coffee_temp_c, set_coffee_temp_c) = create_signal(saved.coffee_temp_c);
    let (room_temp_c, set_room_temp_c) = create_signal(saved.room_temp_c);
    let (min_drinkable_c, set_min_drinkable_c) = create_signal(saved.min_drinkable_c);
    let (use_fahrenheit, set_use_fahrenheit) = create_signal(saved.use_fahrenheit);
    let (duration_minutes, set_duration_minutes) = create_signal(saved.duration_minutes);
    let (preheat_temp_c, set_preheat_temp_c) = create_signal(saved.preheat_temp_c);

    // Save settings whenever any value changes
    create_effect(move |_| {
        let settings = StoredSettings {
            material: material_to_string(material.get()),
            volume_ml: volume_ml.get(),
            wall_thickness_mm: wall_thickness_mm.get(),
            coffee_temp_c: coffee_temp_c.get(),
            room_temp_c: room_temp_c.get(),
            min_drinkable_c: min_drinkable_c.get(),
            use_fahrenheit: use_fahrenheit.get(),
            duration_minutes: duration_minutes.get(),
            preheat_temp_c: preheat_temp_c.get(),
        };
        save_settings(&settings);
    });

    // Compute temperature curves for comparison
    // Always simulate for max duration (120 min) so metrics aren't clamped by chart duration
    const MAX_SIMULATION_MINUTES: usize = 120;

    let cooling_data_cold = create_memo(move |_| {
        let params = MugParameters {
            material: material.get(),
            volume_ml: volume_ml.get(),
            wall_thickness_mm: wall_thickness_mm.get(),
            coffee_temp_c: coffee_temp_c.get(),
            room_temp_c: room_temp_c.get(),
            preheated: false,
            preheat_temp_c: preheat_temp_c.get(),
        };
        simulate_cooling(&params, MAX_SIMULATION_MINUTES)
    });

    let cooling_data_hot = create_memo(move |_| {
        let params = MugParameters {
            material: material.get(),
            volume_ml: volume_ml.get(),
            wall_thickness_mm: wall_thickness_mm.get(),
            coffee_temp_c: coffee_temp_c.get(),
            room_temp_c: room_temp_c.get(),
            preheated: true,
            preheat_temp_c: preheat_temp_c.get(),
        };
        simulate_cooling(&params, MAX_SIMULATION_MINUTES)
    });

    // Calculate crossing times for mobile summary bar
    let crossing_time_cold = create_memo(move |_| {
        find_crossing_time(&cooling_data_cold.get(), min_drinkable_c.get()).unwrap_or(0.0)
    });

    let crossing_time_hot = create_memo(move |_| {
        find_crossing_time(&cooling_data_hot.get(), min_drinkable_c.get()).unwrap_or(0.0)
    });

    view! {
        <main class="container">
            <header>
                <h1>"Prewarm"</h1>
                <p class="tagline">"Why preheating your mug keeps coffee warmer"</p>
            </header>

            <section class="intro">
                <h2>"The Science"</h2>
                <p>
                    "When you pour hot coffee into a cold mug, the mug absorbs heat from your coffee. "
                    "A ceramic mug at room temperature can pull several degrees from your coffee "
                    "in the first few minutes as it warms up."
                </p>
                <p>
                    "By preheating your mug with hot water first, the mug starts at a higher temperature "
                    "and steals less heat from your coffee. The result? Your coffee stays warmer, longer—"
                    "and you get more time in the ideal drinking range."
                </p>
            </section>

            <section class="interactive">
                <h2>"Try It Yourself"</h2>
                <div class="model-container">
                    <MugControls
                        material=material
                        set_material=set_material
                        volume_ml=volume_ml
                        set_volume_ml=set_volume_ml
                        wall_thickness_mm=wall_thickness_mm
                        set_wall_thickness_mm=set_wall_thickness_mm
                        coffee_temp_c=coffee_temp_c
                        set_coffee_temp_c=set_coffee_temp_c
                        preheat_temp_c=preheat_temp_c
                        set_preheat_temp_c=set_preheat_temp_c
                        room_temp_c=room_temp_c
                        set_room_temp_c=set_room_temp_c
                        min_drinkable_c=min_drinkable_c
                        set_min_drinkable_c=set_min_drinkable_c
                        use_fahrenheit=use_fahrenheit
                        set_use_fahrenheit=set_use_fahrenheit
                        set_duration_minutes=set_duration_minutes
                    />
                    <TemperatureDisplay
                        cold_mug_data=cooling_data_cold
                        hot_mug_data=cooling_data_hot
                        room_temp=room_temp_c
                        min_drinkable=min_drinkable_c
                        use_fahrenheit=use_fahrenheit
                        duration_minutes=duration_minutes
                        set_duration_minutes=set_duration_minutes
                    />
                </div>
            </section>

            <section class="explanation">
                <h2>"How It Works"</h2>
                <div class="physics-cards">
                    <div class="card">
                        <h3>"Thermal Mass"</h3>
                        <p>
                            "Your mug has thermal mass—it takes energy to change its temperature. "
                            "A thick ceramic mug has more thermal mass than a thin glass one, "
                            "meaning it absorbs more heat from your coffee."
                        </p>
                    </div>
                    <div class="card">
                        <h3>"Heat Transfer"</h3>
                        <p>
                            "Heat flows from hot to cold. The bigger the temperature difference "
                            "between your coffee and the mug, the faster heat transfers. "
                            "A preheated mug has a smaller temperature gap."
                        </p>
                    </div>
                    <div class="card">
                        <h3>"Newton's Law of Cooling"</h3>
                        <p>
                            "The rate of heat loss is proportional to the temperature difference "
                            "with the environment. As your coffee approaches room temperature, "
                            "it cools more slowly."
                        </p>
                    </div>
                </div>
            </section>

            <footer>
                <p class="disclaimer">
                    "This simulation uses simplified physics models and approximations. "
                    "Real-world results may vary based on factors not modeled here (evaporation, lid usage, mug shape, etc.). "
                    <a href="https://github.com/tomage/prewarm-coffee/issues" target="_blank" rel="noopener">"Report inaccuracies"</a>
                </p>
                <p class="analytics-note">
                    "Analytics by "
                    <a href="https://umami.is" target="_blank" rel="noopener">"Umami"</a>
                    " · No cookies · No personal data"
                </p>
                <p>
                    <a href="https://github.com/tomage/prewarm-coffee" target="_blank" rel="noopener">"Source code"</a>
                    " · AGPL-3.0 · "
                    "© 2025"
                </p>
            </footer>

            <MobileSummaryBar
                time_cold=crossing_time_cold
                time_hot=crossing_time_hot
            />
        </main>
    }
}

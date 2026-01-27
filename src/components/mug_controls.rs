use crate::analytics::track_event;
use crate::haptics::vibrate_tick;
use crate::physics::MugMaterial;
use crate::storage::{c_to_f, f_to_c};
use leptos::*;
use std::cell::Cell;

thread_local! {
    static INTERACTION_TRACKED: Cell<bool> = const { Cell::new(false) };
}

fn track_first_interaction() {
    INTERACTION_TRACKED.with(|tracked| {
        if !tracked.get() {
            tracked.set(true);
            track_event("controls-used");
        }
    });
}

#[component]
pub fn MugControls(
    material: ReadSignal<MugMaterial>,
    set_material: WriteSignal<MugMaterial>,
    volume_ml: ReadSignal<f64>,
    set_volume_ml: WriteSignal<f64>,
    wall_thickness_mm: ReadSignal<f64>,
    set_wall_thickness_mm: WriteSignal<f64>,
    coffee_temp_c: ReadSignal<f64>,
    set_coffee_temp_c: WriteSignal<f64>,
    preheat_temp_c: ReadSignal<f64>,
    set_preheat_temp_c: WriteSignal<f64>,
    room_temp_c: ReadSignal<f64>,
    set_room_temp_c: WriteSignal<f64>,
    min_drinkable_c: ReadSignal<f64>,
    set_min_drinkable_c: WriteSignal<f64>,
    use_fahrenheit: ReadSignal<bool>,
    set_use_fahrenheit: WriteSignal<bool>,
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

    view! {
        <div class="controls">
            <h3>"Mug Properties"</h3>

            <div class="control-group">
                <label for="material">"Material"</label>
                <select
                    id="material"
                    on:change=move |ev| {
                        let value = event_target_value(&ev);
                        let mat = match value.as_str() {
                            "ceramic" => MugMaterial::Ceramic,
                            "glass" => MugMaterial::Glass,
                            "stainless" => MugMaterial::StainlessSteel,
                            "plastic" => MugMaterial::Plastic,
                            _ => MugMaterial::Ceramic,
                        };
                        set_material.set(mat);
                        track_first_interaction();
                    }
                >
                    <option value="ceramic" selected=move || matches!(material.get(), MugMaterial::Ceramic)>
                        "Ceramic"
                    </option>
                    <option value="glass" selected=move || matches!(material.get(), MugMaterial::Glass)>
                        "Glass"
                    </option>
                    <option value="stainless" selected=move || matches!(material.get(), MugMaterial::StainlessSteel)>
                        "Stainless Steel"
                    </option>
                    <option value="plastic" selected=move || matches!(material.get(), MugMaterial::Plastic)>
                        "Plastic"
                    </option>
                </select>
            </div>

            <div class="control-group">
                <label for="volume">"Volume: " {move || format!("{:.0} ml", volume_ml.get())}</label>
                <input
                    type="range"
                    id="volume"
                    min="40"
                    max="500"
                    step="10"
                    prop:value=move || volume_ml.get()
                    on:input=move |ev| {
                        if let Ok(v) = event_target_value(&ev).parse::<f64>() {
                            set_volume_ml.set(v);
                            vibrate_tick();
                            track_first_interaction();
                        }
                    }
                />
            </div>

            <div class="control-group">
                <label for="thickness">"Wall Thickness: " {move || format!("{:.1} mm", wall_thickness_mm.get())}</label>
                <input
                    type="range"
                    id="thickness"
                    min="1"
                    max="10"
                    step="0.5"
                    prop:value=move || wall_thickness_mm.get()
                    on:input=move |ev| {
                        if let Ok(v) = event_target_value(&ev).parse::<f64>() {
                            set_wall_thickness_mm.set(v);
                            vibrate_tick();
                        }
                    }
                />
            </div>

            <div class="control-group-header">
                <h3>"Temperatures"</h3>
                <button
                    class="unit-toggle"
                    on:click=move |_| {
                        set_use_fahrenheit.set(!use_fahrenheit.get());
                        track_first_interaction();
                    }
                >
                    {move || if use_fahrenheit.get() { "°F → °C" } else { "°C → °F" }}
                </button>
            </div>

            <div class="control-group">
                <label for="coffee-temp">"Coffee: " {move || display_temp(coffee_temp_c.get())}</label>
                <input
                    type="range"
                    id="coffee-temp"
                    min=move || if use_fahrenheit.get() { "122" } else { "50" }
                    max=move || if use_fahrenheit.get() { "212" } else { "100" }
                    step="1"
                    prop:value=move || {
                        if use_fahrenheit.get() {
                            c_to_f(coffee_temp_c.get()).round()
                        } else {
                            coffee_temp_c.get().round()
                        }
                    }
                    on:input=move |ev| {
                        if let Ok(v) = event_target_value(&ev).parse::<f64>() {
                            let c = if use_fahrenheit.get() { f_to_c(v).round() } else { v };
                            set_coffee_temp_c.set(c);
                            vibrate_tick();
                        }
                    }
                />
            </div>

            <div class="control-group">
                <label for="preheat-temp">"Preheat mug: " {move || display_temp(preheat_temp_c.get())}</label>
                <input
                    type="range"
                    id="preheat-temp"
                    min=move || if use_fahrenheit.get() { "122" } else { "50" }
                    max=move || if use_fahrenheit.get() { "212" } else { "100" }
                    step="1"
                    prop:value=move || {
                        if use_fahrenheit.get() {
                            c_to_f(preheat_temp_c.get()).round()
                        } else {
                            preheat_temp_c.get().round()
                        }
                    }
                    on:input=move |ev| {
                        if let Ok(v) = event_target_value(&ev).parse::<f64>() {
                            let c = if use_fahrenheit.get() { f_to_c(v).round() } else { v };
                            set_preheat_temp_c.set(c);
                            vibrate_tick();
                        }
                    }
                />
            </div>

            <div class="control-group">
                <label for="room-temp">"Room: " {move || display_temp(room_temp_c.get())}</label>
                <input
                    type="range"
                    id="room-temp"
                    min=move || if use_fahrenheit.get() { "14" } else { "-10" }
                    max=move || if use_fahrenheit.get() { "122" } else { "50" }
                    step="1"
                    prop:value=move || {
                        if use_fahrenheit.get() {
                            c_to_f(room_temp_c.get()).round()
                        } else {
                            room_temp_c.get().round()
                        }
                    }
                    on:input=move |ev| {
                        if let Ok(v) = event_target_value(&ev).parse::<f64>() {
                            let c = if use_fahrenheit.get() { f_to_c(v).round() } else { v };
                            set_room_temp_c.set(c);
                            vibrate_tick();
                        }
                    }
                />
            </div>

            <div class="control-group">
                <label for="min-drinkable">"Too cold below: " {move || display_temp(min_drinkable_c.get())}</label>
                <input
                    type="range"
                    id="min-drinkable"
                    min=move || if use_fahrenheit.get() { "122" } else { "50" }
                    max=move || if use_fahrenheit.get() { "158" } else { "70" }
                    step="1"
                    prop:value=move || {
                        if use_fahrenheit.get() {
                            c_to_f(min_drinkable_c.get()).round()
                        } else {
                            min_drinkable_c.get().round()
                        }
                    }
                    on:input=move |ev| {
                        if let Ok(v) = event_target_value(&ev).parse::<f64>() {
                            let c = if use_fahrenheit.get() { f_to_c(v).round() } else { v };
                            set_min_drinkable_c.set(c);
                            vibrate_tick();
                        }
                    }
                />
            </div>

            <button
                class="reset-button"
                on:click=move |_| {
                    set_material.set(MugMaterial::Ceramic);
                    set_volume_ml.set(180.0);
                    set_wall_thickness_mm.set(5.0);
                    set_coffee_temp_c.set(90.0);
                    set_preheat_temp_c.set(90.0);
                    set_room_temp_c.set(22.0);
                    set_min_drinkable_c.set(60.0);
                    set_duration_minutes.set(60);
                }
            >
                "Reset to defaults"
            </button>
        </div>
    }
}

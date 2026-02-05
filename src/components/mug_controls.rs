use crate::analytics::track_event;
use crate::haptics::vibrate_tick;
use crate::physics::MugMaterial;
use crate::storage::{c_to_f, display_temp as fmt_temp, f_to_c, SharedSettings, StoredSettings};
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
    settings: SharedSettings,
    #[prop(optional)] children: Option<Children>,
    #[prop(optional, into)] on_reset: Option<Callback<()>>,
) -> impl IntoView {
    let (material, set_material) = settings.material.split();
    let (volume_ml, set_volume_ml) = settings.volume_ml.split();
    let (wall_thickness_mm, set_wall_thickness_mm) = settings.wall_thickness_mm.split();
    let (coffee_temp_c, set_coffee_temp_c) = settings.coffee_temp_c.split();
    let (preheat_temp_c, set_preheat_temp_c) = settings.preheat_temp_c.split();
    let (room_temp_c, set_room_temp_c) = settings.room_temp_c.split();
    let (min_drinkable_c, set_min_drinkable_c) = settings.min_drinkable_c.split();
    let (use_fahrenheit, set_use_fahrenheit) = settings.use_fahrenheit.split();
    let (h_conv, set_h_conv) = settings.h_conv.split();
    let (emissivity, set_emissivity) = settings.emissivity.split();
    let (exposure_pct, set_exposure_pct) = settings.exposure_pct.split();
    let (duration_minutes, set_duration_minutes) = settings.duration_minutes.split();

    let display_temp = move |c: f64| fmt_temp(c, use_fahrenheit.get());

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
                        set_emissivity.set(mat.emissivity());
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

            <div class="control-group">
                <label for="h-conv">"Convection: " {move || format!("{:.0} W/(m\u{00B2}\u{00B7}K)", h_conv.get())}</label>
                <input
                    type="range"
                    id="h-conv"
                    min="2"
                    max="50"
                    step="1"
                    prop:value=move || h_conv.get()
                    on:input=move |ev| {
                        if let Ok(v) = event_target_value(&ev).parse::<f64>() {
                            set_h_conv.set(v);
                            vibrate_tick();
                        }
                    }
                />
            </div>

            <div class="control-group">
                <label for="emissivity">"Emissivity: " {move || format!("{:.2}", emissivity.get())}</label>
                <input
                    type="range"
                    id="emissivity"
                    min="0"
                    max="1"
                    step="0.01"
                    prop:value=move || emissivity.get()
                    on:input=move |ev| {
                        if let Ok(v) = event_target_value(&ev).parse::<f64>() {
                            set_emissivity.set(v);
                            vibrate_tick();
                        }
                    }
                />
                <span class="control-hint">"0 = no radiation, 1 = perfect blackbody"</span>
            </div>

            <div class="control-group">
                <label for="exposure">"Lid opening: " {move || format!("{:.0}%", exposure_pct.get())}</label>
                <input
                    type="range"
                    id="exposure"
                    min="0"
                    max="100"
                    step="5"
                    prop:value=move || exposure_pct.get()
                    on:input=move |ev| {
                        if let Ok(v) = event_target_value(&ev).parse::<f64>() {
                            set_exposure_pct.set(v);
                            vibrate_tick();
                        }
                    }
                />
                <span class="control-hint">"0% = lid on, 100% = fully open"</span>
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
                    min=move || if use_fahrenheit.get() { "32" } else { "0" }
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
                    min=move || if use_fahrenheit.get() { "14" } else { "-10" }
                    max=move || if use_fahrenheit.get() { "230" } else { "110" }
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

            <h3>"Simulation"</h3>

            <div class="control-group">
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
                            vibrate_tick();
                            track_first_interaction();
                        }
                    }
                />
            </div>

            {children.map(|c| c())}

            <button
                class="reset-button"
                on:click=move |_| {
                    let d = StoredSettings::default();
                    set_material.set(d.material_enum());
                    set_volume_ml.set(d.volume_ml);
                    set_wall_thickness_mm.set(d.wall_thickness_mm);
                    set_coffee_temp_c.set(d.coffee_temp_c);
                    set_preheat_temp_c.set(d.preheat_temp_c);
                    set_room_temp_c.set(d.room_temp_c);
                    set_min_drinkable_c.set(d.min_drinkable_c);
                    set_h_conv.set(d.h_conv);
                    set_emissivity.set(d.emissivity);
                    set_exposure_pct.set(d.exposure_pct);
                    set_duration_minutes.set(d.duration_minutes);
                    // Preserve the user's F/C preference across resets
                    if let Some(cb) = on_reset {
                        cb.call(());
                    }
                }
            >
                "Reset to defaults"
            </button>
        </div>
    }
}

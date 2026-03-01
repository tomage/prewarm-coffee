use crate::components::AccordionGroup;
use crate::haptics::vibrate_tick;
use crate::physics::MugMaterial;
use crate::storage::{c_to_f, display_temp as fmt_temp, f_to_c, SharedSettings, StoredSettings};
use leptos::*;

/// Lock/unlock body scroll by toggling a CSS class
fn set_body_scroll_lock(locked: bool) {
    if let Some(body) = document().body() {
        let class_list = body.class_list();
        if locked {
            let _ = class_list.add_1("sim-overlay-open");
        } else {
            let _ = class_list.remove_1("sim-overlay-open");
        }
    }
}

#[component]
pub fn SimulatorOverlay(
    show: ReadSignal<bool>,
    set_show: WriteSignal<bool>,
    active_tab: RwSignal<String>,
    chart: ChildrenFn,
) -> impl IntoView {
    let settings = SharedSettings::instance();

    // Lock/unlock body scroll when overlay opens/closes
    create_effect(move |_| {
        set_body_scroll_lock(show.get());
    });

    // Cleanup on unmount
    on_cleanup(move || set_body_scroll_lock(false));

    let close = move |_| set_show.set(false);

    // Destructure settings for sliders
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
        <Show when=move || show.get()>
            <div class="simulator-overlay">
                <div class="overlay-header">
                    <div class="overlay-tabs">
                        <button
                            class:active=move || active_tab.get() == "lumped"
                            on:click=move |_| active_tab.set("lumped".to_string())
                        >"Lumped"</button>
                        <button
                            class:active=move || active_tab.get() == "heat_equation"
                            on:click=move |_| active_tab.set("heat_equation".to_string())
                        >"Heat Eq"</button>
                        <button
                            class:active=move || active_tab.get() == "compare"
                            on:click=move |_| active_tab.set("compare".to_string())
                        >"Compare"</button>
                    </div>
                    <button class="close-btn" on:click=close>"Done"</button>
                </div>

                <div class="overlay-chart">
                    {chart()}
                </div>

                <div class="overlay-controls">
                    <AccordionGroup title="Mug Properties">
                        <div class="control-group">
                            <label for="m-material">"Material"</label>
                            <select
                                id="m-material"
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
                                }
                            >
                                <option value="ceramic" selected=move || matches!(material.get(), MugMaterial::Ceramic)>"Ceramic"</option>
                                <option value="glass" selected=move || matches!(material.get(), MugMaterial::Glass)>"Glass"</option>
                                <option value="stainless" selected=move || matches!(material.get(), MugMaterial::StainlessSteel)>"Stainless Steel"</option>
                                <option value="plastic" selected=move || matches!(material.get(), MugMaterial::Plastic)>"Plastic"</option>
                            </select>
                        </div>
                        <div class="control-group">
                            <label for="m-volume">"Volume: " {move || format!("{:.0} ml", volume_ml.get())}</label>
                            <input type="range" id="m-volume" min="40" max="500" step="10"
                                prop:value=move || volume_ml.get()
                                on:input=move |ev| {
                                    if let Ok(v) = event_target_value(&ev).parse::<f64>() {
                                        set_volume_ml.set(v);
                                        vibrate_tick();
                                    }
                                }
                            />
                        </div>
                        <div class="control-group">
                            <label for="m-thickness">"Wall Thickness: " {move || format!("{:.1} mm", wall_thickness_mm.get())}</label>
                            <input type="range" id="m-thickness" min="1" max="10" step="0.5"
                                prop:value=move || wall_thickness_mm.get()
                                on:input=move |ev| {
                                    if let Ok(v) = event_target_value(&ev).parse::<f64>() {
                                        set_wall_thickness_mm.set(v);
                                        vibrate_tick();
                                    }
                                }
                            />
                        </div>
                    </AccordionGroup>

                    <AccordionGroup title="Temperatures" initially_open=true>
                        <div class="control-group">
                            <label for="m-coffee-temp">"Coffee: " {move || display_temp(coffee_temp_c.get())}</label>
                            <input type="range" id="m-coffee-temp"
                                min=move || if use_fahrenheit.get() { "32" } else { "0" }
                                max=move || if use_fahrenheit.get() { "212" } else { "100" }
                                step="1"
                                prop:value=move || if use_fahrenheit.get() { c_to_f(coffee_temp_c.get()).round() } else { coffee_temp_c.get().round() }
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
                            <label for="m-preheat-temp">"Preheat mug: " {move || display_temp(preheat_temp_c.get())}</label>
                            <input type="range" id="m-preheat-temp"
                                min=move || if use_fahrenheit.get() { "14" } else { "-10" }
                                max=move || if use_fahrenheit.get() { "230" } else { "110" }
                                step="1"
                                prop:value=move || if use_fahrenheit.get() { c_to_f(preheat_temp_c.get()).round() } else { preheat_temp_c.get().round() }
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
                            <label for="m-room-temp">"Room: " {move || display_temp(room_temp_c.get())}</label>
                            <input type="range" id="m-room-temp"
                                min=move || if use_fahrenheit.get() { "14" } else { "-10" }
                                max=move || if use_fahrenheit.get() { "122" } else { "50" }
                                step="1"
                                prop:value=move || if use_fahrenheit.get() { c_to_f(room_temp_c.get()).round() } else { room_temp_c.get().round() }
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
                            <label for="m-min-drinkable">"Too cold below: " {move || display_temp(min_drinkable_c.get())}</label>
                            <input type="range" id="m-min-drinkable"
                                min=move || if use_fahrenheit.get() { "122" } else { "50" }
                                max=move || if use_fahrenheit.get() { "158" } else { "70" }
                                step="1"
                                prop:value=move || if use_fahrenheit.get() { c_to_f(min_drinkable_c.get()).round() } else { min_drinkable_c.get().round() }
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
                            class="unit-toggle"
                            on:click=move |_| set_use_fahrenheit.set(!use_fahrenheit.get())
                        >
                            {move || if use_fahrenheit.get() { "\u{00B0}F \u{2192} \u{00B0}C" } else { "\u{00B0}C \u{2192} \u{00B0}F" }}
                        </button>
                    </AccordionGroup>

                    <AccordionGroup title="Simulation">
                        <div class="control-group">
                            <label for="m-duration">"Duration: " {move || format!("{} min", duration_minutes.get())}</label>
                            <input type="range" id="m-duration" min="10" max="120" step="5"
                                prop:value=move || duration_minutes.get()
                                on:input=move |ev| {
                                    if let Ok(v) = event_target_value(&ev).parse::<usize>() {
                                        set_duration_minutes.set(v);
                                        vibrate_tick();
                                    }
                                }
                            />
                        </div>
                    </AccordionGroup>

                    <AccordionGroup title="Advanced">
                        <div class="control-group">
                            <label for="m-h-conv">"Convection: " {move || format!("{:.0} W/(m\u{00B2}\u{00B7}K)", h_conv.get())}</label>
                            <input type="range" id="m-h-conv" min="2" max="50" step="1"
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
                            <label for="m-emissivity">"Emissivity: " {move || format!("{:.2}", emissivity.get())}</label>
                            <input type="range" id="m-emissivity" min="0" max="1" step="0.01"
                                prop:value=move || emissivity.get()
                                on:input=move |ev| {
                                    if let Ok(v) = event_target_value(&ev).parse::<f64>() {
                                        set_emissivity.set(v);
                                        vibrate_tick();
                                    }
                                }
                            />
                        </div>
                        <div class="control-group">
                            <label for="m-exposure">"Lid opening: " {move || format!("{:.0}%", exposure_pct.get())}</label>
                            <input type="range" id="m-exposure" min="0" max="100" step="5"
                                prop:value=move || exposure_pct.get()
                                on:input=move |ev| {
                                    if let Ok(v) = event_target_value(&ev).parse::<f64>() {
                                        set_exposure_pct.set(v);
                                        vibrate_tick();
                                    }
                                }
                            />
                        </div>
                    </AccordionGroup>

                    <div style="padding: 0.75rem 1rem;">
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
                            }
                        >"Reset to defaults"</button>
                    </div>
                </div>
            </div>
        </Show>
    }
}

use crate::analytics::track_event;
use crate::components::AccordionGroup;
use crate::haptics::vibrate_tick;
use crate::physics::MugMaterial;
use crate::storage::{
    c_to_f, display_temp as fmt_temp, f_to_c, material_to_string, SharedSettings,
};
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

/// How the controls are arranged.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ControlsLayout {
    /// One card with section headings, beside the chart (desktop, theory pages).
    #[default]
    Flat,
    /// Collapsible groups under the chart (mobile overlay).
    Accordion,
}

impl ControlsLayout {
    /// Prefix for element ids, so ids name their layout. Only one layout is
    /// mounted at a time; the overlay's ids keep their `m-` prefix so they
    /// stay stable for tests and anything else that names them.
    pub fn id_prefix(self) -> &'static str {
        match self {
            ControlsLayout::Flat => "",
            ControlsLayout::Accordion => "m-",
        }
    }
}

/// The mug, temperature and simulation controls plus Reset.
///
/// `children` are the active panel's own controls: after Duration in the
/// flat layout, in the Advanced group in the accordion layout. `on_reset`
/// runs in the same batch as the shared reset.
#[component]
pub fn MugControls(
    settings: SharedSettings,
    #[prop(optional)] layout: ControlsLayout,
    #[prop(optional)] children: Option<Children>,
    #[prop(optional, into)] on_reset: Option<Callback<()>>,
) -> impl IntoView {
    let material = settings.material.read_only();
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
    let prefix = layout.id_prefix();
    let id = move |name: &str| format!("{prefix}{name}");

    let material_control = move || {
        view! {
            <div class="control-group">
                <label for=id("material")>"Material"</label>
                // prop:value (not only the options' `selected`) so the shown
                // choice follows every change of the setting: Reset, sync.
                <select
                    id=id("material")
                    prop:value=move || material_to_string(material.get())
                    on:change=move |ev| {
                        let value = event_target_value(&ev);
                        let mat = match value.as_str() {
                            "ceramic" => MugMaterial::Ceramic,
                            "glass" => MugMaterial::Glass,
                            "stainless" => MugMaterial::StainlessSteel,
                            "plastic" => MugMaterial::Plastic,
                            _ => MugMaterial::Ceramic,
                        };
                        settings.set_material(mat);
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
        }
    };

    let volume_control = move || {
        view! {
            <div class="control-group">
                <label for=id("volume")>"Volume: " {move || format!("{:.0} ml", volume_ml.get())}</label>
                <input
                    type="range"
                    id=id("volume")
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
        }
    };

    let thickness_control = move || {
        view! {
            <div class="control-group">
                <label for=id("thickness")>"Wall Thickness: " {move || format!("{:.1} mm", wall_thickness_mm.get())}</label>
                <input
                    type="range"
                    id=id("thickness")
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
        }
    };

    let h_conv_control = move || {
        view! {
            <div class="control-group">
                <label for=id("h-conv")>"Convection: " {move || format!("{:.0} W/(m\u{00B2}\u{00B7}K)", h_conv.get())}</label>
                <input
                    type="range"
                    id=id("h-conv")
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
        }
    };

    let emissivity_control = move || {
        view! {
            <div class="control-group">
                <label for=id("emissivity")>"Emissivity: " {move || format!("{:.2}", emissivity.get())}</label>
                <input
                    type="range"
                    id=id("emissivity")
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
        }
    };

    let exposure_control = move || {
        view! {
            <div class="control-group">
                <label for=id("exposure")>"Lid opening: " {move || format!("{:.0}%", exposure_pct.get())}</label>
                <input
                    type="range"
                    id=id("exposure")
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
        }
    };

    let unit_toggle = move || {
        view! {
            <button
                class="unit-toggle"
                on:click=move |_| {
                    set_use_fahrenheit.set(!use_fahrenheit.get());
                    track_first_interaction();
                }
            >
                {move || if use_fahrenheit.get() { "°F → °C" } else { "°C → °F" }}
            </button>
        }
    };

    // One temperature slider: the setting is stored in °C, the slider shows
    // the current unit; `range_c`/`range_f` are its (min, max) per unit.
    let temp_control = move |name: &'static str,
                             label: &'static str,
                             (temp_c, set_temp_c): (ReadSignal<f64>, WriteSignal<f64>),
                             range_c: (&'static str, &'static str),
                             range_f: (&'static str, &'static str)| {
        view! {
            <div class="control-group">
                <label for=id(name)>{label} {move || display_temp(temp_c.get())}</label>
                <input
                    type="range"
                    id=id(name)
                    min=move || if use_fahrenheit.get() { range_f.0 } else { range_c.0 }
                    max=move || if use_fahrenheit.get() { range_f.1 } else { range_c.1 }
                    step="1"
                    prop:value=move || {
                        if use_fahrenheit.get() {
                            c_to_f(temp_c.get()).round()
                        } else {
                            temp_c.get().round()
                        }
                    }
                    on:input=move |ev| {
                        if let Ok(v) = event_target_value(&ev).parse::<f64>() {
                            let c = if use_fahrenheit.get() { f_to_c(v).round() } else { v };
                            set_temp_c.set(c);
                            vibrate_tick();
                        }
                    }
                />
            </div>
        }
    };
    let temperature_controls = move || {
        view! {
            {temp_control("coffee-temp", "Coffee: ", (coffee_temp_c, set_coffee_temp_c), ("0", "100"), ("32", "212"))}
            {temp_control("preheat-temp", "Preheat mug: ", (preheat_temp_c, set_preheat_temp_c), ("-10", "110"), ("14", "230"))}
            {temp_control("room-temp", "Room: ", (room_temp_c, set_room_temp_c), ("-10", "50"), ("14", "122"))}
            {temp_control("min-drinkable", "Too cold below: ", (min_drinkable_c, set_min_drinkable_c), ("50", "70"), ("122", "158"))}
        }
    };

    let duration_control = move || {
        view! {
            <div class="control-group">
                <label for=id("duration")>"Duration: " {move || format!("{} min", duration_minutes.get())}</label>
                <input
                    type="range"
                    id=id("duration")
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
        }
    };

    let reset_button = move || {
        view! {
            <button
                class="reset-button"
                on:click=move |_| {
                    // One batch: the shared reset and the panel's own reset
                    // land as a single update, so dependents re-solve once.
                    batch(|| {
                        settings.reset_to_defaults();
                        if let Some(cb) = on_reset {
                            cb.call(());
                        }
                    });
                }
            >
                "Reset to defaults"
            </button>
        }
    };

    match layout {
        ControlsLayout::Flat => view! {
            <div class="controls">
                <h3>"Mug Properties"</h3>
                {material_control()}
                {volume_control()}
                {thickness_control()}
                {h_conv_control()}
                {emissivity_control()}
                {exposure_control()}

                <div class="control-group-header">
                    <h3>"Temperatures"</h3>
                    {unit_toggle()}
                </div>
                {temperature_controls()}

                <h3>"Simulation"</h3>
                {duration_control()}

                {children.map(|c| c())}

                {reset_button()}
            </div>
        }
        .into_view(),
        ControlsLayout::Accordion => view! {
            <div class="controls-accordion">
                <AccordionGroup title="Mug Properties">
                    {material_control()}
                    {volume_control()}
                    {thickness_control()}
                </AccordionGroup>

                <AccordionGroup title="Temperatures" initially_open=true>
                    {temperature_controls()}
                    {unit_toggle()}
                </AccordionGroup>

                <AccordionGroup title="Simulation">
                    {duration_control()}
                </AccordionGroup>

                <AccordionGroup title="Advanced">
                    {h_conv_control()}
                    {emissivity_control()}
                    {exposure_control()}
                    {children.map(|c| c())}
                </AccordionGroup>

                <div style="padding: 0.75rem 1rem;">
                    {reset_button()}
                </div>
            </div>
        }
        .into_view(),
    }
}

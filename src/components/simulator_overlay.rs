use crate::components::overlay_history::{
    pop_overlay_entry, push_overlay_entry, set_overlay_closer,
};
use crate::components::{ControlsLayout, MugControls, SimTab};
use crate::storage::SharedSettings;
use leptos::*;
use wasm_bindgen::JsCast;

/// Elements Tab can reach inside the overlay (the e2e spec uses the same).
const FOCUSABLE: &str = "button:not([disabled]), a[href], input:not([disabled]), \
    select:not([disabled]), textarea:not([disabled]), [tabindex]:not([tabindex=\"-1\"])";

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

/// The overlay's Tab stops in order, leaving out those not rendered (e.g.
/// inside a collapsed accordion group).
fn focusable_elements(dialog: &web_sys::Element) -> Vec<web_sys::HtmlElement> {
    let Ok(nodes) = dialog.query_selector_all(FOCUSABLE) else {
        return Vec::new();
    };
    (0..nodes.length())
        .filter_map(|i| nodes.item(i))
        .filter_map(|node| node.dyn_into::<web_sys::HtmlElement>().ok())
        .filter(|el| el.get_client_rects().length() > 0)
        .collect()
}

/// Keep Tab and Shift+Tab inside the dialog: wrap from the last stop to the
/// first and back, and pull focus in if it is on the dialog itself.
fn trap_tab(ev: &ev::KeyboardEvent, dialog: &web_sys::Element) {
    let stops = focusable_elements(dialog);
    let (Some(first), Some(last)) = (stops.first(), stops.last()) else {
        ev.prevent_default();
        return;
    };
    let active = document().active_element();
    let at = |el: &web_sys::HtmlElement| active.as_ref() == Some(el.unchecked_ref());
    let inside_stops = stops.iter().any(at);
    let target = if ev.shift_key() {
        (at(first) || !inside_stops).then_some(last)
    } else {
        (at(last) || !inside_stops).then_some(first)
    };
    if let Some(target) = target {
        ev.prevent_default();
        let _ = target.focus();
    }
}

/// The mobile simulator: a full-screen modal dialog. Opening it moves focus
/// in and pushes a history entry; Done, Escape and the Back button close it,
/// and focus returns to `return_focus` (the "Open Simulator" button).
#[component]
pub fn SimulatorOverlay(
    show: ReadSignal<bool>,
    set_show: WriteSignal<bool>,
    active_tab: RwSignal<String>,
    return_focus: NodeRef<html::Button>,
) -> impl IntoView {
    let settings = SharedSettings::instance();

    let dialog_ref = create_node_ref::<html::Div>();

    // Open: lock body scroll, push the history entry Back pops, and move focus
    // into the dialog once it is in the document. Close: unlock and return
    // focus to the opener. Acts on transitions only, since the close paths can
    // set `show` to false more than once.
    create_effect(move |was_open: Option<bool>| {
        let open = show.get();
        if was_open.unwrap_or(false) != open {
            set_body_scroll_lock(open);
            if open {
                push_overlay_entry();
                request_animation_frame(move || {
                    if let Some(dialog) = dialog_ref.get_untracked() {
                        let _ = dialog.focus();
                    }
                });
            } else if let Some(button) = return_focus.get_untracked() {
                let _ = button.focus();
            }
        }
        open
    });

    // Back pops the entry: the popstate listener closes the overlay through
    // this. Unmounted (left the mobile layout), the entry stays and the next
    // Back pops it without closing anything.
    set_overlay_closer(Some(Callback::new(move |()| set_show.set(false))));
    on_cleanup(move || {
        set_overlay_closer(None);
        set_body_scroll_lock(false);
    });

    // Done and Escape close here and step back over the history entry; the
    // resulting popstate is consumed without a reload.
    let close = move || {
        set_show.set(false);
        pop_overlay_entry();
    };
    let on_keydown = move |ev: ev::KeyboardEvent| match ev.key().as_str() {
        "Escape" => {
            ev.prevent_default();
            close();
        }
        "Tab" => {
            if let Some(dialog) = dialog_ref.get_untracked() {
                trap_tab(&ev, &dialog);
            }
        }
        _ => {}
    };

    // The shell and the shared controls are built once per opening; only the
    // panel and the active tab's own controls follow the tab, so a tab switch
    // keeps the accordion state and the controls' scroll position.
    let panel = move || SimTab::from_key(&active_tab.get()).panel(false);
    let extra_controls =
        move || SimTab::from_key(&active_tab.get()).extra_controls(ControlsLayout::Accordion);
    let reset_extras = move |_| SimTab::from_key(&active_tab.get_untracked()).reset_own_settings();

    view! {
        <Show when=move || show.get()>
            <div
                class="simulator-overlay"
                role="dialog"
                aria-modal="true"
                aria-labelledby="sim-overlay-title"
                tabindex="-1"
                node_ref=dialog_ref
                on:keydown=on_keydown
            >
                <div class="overlay-header">
                    <h2 id="sim-overlay-title" class="overlay-title">"Simulator"</h2>
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
                    <button class="close-btn" on:click=move |_| close()>"Done"</button>
                </div>

                <div class="overlay-chart">
                    {panel}
                </div>

                <div class="overlay-controls">
                    <MugControls
                        settings=settings
                        layout=ControlsLayout::Accordion
                        on_reset=Callback::new(reset_extras)
                    >
                        {extra_controls}
                    </MugControls>
                </div>
            </div>
        </Show>
    }
}

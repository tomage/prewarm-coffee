pub mod analytics;
pub mod components;
pub mod haptics;
pub mod pages;
pub mod physics;
pub mod storage;
pub mod timing;

#[cfg(test)]
mod reactive_test;

use components::{adopt_overlay_entry, overlay_popstate_handled, App, NavBar, PanelSettings};
use leptos::prelude::*;
use leptos_router::components::{Route, Router, Routes};
use leptos_router::hooks::use_location;
use leptos_router::path;
use pages::{AboutPage, NotFoundPage, TheoryHeatEquationPage, TheoryLumpedPage, TheoryPage};
use storage::{
    flush_pending_save_on_pagehide, listen_for_cross_tab_sync, needs_save, schedule_save,
    SharedSettings,
};
use wasm_bindgen::closure::Closure;
use wasm_bindgen::prelude::*;

/// Browser Back/forward. The router follows `popstate` itself; this listener
/// runs first (it is registered before the router mounts) and handles what
/// the router does not know about.
///
/// A popstate that pops the mobile simulator overlay's own history entry is
/// not a navigation: it closes the overlay, and the event stops here. The
/// router never saw that entry pushed; letting it see the pop would count it
/// as a Back in its own path stack. Any other popstate is a navigation:
/// scroll restoration saves the page being left and restores the
/// destination's position, and an overlay entry it lands on is adopted, as
/// one a page loads onto is.
fn setup_history_navigation() {
    adopt_overlay_entry();
    let closure = Closure::wrap(Box::new(move |ev: web_sys::Event| {
        if overlay_popstate_handled() {
            ev.stop_immediate_propagation();
            return;
        }
        adopt_overlay_entry();
        storage::on_history_navigation();
    }) as Box<dyn Fn(web_sys::Event)>);

    if let Some(window) = web_sys::window() {
        let _ =
            window.add_event_listener_with_callback("popstate", closure.as_ref().unchecked_ref());
    }
    closure.forget();
}

/// Tells scroll restoration about client-side route changes.
/// Inside `<Router>`, which `use_location` needs.
#[component]
fn ScrollTracker() -> impl IntoView {
    let location = use_location();
    Effect::new(move |prev: Option<String>| {
        let path = location.pathname.get();
        if prev.is_some_and(|prev| prev != path) {
            storage::on_route_change(&path);
        }
        path
    });
}

/// Root component: initializes the shared settings, persists them, and
/// declares the routes.
#[component]
fn Root() -> impl IntoView {
    // Shared settings singleton — same RwSignals for all routes
    let settings = SharedSettings::instance();
    // Created here so their signals belong to the app root: created inside a
    // panel, they would be disposed when the route unmounts that panel.
    PanelSettings::instance();

    // Save to localStorage when the settings differ from the last synced
    // snapshot (loaded, received from another tab, or saved), debounced so a
    // slider drag is one write. Values read from storage therefore never echo
    // back, nothing is written at startup (stored data, possibly from a newer
    // version, stays intact until a setting changes), and a `?tab=` link
    // selecting a different tab than the stored one counts as a change.
    Effect::new(move |_| {
        if needs_save(&settings.to_stored()) {
            schedule_save();
        }
    });
    flush_pending_save_on_pagehide();

    // Cross-tab sync: storage events, plus a re-read on focus and
    // visibilitychange for events missed while the tab was in the background
    listen_for_cross_tab_sync();

    view! {
        <Router>
            <ScrollTracker/>
            <NavBar/>
            <Routes fallback=NotFoundPage>
                <Route path=path!("/") view=App/>
                <Route path=path!("/about") view=AboutPage/>
                <Route path=path!("/theory") view=TheoryPage/>
                <Route path=path!("/theory/lumped") view=TheoryLumpedPage/>
                <Route path=path!("/theory/heat-equation") view=TheoryHeatEquationPage/>
            </Routes>
        </Router>
    }
}

/// Mount the application to the DOM
#[wasm_bindgen(start)]
pub fn main() {
    console_error_panic_hook::set_once();
    storage::setup_scroll_restoration();
    // Before mounting, so its popstate listener runs before the router's.
    setup_history_navigation();
    mount_to_body(Root);
    storage::restore_scroll_on_load();
}

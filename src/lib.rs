pub mod analytics;
pub mod components;
pub mod haptics;
pub mod pages;
pub mod physics;
pub mod storage;
pub mod timing;

use components::{
    adopt_overlay_entry_on_load, overlay_popstate_handled, App, NavBar, PanelSettings,
};
use leptos::*;
use leptos_router::*;
use pages::{AboutPage, NotFoundPage, TheoryHeatEquationPage, TheoryLumpedPage, TheoryPage};
use storage::{
    flush_pending_save_on_pagehide, is_applying_external, listen_for_cross_tab_sync, schedule_save,
    SharedSettings,
};
use wasm_bindgen::closure::Closure;
use wasm_bindgen::prelude::*;

/// Workaround for Leptos 0.6 router not re-rendering on browser back/forward.
///
/// On `popstate`, the router updates its internal location signal but doesn't
/// always trigger the `<Routes>` component to re-evaluate which view to show.
/// We work around this by forcing a full page reload on popstate, which
/// re-initializes the WASM app at the correct URL. State is preserved via
/// localStorage (SharedSettings), and the page's scroll position is restored
/// after the reload (see `storage::prepare_popstate_reload`).
///
/// A popstate that pops the mobile simulator overlay's own history entry is
/// not a navigation: it closes the overlay and the page stays. A page loaded
/// onto that entry (a reload with the overlay open) adopts it first.
fn setup_popstate_reload() {
    adopt_overlay_entry_on_load();
    let closure = Closure::wrap(Box::new(move |_: web_sys::Event| {
        if overlay_popstate_handled() {
            return;
        }
        storage::prepare_popstate_reload();
        if let Some(window) = web_sys::window() {
            let _ = window.location().reload();
        }
    }) as Box<dyn Fn(web_sys::Event)>);

    if let Some(window) = web_sys::window() {
        let _ =
            window.add_event_listener_with_callback("popstate", closure.as_ref().unchecked_ref());
    }
    closure.forget();
}

/// Tells scroll restoration about client-side route changes (link clicks).
/// Inside `<Router>`, which `use_location` needs.
#[component]
fn ScrollTracker() -> impl IntoView {
    let location = use_location();
    create_effect(move |prev: Option<String>| {
        let path = location.pathname.get();
        if prev.is_some_and(|prev| prev != path) {
            storage::on_route_change(&path);
        }
        path
    });
}

/// Root component with routing
#[component]
fn Root() -> impl IntoView {
    view! {
        <ErrorBoundary fallback=|errors| view! {
            <main class="container">
                <div class="error-container">
                    <h2>"Something went wrong"</h2>
                    <p>"The simulation encountered an error. Try refreshing the page or resetting to defaults."</p>
                    <ul>
                        {move || errors.get()
                            .into_iter()
                            .map(|(_, e)| view! { <li>{e.to_string()}</li> })
                            .collect_view()
                        }
                    </ul>
                    <button on:click=move |_| {
                        // Clear stored settings and reload
                        if let Some(storage) = web_sys::window()
                            .and_then(|w| w.local_storage().ok().flatten())
                        {
                            let _ = storage.remove_item("prewarm_settings");
                        }
                        if let Some(window) = web_sys::window() {
                            let _ = window.location().reload();
                        }
                    }>"Reset & Reload"</button>
                </div>
            </main>
        }>
            <RootInner/>
        </ErrorBoundary>
    }
}

/// Inner root that initializes settings and routes.
/// Wrapped by ErrorBoundary so initialization panics are caught.
#[component]
fn RootInner() -> impl IntoView {
    // Shared settings singleton — same RwSignals for all routes
    let settings = SharedSettings::instance();
    // Created here so their signals belong to the app root: created inside a
    // panel, they would be disposed when the route unmounts that panel.
    PanelSettings::instance();

    // Save to localStorage whenever any signal changes, debounced so a slider
    // drag is one write (skip when applying changes received from another tab).
    // The effect runs once on creation; skip that run so nothing is written
    // until a setting changes, and stored data (possibly from a newer version)
    // stays intact until then; the first save merges into it. A `?tab=` link
    // selecting a different tab than the stored one is such a change.
    create_effect(move |prev: Option<()>| {
        settings.to_stored(); // track every setting
        if prev.is_some() && !is_applying_external() {
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
            <Routes>
                <Route path="/" view=App/>
                <Route path="/about" view=AboutPage/>
                <Route path="/theory" view=TheoryPage/>
                <Route path="/theory/lumped" view=TheoryLumpedPage/>
                <Route path="/theory/heat-equation" view=TheoryHeatEquationPage/>
                <Route path="/*" view=NotFoundPage/>
            </Routes>
        </Router>
    }
}

/// Mount the application to the DOM
#[wasm_bindgen(start)]
pub fn main() {
    console_error_panic_hook::set_once();
    storage::setup_scroll_restoration();
    setup_popstate_reload();
    mount_to_body(Root);
    storage::restore_scroll_on_load();
}

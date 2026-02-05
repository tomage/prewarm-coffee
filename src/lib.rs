pub mod analytics;
pub mod components;
pub mod haptics;
pub mod pages;
pub mod physics;
pub mod storage;

use components::{App, NavBar};
use leptos::*;
use leptos_router::*;
use pages::{AboutPage, NotFoundPage, TheoryLumpedPage, TheoryHeatEquationPage, TheoryPage};
use storage::{is_applying_external, listen_for_cross_tab_sync, save_settings, SharedSettings};
use wasm_bindgen::prelude::*;
use wasm_bindgen::closure::Closure;

/// Workaround for Leptos 0.6 router not re-rendering on browser back/forward.
///
/// On `popstate`, the router updates its internal location signal but doesn't
/// always trigger the `<Routes>` component to re-evaluate which view to show.
/// We work around this by forcing a full page reload on popstate, which
/// re-initializes the WASM app at the correct URL. State is preserved via
/// localStorage (SharedSettings) and scroll restoration.
fn setup_popstate_reload() {
    let closure = Closure::wrap(Box::new(move |_: web_sys::Event| {
        if let Some(window) = web_sys::window() {
            let _ = window.location().reload();
        }
    }) as Box<dyn Fn(web_sys::Event)>);

    if let Some(window) = web_sys::window() {
        let _ = window.add_event_listener_with_callback(
            "popstate",
            closure.as_ref().unchecked_ref(),
        );
    }
    closure.forget();
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

    // Save to localStorage whenever any signal changes
    // (skip when applying changes received from another tab)
    create_effect(move |_| {
        let stored = settings.to_stored();
        if !is_applying_external() {
            save_settings(&stored);
        }
    });

    // Cross-tab sync: storage events + visibilitychange fallback
    listen_for_cross_tab_sync();

    view! {
        <Router>
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
    // Restore scroll after a brief delay to ensure content has rendered
    storage::restore_scroll_after_delay(50);
}

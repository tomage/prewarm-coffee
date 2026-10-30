use leptos::*;
use wasm_bindgen::closure::Closure;
use wasm_bindgen::JsCast;

/// Viewports that get the mobile simulator (entry card + fullscreen
/// overlay) instead of the inline one. Same breakpoint as the mobile rules in
/// `static/styles.css`.
pub const MOBILE_QUERY: &str = "(max-width: 700px)";

/// Whether the viewport matches the media `query`, following changes (resize,
/// rotation). The listener is removed when the calling scope is disposed.
pub fn use_media_query(query: &'static str) -> Signal<bool> {
    let list = window()
        .match_media(query)
        .expect("matchMedia accepts the query")
        .expect("the browser supports matchMedia");
    let (matches, set_matches) = create_signal(list.matches());

    let on_change = {
        let list = list.clone();
        Closure::<dyn Fn()>::new(move || set_matches.set(list.matches()))
    };
    list.add_event_listener_with_callback("change", on_change.as_ref().unchecked_ref())
        .expect("a MediaQueryList accepts a change listener");
    on_cleanup(move || {
        // Remove the listener before dropping the closure: a change event
        // reaching a dropped closure throws.
        let _ =
            list.remove_event_listener_with_callback("change", on_change.as_ref().unchecked_ref());
        drop(on_change);
    });

    matches.into()
}

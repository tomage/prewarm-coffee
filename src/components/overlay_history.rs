//! The mobile simulator overlay's browser-history entry, so the Back button
//! closes the overlay instead of leaving the page.
//!
//! Opening the overlay pushes one entry with state `{overlay: true}` and the
//! same URL, which the router never sees. The `popstate` that pops it is
//! consumed by [`overlay_popstate_handled`]: it closes the overlay, and the
//! listener in `lib.rs` keeps the event from the router. Closing from inside
//! the overlay (Done, Escape) goes `history.back()` over the entry, so the
//! history stack matches what is on screen.
//!
//! "Entry on the stack" is tracked separately from "overlay open": a resize
//! across the mobile breakpoint unmounts the overlay without touching
//! history, and the next Back then pops the leftover entry silently.

use leptos::prelude::*;
use std::cell::RefCell;

thread_local! {
    /// The URL (path and query) of the overlay's entry while it is on the
    /// history stack above its page's entry; `None` otherwise.
    static ENTRY_URL: RefCell<Option<String>> = const { RefCell::new(None) };
    /// Closes the mounted overlay; `None` while no overlay is mounted.
    static CLOSE_OVERLAY: RefCell<Option<Callback<()>>> = const { RefCell::new(None) };
}

/// Push the overlay's history entry, unless it is still on the stack from an
/// earlier opening (then that entry is reused, keeping one entry per overlay).
pub fn push_overlay_entry() {
    if entry_pushed() {
        return;
    }
    let state = js_sys::Object::new();
    let _ = js_sys::Reflect::set(&state, &"overlay".into(), &true.into());
    if window()
        .history()
        .and_then(|h| h.push_state(&state, ""))
        .is_ok()
    {
        set_entry_url(Some(current_url()));
    }
}

/// Adopt the overlay entry the history is on, when the page arrives at it
/// with the overlay closed (a reload with the overlay open, or Back/forward
/// onto an entry whose overlay is gone): the next Back pops that entry
/// silently and the next opening reuses it.
pub fn adopt_overlay_entry() {
    if on_overlay_entry() {
        set_entry_url(Some(current_url()));
    }
}

/// Step back over the overlay's history entry, if it is on the stack. The
/// resulting `popstate` is consumed by [`overlay_popstate_handled`].
pub fn pop_overlay_entry() {
    if entry_pushed() {
        let _ = window().history().and_then(|h| h.back());
    }
}

/// Forget the overlay's entry without popping it: the page that pushed it is
/// gone, so a later `popstate` onto that entry is a real navigation.
pub fn forget_overlay_entry() {
    set_entry_url(None);
}

/// Register how to close the mounted overlay; `None` when it unmounts.
pub fn set_overlay_closer(close: Option<Callback<()>>) {
    CLOSE_OVERLAY.with_borrow_mut(|slot| *slot = close);
}

/// Called first by the `popstate` listener. While the overlay's entry is on
/// the stack, any `popstate` leaves it, and the entry is forgotten. Leaving
/// it for its page's entry (same URL, not an overlay entry) pops the overlay:
/// close it (if still mounted) and return `true`, so the page stays and the
/// router must not see the event. Leaving it anywhere else (Forward past a
/// leftover entry, or a jump of several entries) is a real navigation:
/// `false`, as for every `popstate` while no entry is on the stack.
pub fn overlay_popstate_handled() -> bool {
    let Some(entry_url) = ENTRY_URL.with_borrow_mut(Option::take) else {
        return false;
    };
    if on_overlay_entry() || current_url() != entry_url {
        return false;
    }
    if let Some(close) = CLOSE_OVERLAY.with_borrow(|slot| *slot) {
        close.run(());
    }
    true
}

fn entry_pushed() -> bool {
    ENTRY_URL.with_borrow(Option::is_some)
}

fn set_entry_url(url: Option<String>) {
    ENTRY_URL.with_borrow_mut(|slot| *slot = url);
}

/// Whether the current history entry is an overlay entry.
fn on_overlay_entry() -> bool {
    window()
        .history()
        .and_then(|h| h.state())
        .and_then(|state| js_sys::Reflect::get(&state, &"overlay".into()))
        .is_ok_and(|v| v.as_bool() == Some(true))
}

/// The current URL's path and query: what the overlay's entry shares with
/// its page's entry.
fn current_url() -> String {
    let location = window().location();
    let path = location.pathname().unwrap_or_default();
    let search = location.search().unwrap_or_default();
    format!("{path}{search}")
}

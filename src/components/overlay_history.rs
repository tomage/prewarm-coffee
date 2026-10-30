//! The mobile simulator overlay's browser-history entry, so the Back button
//! closes the overlay instead of leaving (or reloading) the page.
//!
//! Opening the overlay pushes one entry with state `{overlay: true}` and the
//! same URL. The `popstate` that pops it is consumed by
//! [`overlay_popstate_handled`]: it closes the overlay and the reload
//! workaround in `lib.rs` does not run. Closing from inside the overlay (Done,
//! Escape) goes `history.back()` over the entry, so the history stack matches
//! what is on screen.
//!
//! "Entry pushed" is tracked separately from "overlay open": a resize across
//! the mobile breakpoint unmounts the overlay without touching history, and
//! the next Back then pops the leftover entry silently.

use leptos::*;
use std::cell::{Cell, RefCell};

thread_local! {
    /// Whether the overlay's entry is on the history stack above the page's.
    static ENTRY_PUSHED: Cell<bool> = const { Cell::new(false) };
    /// Closes the mounted overlay; `None` while no overlay is mounted.
    static CLOSE_OVERLAY: RefCell<Option<Callback<()>>> = const { RefCell::new(None) };
}

/// Push the overlay's history entry, unless it is still on the stack from an
/// earlier opening (then that entry is reused, keeping one entry per overlay).
pub fn push_overlay_entry() {
    if ENTRY_PUSHED.get() {
        return;
    }
    let state = js_sys::Object::new();
    let _ = js_sys::Reflect::set(&state, &"overlay".into(), &true.into());
    if window()
        .history()
        .and_then(|h| h.push_state(&state, ""))
        .is_ok()
    {
        ENTRY_PUSHED.set(true);
    }
}

/// Adopt an overlay entry the page loaded onto (a reload with the overlay
/// open): the overlay starts closed, so the next Back pops that entry silently
/// and the next opening reuses it.
pub fn adopt_overlay_entry_on_load() {
    let on_overlay_entry = window()
        .history()
        .and_then(|h| h.state())
        .and_then(|state| js_sys::Reflect::get(&state, &"overlay".into()))
        .is_ok_and(|v| v.as_bool() == Some(true));
    if on_overlay_entry {
        ENTRY_PUSHED.set(true);
    }
}

/// Step back over the overlay's history entry, if it is on the stack. The
/// resulting `popstate` is consumed by [`overlay_popstate_handled`].
pub fn pop_overlay_entry() {
    if ENTRY_PUSHED.get() {
        let _ = window().history().and_then(|h| h.back());
    }
}

/// Forget the overlay's entry without popping it: the page that pushed it is
/// gone, so a later `popstate` onto that entry is a real navigation and must
/// take the reload path.
pub fn forget_overlay_entry() {
    ENTRY_PUSHED.set(false);
}

/// Register how to close the mounted overlay; `None` when it unmounts.
pub fn set_overlay_closer(close: Option<Callback<()>>) {
    CLOSE_OVERLAY.with_borrow_mut(|slot| *slot = close);
}

/// Called first by the `popstate` listener. When the overlay's entry was
/// pushed, this popstate pops it: close the overlay (if still mounted) and
/// return `true` so no reload happens. Otherwise `false`.
pub fn overlay_popstate_handled() -> bool {
    if !ENTRY_PUSHED.replace(false) {
        return false;
    }
    if let Some(close) = CLOSE_OVERLAY.with_borrow(|slot| *slot) {
        close.call(());
    }
    true
}

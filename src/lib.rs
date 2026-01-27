pub mod analytics;
pub mod components;
pub mod haptics;
pub mod physics;
pub mod storage;

use components::App;
use leptos::*;
use wasm_bindgen::prelude::*;

/// Mount the application to the DOM
#[wasm_bindgen(start)]
pub fn main() {
    console_error_panic_hook::set_once();
    storage::setup_scroll_restoration();
    mount_to_body(App);
    // Restore scroll after a brief delay to ensure content has rendered
    storage::restore_scroll_after_delay(50);
}

use wasm_bindgen::prelude::*;

#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(js_namespace = umami, js_name = track)]
    fn umami_track(event: &str);
}

/// Track a custom event in Umami analytics
/// Fails silently if Umami is not loaded (e.g., blocked by adblocker)
pub fn track_event(event: &str) {
    // Check if umami exists before calling
    let result = js_sys::eval("typeof umami !== 'undefined'");
    if let Ok(val) = result {
        if val.as_bool().unwrap_or(false) {
            umami_track(event);
        }
    }
}

use web_sys::window;

/// Trigger a short haptic feedback vibration (if supported)
pub fn vibrate_tick() {
    if let Some(window) = window() {
        let _ = window.navigator().vibrate_with_duration(10);
    }
}

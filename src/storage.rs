use crate::physics::MugMaterial;
use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

const STORAGE_KEY: &str = "prewarm_settings";
const SCROLL_KEY: &str = "prewarm_scroll";

#[derive(Serialize, Deserialize)]
pub struct StoredSettings {
    pub material: String,
    pub volume_ml: f64,
    pub wall_thickness_mm: f64,
    pub coffee_temp_c: f64,
    pub room_temp_c: f64,
    #[serde(default = "default_min_drinkable")]
    pub min_drinkable_c: f64,
    #[serde(default = "default_use_fahrenheit")]
    pub use_fahrenheit: bool,
    #[serde(default = "default_duration_minutes")]
    pub duration_minutes: usize,
    #[serde(default = "default_preheat_temp")]
    pub preheat_temp_c: f64,
}

fn default_use_fahrenheit() -> bool {
    true
}

fn default_duration_minutes() -> usize {
    60
}

fn default_min_drinkable() -> f64 {
    60.0
}

fn default_preheat_temp() -> f64 {
    90.0
}

impl Default for StoredSettings {
    fn default() -> Self {
        Self {
            material: "ceramic".to_string(),
            volume_ml: 180.0,
            wall_thickness_mm: 5.0,
            coffee_temp_c: 90.0,
            room_temp_c: 22.0,
            min_drinkable_c: 60.0,
            use_fahrenheit: true,
            duration_minutes: 60,
            preheat_temp_c: 90.0,
        }
    }
}

/// Convert Celsius to Fahrenheit
pub fn c_to_f(c: f64) -> f64 {
    c * 9.0 / 5.0 + 32.0
}

/// Convert Fahrenheit to Celsius
pub fn f_to_c(f: f64) -> f64 {
    (f - 32.0) * 5.0 / 9.0
}

impl StoredSettings {
    pub fn material_enum(&self) -> MugMaterial {
        match self.material.as_str() {
            "ceramic" => MugMaterial::Ceramic,
            "glass" => MugMaterial::Glass,
            "stainless" => MugMaterial::StainlessSteel,
            "plastic" => MugMaterial::Plastic,
            _ => MugMaterial::Ceramic,
        }
    }
}

pub fn material_to_string(mat: MugMaterial) -> String {
    match mat {
        MugMaterial::Ceramic => "ceramic",
        MugMaterial::Glass => "glass",
        MugMaterial::StainlessSteel => "stainless",
        MugMaterial::Plastic => "plastic",
    }
    .to_string()
}

/// Attempts to get the browser's localStorage.
///
/// Returns `None` if:
/// - Running outside a browser environment
/// - localStorage is disabled (e.g., private/incognito browsing mode)
/// - The browser's storage quota has been exceeded
/// - A SecurityError occurs (e.g., third-party iframe with restricted access)
fn get_storage() -> Option<web_sys::Storage> {
    web_sys::window()?.local_storage().ok()?
}

/// Log a warning message to the browser console.
///
/// Used for non-fatal storage errors that the user might want to know about.
fn log_warning(msg: &str) {
    web_sys::console::warn_1(&msg.into());
}

/// Load user settings from localStorage.
///
/// # Fallback Behavior
///
/// Returns default settings if:
/// - localStorage is unavailable (private browsing, disabled, etc.)
/// - The stored JSON is corrupted or from an incompatible version
/// - No settings have been saved yet
///
/// This graceful degradation ensures the app always works, even without persistence.
pub fn load_settings() -> StoredSettings {
    match get_storage() {
        Some(storage) => {
            match storage.get_item(STORAGE_KEY) {
                Ok(Some(json)) => {
                    serde_json::from_str(&json).unwrap_or_else(|e| {
                        log_warning(&format!(
                            "Prewarm: Failed to parse stored settings (using defaults): {}",
                            e
                        ));
                        StoredSettings::default()
                    })
                }
                Ok(None) => StoredSettings::default(), // No saved settings yet
                Err(_) => {
                    log_warning("Prewarm: Could not read from localStorage (using defaults)");
                    StoredSettings::default()
                }
            }
        }
        None => {
            // localStorage unavailable - likely private browsing mode
            // This is expected behavior, so we don't warn
            StoredSettings::default()
        }
    }
}

/// Save user settings to localStorage.
///
/// # Error Handling
///
/// Silently fails (with console warning) if:
/// - localStorage is unavailable (private/incognito browsing mode)
/// - Storage quota has been exceeded
/// - Serialization fails (should not happen with valid settings)
///
/// The app continues to function without persistence in these cases.
/// Users in private browsing mode will simply start with defaults each visit.
pub fn save_settings(settings: &StoredSettings) {
    let storage = match get_storage() {
        Some(s) => s,
        None => return, // localStorage unavailable, silently skip (expected in private mode)
    };

    let json = match serde_json::to_string(settings) {
        Ok(j) => j,
        Err(e) => {
            log_warning(&format!("Prewarm: Failed to serialize settings: {}", e));
            return;
        }
    };

    if storage.set_item(STORAGE_KEY, &json).is_err() {
        log_warning(
            "Prewarm: Could not save settings to localStorage. \
             You may be in private browsing mode or storage quota exceeded.",
        );
    }
}

/// Set up scroll restoration - call this once on app init
pub fn setup_scroll_restoration() {
    // Disable browser's automatic scroll restoration
    if let Some(window) = web_sys::window() {
        if let Ok(history) = window.history() {
            let _ = history.set_scroll_restoration(web_sys::ScrollRestoration::Manual);
        }
    }

    // Save scroll position before unload
    if let Some(window) = web_sys::window() {
        let closure = Closure::wrap(Box::new(move || {
            if let Some(win) = web_sys::window() {
                let scroll_y = win.scroll_y().unwrap_or(0.0);
                if let Some(storage) = win.session_storage().ok().flatten() {
                    let _ = storage.set_item(SCROLL_KEY, &scroll_y.to_string());
                }
            }
        }) as Box<dyn Fn()>);

        let _ = window
            .add_event_listener_with_callback("beforeunload", closure.as_ref().unchecked_ref());
        closure.forget(); // Keep the closure alive
    }
}

/// Restore scroll position - call this after app has rendered
pub fn restore_scroll_position() {
    if let Some(window) = web_sys::window() {
        if let Ok(Some(storage)) = window.session_storage() {
            if let Ok(Some(scroll_str)) = storage.get_item(SCROLL_KEY) {
                if let Ok(scroll_y) = scroll_str.parse::<f64>() {
                    window.scroll_to_with_x_and_y(0.0, scroll_y);
                }
            }
        }
    }
}

/// Restore scroll position after a delay (in ms)
pub fn restore_scroll_after_delay(delay_ms: i32) {
    if let Some(window) = web_sys::window() {
        let closure = Closure::once(Box::new(|| {
            restore_scroll_position();
        }) as Box<dyn FnOnce()>);

        let _ = window.set_timeout_with_callback_and_timeout_and_arguments_0(
            closure.as_ref().unchecked_ref(),
            delay_ms,
        );
        closure.forget();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c_to_f_known_values() {
        // Freezing point of water
        assert!((c_to_f(0.0) - 32.0).abs() < 0.001);
        // Boiling point of water
        assert!((c_to_f(100.0) - 212.0).abs() < 0.001);
        // Body temperature
        assert!((c_to_f(37.0) - 98.6).abs() < 0.1);
        // Room temperature
        assert!((c_to_f(22.0) - 71.6).abs() < 0.1);
    }

    #[test]
    fn test_f_to_c_known_values() {
        // Freezing point of water
        assert!((f_to_c(32.0) - 0.0).abs() < 0.001);
        // Boiling point of water
        assert!((f_to_c(212.0) - 100.0).abs() < 0.001);
        // Body temperature
        assert!((f_to_c(98.6) - 37.0).abs() < 0.1);
    }

    #[test]
    fn test_temperature_conversion_roundtrip() {
        // Converting C -> F -> C should give original value
        let test_temps = [-40.0, 0.0, 22.0, 37.0, 60.0, 90.0, 100.0];
        for temp in test_temps {
            let roundtrip = f_to_c(c_to_f(temp));
            assert!(
                (roundtrip - temp).abs() < 0.0001,
                "Roundtrip failed for {}°C: got {}°C",
                temp,
                roundtrip
            );
        }
    }

    #[test]
    fn test_temperature_conversion_symmetry() {
        // -40 is the same in both scales
        assert!((c_to_f(-40.0) - (-40.0)).abs() < 0.001);
        assert!((f_to_c(-40.0) - (-40.0)).abs() < 0.001);
    }

    #[test]
    fn test_material_enum_conversion() {
        let settings = StoredSettings {
            material: "ceramic".to_string(),
            ..Default::default()
        };
        assert!(matches!(settings.material_enum(), MugMaterial::Ceramic));

        let settings = StoredSettings {
            material: "glass".to_string(),
            ..Default::default()
        };
        assert!(matches!(settings.material_enum(), MugMaterial::Glass));

        let settings = StoredSettings {
            material: "stainless".to_string(),
            ..Default::default()
        };
        assert!(matches!(
            settings.material_enum(),
            MugMaterial::StainlessSteel
        ));

        let settings = StoredSettings {
            material: "plastic".to_string(),
            ..Default::default()
        };
        assert!(matches!(settings.material_enum(), MugMaterial::Plastic));

        // Unknown material defaults to ceramic
        let settings = StoredSettings {
            material: "unknown".to_string(),
            ..Default::default()
        };
        assert!(matches!(settings.material_enum(), MugMaterial::Ceramic));
    }

    #[test]
    fn test_material_to_string_roundtrip() {
        let materials = [
            MugMaterial::Ceramic,
            MugMaterial::Glass,
            MugMaterial::StainlessSteel,
            MugMaterial::Plastic,
        ];
        for mat in materials {
            let s = material_to_string(mat);
            let settings = StoredSettings {
                material: s,
                ..Default::default()
            };
            assert_eq!(
                std::mem::discriminant(&settings.material_enum()),
                std::mem::discriminant(&mat),
                "Roundtrip failed for {:?}",
                mat
            );
        }
    }
}

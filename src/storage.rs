use crate::physics::MugMaterial;
use leptos::*;
use serde::{Deserialize, Serialize};
use std::cell::Cell;
use wasm_bindgen::prelude::*;

const STORAGE_KEY: &str = "prewarm_settings";
const SCROLL_KEY: &str = "prewarm_scroll";

#[derive(Serialize, Deserialize)]
#[serde(default)]
pub struct StoredSettings {
    pub material: String,
    pub volume_ml: f64,
    pub wall_thickness_mm: f64,
    pub coffee_temp_c: f64,
    pub room_temp_c: f64,
    pub min_drinkable_c: f64,
    pub use_fahrenheit: bool,
    pub duration_minutes: usize,
    pub preheat_temp_c: f64,
    pub h_conv: f64,
    pub exposure_pct: f64,
    pub emissivity: f64,
    pub active_tab: String,
    pub use_analytical_lumped: bool,
    pub use_analytical_heat_eq: bool,
    pub use_near_wall: bool,
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
            duration_minutes: 30,
            preheat_temp_c: 80.0,
            h_conv: 20.0,
            exposure_pct: 20.0,
            emissivity: 0.90,
            active_tab: "lumped".to_string(),
            use_analytical_lumped: true,
            use_analytical_heat_eq: true,
            use_near_wall: false,
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

/// Format a temperature in Celsius or Fahrenheit depending on the flag.
///
/// Produces strings like "90°C" or "194°F".
pub fn display_temp(celsius: f64, use_fahrenheit: bool) -> String {
    if use_fahrenheit {
        format!("{:.0}\u{00B0}F", c_to_f(celsius))
    } else {
        format!("{:.0}\u{00B0}C", celsius)
    }
}

/// Format a duration given in minutes.  When the value is below 5 min the
/// seconds are included: "3 min 24 sec".  At or above 5 min only whole
/// minutes are shown: "42 min".
pub fn format_minutes(minutes: f64) -> String {
    if minutes < 5.0 {
        let total_secs = (minutes * 60.0).round() as u64;
        let m = total_secs / 60;
        let s = total_secs % 60;
        if m > 0 {
            format!("{} min {} sec", m, s)
        } else {
            format!("{} sec", s)
        }
    } else {
        format!("{:.0} min", minutes)
    }
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

    /// Clamp all values to physically valid ranges.
    ///
    /// Prevents invalid inputs from corrupted localStorage or manual edits
    /// from causing simulation panics or nonsensical results.
    pub fn validate(&mut self) {
        // Replace NaN/Infinity with defaults before clamping (f64::clamp returns NaN for NaN)
        let defaults = StoredSettings::default();
        if !self.volume_ml.is_finite() { self.volume_ml = defaults.volume_ml; }
        if !self.wall_thickness_mm.is_finite() { self.wall_thickness_mm = defaults.wall_thickness_mm; }
        if !self.coffee_temp_c.is_finite() { self.coffee_temp_c = defaults.coffee_temp_c; }
        if !self.room_temp_c.is_finite() { self.room_temp_c = defaults.room_temp_c; }
        if !self.min_drinkable_c.is_finite() { self.min_drinkable_c = defaults.min_drinkable_c; }
        if !self.preheat_temp_c.is_finite() { self.preheat_temp_c = defaults.preheat_temp_c; }
        if !self.h_conv.is_finite() { self.h_conv = defaults.h_conv; }
        if !self.exposure_pct.is_finite() { self.exposure_pct = defaults.exposure_pct; }
        if !self.emissivity.is_finite() { self.emissivity = defaults.emissivity; }

        self.volume_ml = self.volume_ml.clamp(40.0, 500.0);
        self.wall_thickness_mm = self.wall_thickness_mm.clamp(1.0, 10.0);
        self.coffee_temp_c = self.coffee_temp_c.clamp(0.0, 100.0);
        self.room_temp_c = self.room_temp_c.clamp(-10.0, 50.0);
        self.min_drinkable_c = self.min_drinkable_c.clamp(50.0, 70.0);
        self.duration_minutes = self.duration_minutes.clamp(10, 120);
        self.preheat_temp_c = self.preheat_temp_c.clamp(-10.0, 110.0);
        self.h_conv = self.h_conv.clamp(2.0, 50.0);
        self.exposure_pct = self.exposure_pct.clamp(0.0, 100.0);
        self.emissivity = self.emissivity.clamp(0.0, 1.0);

        // Ensure room temp < coffee temp (otherwise simulation is nonsensical)
        if self.room_temp_c >= self.coffee_temp_c {
            self.room_temp_c = self.coffee_temp_c - 10.0;
        }
        // Ensure min_drinkable is between room and coffee
        if self.min_drinkable_c <= self.room_temp_c {
            self.min_drinkable_c = self.room_temp_c + 5.0;
        }
        if self.min_drinkable_c >= self.coffee_temp_c {
            self.min_drinkable_c = self.coffee_temp_c - 5.0;
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

/// Shared reactive settings — singleton backed by a thread-local cache.
///
/// Each field is an `RwSignal`, which is `Copy` in Leptos 0.6, making
/// the whole struct `Copy`. Signals are created once (from localStorage)
/// and reused across all route components.
#[derive(Clone, Copy)]
pub struct SharedSettings {
    pub material: RwSignal<MugMaterial>,
    pub volume_ml: RwSignal<f64>,
    pub wall_thickness_mm: RwSignal<f64>,
    pub coffee_temp_c: RwSignal<f64>,
    pub room_temp_c: RwSignal<f64>,
    pub min_drinkable_c: RwSignal<f64>,
    pub use_fahrenheit: RwSignal<bool>,
    pub duration_minutes: RwSignal<usize>,
    pub preheat_temp_c: RwSignal<f64>,
    pub h_conv: RwSignal<f64>,
    pub exposure_pct: RwSignal<f64>,
    pub emissivity: RwSignal<f64>,
    pub active_tab: RwSignal<String>,
    pub use_analytical_lumped: RwSignal<bool>,
    pub use_analytical_heat_eq: RwSignal<bool>,
    pub use_near_wall: RwSignal<bool>,
}

thread_local! {
    static SHARED_SETTINGS: Cell<Option<SharedSettings>> = const { Cell::new(None) };
    /// Guard flag: when true, the save effect should skip writing to localStorage.
    /// Set during `listen_for_storage_changes` to break the cross-tab feedback loop.
    static APPLYING_EXTERNAL: Cell<bool> = const { Cell::new(false) };
}

impl SharedSettings {
    /// Get or create the shared settings singleton.
    ///
    /// First call creates `RwSignal`s from localStorage values.
    /// Subsequent calls return the same signal instances.
    pub fn instance() -> Self {
        SHARED_SETTINGS.with(|cell| {
            if let Some(s) = cell.get() {
                s
            } else {
                let saved = load_settings();
                let s = Self {
                    material: create_rw_signal(saved.material_enum()),
                    volume_ml: create_rw_signal(saved.volume_ml),
                    wall_thickness_mm: create_rw_signal(saved.wall_thickness_mm),
                    coffee_temp_c: create_rw_signal(saved.coffee_temp_c),
                    room_temp_c: create_rw_signal(saved.room_temp_c),
                    min_drinkable_c: create_rw_signal(saved.min_drinkable_c),
                    use_fahrenheit: create_rw_signal(saved.use_fahrenheit),
                    duration_minutes: create_rw_signal(saved.duration_minutes),
                    preheat_temp_c: create_rw_signal(saved.preheat_temp_c),
                    h_conv: create_rw_signal(saved.h_conv),
                    exposure_pct: create_rw_signal(saved.exposure_pct),
                    emissivity: create_rw_signal(saved.emissivity),
                    active_tab: create_rw_signal(saved.active_tab.clone()),
                    use_analytical_lumped: create_rw_signal(saved.use_analytical_lumped),
                    use_analytical_heat_eq: create_rw_signal(saved.use_analytical_heat_eq),
                    use_near_wall: create_rw_signal(saved.use_near_wall),
                };
                cell.set(Some(s));
                s
            }
        })
    }

    /// Snapshot all signals into a serializable `StoredSettings`.
    pub fn to_stored(&self) -> StoredSettings {
        StoredSettings {
            material: material_to_string(self.material.get()),
            volume_ml: self.volume_ml.get(),
            wall_thickness_mm: self.wall_thickness_mm.get(),
            coffee_temp_c: self.coffee_temp_c.get(),
            room_temp_c: self.room_temp_c.get(),
            min_drinkable_c: self.min_drinkable_c.get(),
            use_fahrenheit: self.use_fahrenheit.get(),
            duration_minutes: self.duration_minutes.get(),
            preheat_temp_c: self.preheat_temp_c.get(),
            h_conv: self.h_conv.get(),
            exposure_pct: self.exposure_pct.get(),
            emissivity: self.emissivity.get(),
            active_tab: self.active_tab.get(),
            use_analytical_lumped: self.use_analytical_lumped.get(),
            use_analytical_heat_eq: self.use_analytical_heat_eq.get(),
            use_near_wall: self.use_near_wall.get(),
        }
    }

    /// Build `MugParameters` from the current signal values.
    ///
    /// Normalizes `exposure_pct` from 0–100 to 0.0–1.0.
    pub fn make_mug_params(&self, preheated: bool) -> crate::physics::MugParameters {
        crate::physics::MugParameters {
            material: self.material.get(),
            volume_ml: self.volume_ml.get(),
            wall_thickness_mm: self.wall_thickness_mm.get(),
            coffee_temp_c: self.coffee_temp_c.get(),
            room_temp_c: self.room_temp_c.get(),
            preheated,
            preheat_temp_c: self.preheat_temp_c.get(),
            h_conv: self.h_conv.get(),
            emissivity: self.emissivity.get(),
            exposure_pct: self.exposure_pct.get() / 100.0,
        }
    }

    /// Push values from a `StoredSettings` into the shared signals.
    fn update_from_stored(&self, s: &StoredSettings) {
        self.material.set(s.material_enum());
        self.volume_ml.set(s.volume_ml);
        self.wall_thickness_mm.set(s.wall_thickness_mm);
        self.coffee_temp_c.set(s.coffee_temp_c);
        self.room_temp_c.set(s.room_temp_c);
        self.min_drinkable_c.set(s.min_drinkable_c);
        self.use_fahrenheit.set(s.use_fahrenheit);
        self.duration_minutes.set(s.duration_minutes);
        self.preheat_temp_c.set(s.preheat_temp_c);
        self.h_conv.set(s.h_conv);
        self.exposure_pct.set(s.exposure_pct);
        self.emissivity.set(s.emissivity);
        self.active_tab.set(s.active_tab.clone());
        self.use_analytical_lumped.set(s.use_analytical_lumped);
        self.use_analytical_heat_eq.set(s.use_analytical_heat_eq);
        self.use_near_wall.set(s.use_near_wall);
    }
}

/// Returns true when signals are being updated from another tab's storage event.
/// The save effect should skip writing to localStorage in this case.
pub fn is_applying_external() -> bool {
    APPLYING_EXTERNAL.with(|c| c.get())
}

/// Set up cross-tab sync for shared settings.
///
/// Three mechanisms for reliability:
/// 1. `storage` event — fires on other tabs when localStorage changes
/// 2. `visibilitychange` — re-reads localStorage when tab visibility changes
/// 3. `focus` — re-reads localStorage when window/tab gains focus
pub fn listen_for_cross_tab_sync() {
    // Instant sync: storage event fires on OTHER tabs when one tab writes
    let storage_cb = Closure::wrap(Box::new(move |ev: web_sys::StorageEvent| {
        if ev.key().as_deref() != Some(STORAGE_KEY) {
            return;
        }
        if let Some(json) = ev.new_value() {
            if let Ok(stored) = serde_json::from_str::<StoredSettings>(&json) {
                apply_external(&stored);
            }
        }
    }) as Box<dyn Fn(web_sys::StorageEvent)>);

    // Re-read localStorage on tab/window focus changes
    let sync_from_storage = Closure::wrap(Box::new(move || {
        let stored = load_settings();
        apply_external(&stored);
    }) as Box<dyn Fn()>);

    if let Some(window) = web_sys::window() {
        let _ = window
            .add_event_listener_with_callback("storage", storage_cb.as_ref().unchecked_ref());
        // visibilitychange: fires when switching tabs in the same window
        if let Some(doc) = window.document() {
            let _ = doc.add_event_listener_with_callback(
                "visibilitychange",
                sync_from_storage.as_ref().unchecked_ref(),
            );
        }
        // focus: fires when the window gains focus (covers side-by-side windows)
        let _ = window
            .add_event_listener_with_callback("focus", sync_from_storage.as_ref().unchecked_ref());
    }
    storage_cb.forget();
    sync_from_storage.forget();
}

/// Apply settings from an external source (other tab) into the shared signals,
/// suppressing the save effect to avoid a feedback loop.
fn apply_external(stored: &StoredSettings) {
    let settings = SharedSettings::instance();
    APPLYING_EXTERNAL.with(|c| c.set(true));
    settings.update_from_stored(stored);
    APPLYING_EXTERNAL.with(|c| c.set(false));
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
    let mut settings = match get_storage() {
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
                Ok(None) => StoredSettings::default(),
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
    };
    settings.validate();
    settings
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
    fn test_display_temp() {
        assert_eq!(display_temp(90.0, false), "90\u{00B0}C");
        assert_eq!(display_temp(90.0, true), "194\u{00B0}F");
        assert_eq!(display_temp(0.0, false), "0\u{00B0}C");
        assert_eq!(display_temp(0.0, true), "32\u{00B0}F");
        assert_eq!(display_temp(100.0, true), "212\u{00B0}F");
    }

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
    fn test_validate_clamps_out_of_range() {
        let mut s = StoredSettings {
            volume_ml: 1000.0,       // over max 500
            wall_thickness_mm: 0.1,  // under min 1
            coffee_temp_c: 200.0,    // over max 100
            room_temp_c: -50.0,      // under min -10
            h_conv: 100.0,           // over max 50
            exposure_pct: 150.0,     // over max 100
            ..Default::default()
        };
        s.validate();
        assert_eq!(s.volume_ml, 500.0);
        assert_eq!(s.wall_thickness_mm, 1.0);
        assert_eq!(s.coffee_temp_c, 100.0);
        assert_eq!(s.room_temp_c, -10.0);
        assert_eq!(s.h_conv, 50.0);
        assert_eq!(s.exposure_pct, 100.0);
    }

    #[test]
    fn test_validate_fixes_room_temp_above_coffee() {
        let mut s = StoredSettings {
            coffee_temp_c: 80.0,
            room_temp_c: 80.0,  // equal to coffee — invalid
            ..Default::default()
        };
        s.validate();
        assert!(s.room_temp_c < s.coffee_temp_c);
    }

    #[test]
    fn test_validate_fixes_drinkable_outside_range() {
        let mut s = StoredSettings {
            coffee_temp_c: 90.0,
            room_temp_c: 22.0,
            min_drinkable_c: 50.0, // below min clamp of 50 — will clamp to 50
            ..Default::default()
        };
        s.validate();
        assert!(s.min_drinkable_c > s.room_temp_c);
        assert!(s.min_drinkable_c < s.coffee_temp_c);
    }

    #[test]
    fn test_validate_handles_nan_and_infinity() {
        let mut s = StoredSettings {
            volume_ml: f64::NAN,
            wall_thickness_mm: f64::INFINITY,
            coffee_temp_c: f64::NEG_INFINITY,
            h_conv: f64::NAN,
            ..Default::default()
        };
        s.validate();
        let defaults = StoredSettings::default();
        assert_eq!(s.volume_ml, defaults.volume_ml);
        assert_eq!(s.wall_thickness_mm, defaults.wall_thickness_mm);
        assert_eq!(s.coffee_temp_c, defaults.coffee_temp_c);
        assert_eq!(s.h_conv, defaults.h_conv);
    }

    #[test]
    fn test_validate_default_is_valid() {
        let mut s = StoredSettings::default();
        let before = StoredSettings::default();
        s.validate();
        assert_eq!(s.volume_ml, before.volume_ml);
        assert_eq!(s.coffee_temp_c, before.coffee_temp_c);
        assert_eq!(s.room_temp_c, before.room_temp_c);
        assert_eq!(s.min_drinkable_c, before.min_drinkable_c);
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

    #[test]
    fn test_serde_roundtrip() {
        let original = StoredSettings::default();
        let json = serde_json::to_string(&original).unwrap();
        let restored: StoredSettings = serde_json::from_str(&json).unwrap();
        assert_eq!(restored.volume_ml, original.volume_ml);
        assert_eq!(restored.wall_thickness_mm, original.wall_thickness_mm);
        assert_eq!(restored.coffee_temp_c, original.coffee_temp_c);
        assert_eq!(restored.room_temp_c, original.room_temp_c);
        assert_eq!(restored.min_drinkable_c, original.min_drinkable_c);
        assert_eq!(restored.use_fahrenheit, original.use_fahrenheit);
        assert_eq!(restored.duration_minutes, original.duration_minutes);
        assert_eq!(restored.preheat_temp_c, original.preheat_temp_c);
        assert_eq!(restored.h_conv, original.h_conv);
        assert_eq!(restored.exposure_pct, original.exposure_pct);
        assert_eq!(restored.material, original.material);
    }

    #[test]
    fn test_serde_missing_fields_use_defaults() {
        // Simulates loading settings saved by an older version without newer fields
        let minimal_json = r#"{"material":"glass","volume_ml":200.0,"wall_thickness_mm":4.0,"coffee_temp_c":85.0,"room_temp_c":20.0}"#;
        let restored: StoredSettings = serde_json::from_str(minimal_json).unwrap();
        assert_eq!(restored.material, "glass");
        assert_eq!(restored.volume_ml, 200.0);
        // Missing fields should get their defaults from Default impl
        let defaults = StoredSettings::default();
        assert_eq!(restored.min_drinkable_c, defaults.min_drinkable_c);
        assert_eq!(restored.use_fahrenheit, defaults.use_fahrenheit);
        assert_eq!(restored.duration_minutes, defaults.duration_minutes);
        assert_eq!(restored.preheat_temp_c, defaults.preheat_temp_c);
        assert_eq!(restored.h_conv, defaults.h_conv);
        assert_eq!(restored.exposure_pct, defaults.exposure_pct);
    }

    #[test]
    fn test_serde_invalid_json_fails() {
        let bad_json = r#"{"volume_ml": "not a number"}"#;
        assert!(serde_json::from_str::<StoredSettings>(bad_json).is_err());
    }

    #[test]
    fn test_validate_clamps_preheat_and_duration() {
        let mut s = StoredSettings {
            preheat_temp_c: 200.0,   // over max 100
            duration_minutes: 500,    // over max 120
            ..Default::default()
        };
        s.validate();
        assert_eq!(s.preheat_temp_c, 110.0);
        assert_eq!(s.duration_minutes, 120);
    }

    #[test]
    fn test_validate_clamps_below_min() {
        let mut s = StoredSettings {
            preheat_temp_c: -20.0,   // under min -10
            duration_minutes: 1,     // under min 10
            h_conv: 0.5,            // under min 2
            exposure_pct: -10.0,    // under min 0
            ..Default::default()
        };
        s.validate();
        assert_eq!(s.preheat_temp_c, -10.0);
        assert_eq!(s.duration_minutes, 10);
        assert_eq!(s.h_conv, 2.0);
        assert_eq!(s.exposure_pct, 0.0);
    }

    #[test]
    fn test_validate_nan_preheat_and_exposure() {
        let mut s = StoredSettings {
            preheat_temp_c: f64::NAN,
            exposure_pct: f64::NAN,
            ..Default::default()
        };
        s.validate();
        let defaults = StoredSettings::default();
        assert_eq!(s.preheat_temp_c, defaults.preheat_temp_c);
        assert_eq!(s.exposure_pct, defaults.exposure_pct);
    }

    #[test]
    fn test_validate_drinkable_clamped_above_coffee() {
        let mut s = StoredSettings {
            coffee_temp_c: 90.0,
            room_temp_c: 22.0,
            min_drinkable_c: 70.0, // will clamp to 70 (within 50-70 range)
            ..Default::default()
        };
        s.validate();
        // min_drinkable must be < coffee_temp, so 70 → 85 (coffee - 5)
        assert!(s.min_drinkable_c < s.coffee_temp_c);
    }
}

use crate::physics::MugMaterial;
use crate::timing::{Debounce, SAVE_DEBOUNCE_MS, SCROLL_SAVE_DEBOUNCE_MS};
use leptos::*;
use serde::Serialize;
use std::cell::{Cell, RefCell};
use std::thread::LocalKey;
use std::time::Duration;
use wasm_bindgen::prelude::*;

const STORAGE_KEY: &str = "prewarm_settings";

/// Schema version written with every save. Bump when a stored field changes
/// meaning; unknown newer versions are still read best-effort and never wiped.
pub const SETTINGS_VERSION: u32 = 1;

const TAB_IDS: [&str; 3] = ["lumped", "heat_equation", "compare"];

#[derive(Serialize)]
pub struct StoredSettings {
    pub version: u32,
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
            version: SETTINGS_VERSION,
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
            // Numerical: the two heat-equation solvers agree on what preheating
            // buys, and the numerical one is the cheaper of the two.
            use_analytical_heat_eq: false,
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
        if !self.volume_ml.is_finite() {
            self.volume_ml = defaults.volume_ml;
        }
        if !self.wall_thickness_mm.is_finite() {
            self.wall_thickness_mm = defaults.wall_thickness_mm;
        }
        if !self.coffee_temp_c.is_finite() {
            self.coffee_temp_c = defaults.coffee_temp_c;
        }
        if !self.room_temp_c.is_finite() {
            self.room_temp_c = defaults.room_temp_c;
        }
        if !self.min_drinkable_c.is_finite() {
            self.min_drinkable_c = defaults.min_drinkable_c;
        }
        if !self.preheat_temp_c.is_finite() {
            self.preheat_temp_c = defaults.preheat_temp_c;
        }
        if !self.h_conv.is_finite() {
            self.h_conv = defaults.h_conv;
        }
        if !self.exposure_pct.is_finite() {
            self.exposure_pct = defaults.exposure_pct;
        }
        if !self.emissivity.is_finite() {
            self.emissivity = defaults.emissivity;
        }

        if !TAB_IDS.contains(&self.active_tab.as_str()) {
            self.active_tab = defaults.active_tab;
        }

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
    /// Held across `apply_external`: the batched `apply` runs the save effect
    /// before returning, and a tab must not write back settings it just read
    /// from storage (that would echo every change between tabs).
    static APPLYING_EXTERNAL: Cell<bool> = const { Cell::new(false) };
    /// Pending settings save; a burst of changes (a slider drag) becomes one write.
    static SAVE_DEBOUNCE: RefCell<Debounce> = RefCell::new(Debounce::new(SAVE_DEBOUNCE_MS));
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
                let s = Self::from_stored(&load_settings());
                cell.set(Some(s));
                s
            }
        })
    }

    /// Create a fresh set of signals holding the values in `saved`.
    fn from_stored(saved: &StoredSettings) -> Self {
        Self {
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
        }
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
            ..StoredSettings::default()
        }
    }

    /// Build `MugParameters` from the current signal values.
    ///
    /// Converts the stored lid opening `exposure_pct` (0–100) to `exposure` (0.0–1.0).
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
            exposure: self.exposure_pct.get() / 100.0,
        }
    }

    /// Apply `s` to the signals as one update: only fields whose value
    /// differs are set (`f64` compared exactly), all inside one `batch`, so
    /// dependent effects (solves, the save) run at most once per call.
    pub fn apply(&self, s: &StoredSettings) {
        batch(|| {
            set_if_changed(self.material, s.material_enum());
            set_if_changed(self.volume_ml, s.volume_ml);
            set_if_changed(self.wall_thickness_mm, s.wall_thickness_mm);
            set_if_changed(self.coffee_temp_c, s.coffee_temp_c);
            set_if_changed(self.room_temp_c, s.room_temp_c);
            set_if_changed(self.min_drinkable_c, s.min_drinkable_c);
            set_if_changed(self.use_fahrenheit, s.use_fahrenheit);
            set_if_changed(self.duration_minutes, s.duration_minutes);
            set_if_changed(self.preheat_temp_c, s.preheat_temp_c);
            set_if_changed(self.h_conv, s.h_conv);
            set_if_changed(self.exposure_pct, s.exposure_pct);
            set_if_changed(self.emissivity, s.emissivity);
            set_if_changed(self.active_tab, s.active_tab.clone());
            set_if_changed(self.use_analytical_lumped, s.use_analytical_lumped);
            set_if_changed(self.use_analytical_heat_eq, s.use_analytical_heat_eq);
            set_if_changed(self.use_near_wall, s.use_near_wall);
        });
    }
}

impl SharedSettings {
    /// Reset the mug and environment to defaults, keeping the user's
    /// preferences: temperature unit, active tab, solver and near-wall choices.
    pub fn reset_to_defaults(&self) {
        let current = untrack(|| self.to_stored());
        self.apply(&StoredSettings {
            use_fahrenheit: current.use_fahrenheit,
            active_tab: current.active_tab,
            use_analytical_lumped: current.use_analytical_lumped,
            use_analytical_heat_eq: current.use_analytical_heat_eq,
            use_near_wall: current.use_near_wall,
            ..StoredSettings::default()
        });
    }

    /// Change the material together with its typical emissivity, as one update.
    pub fn set_material(&self, material: MugMaterial) {
        batch(|| {
            self.material.set(material);
            self.emissivity.set(material.emissivity());
        });
    }
}

/// Set `signal` to `value` only when it differs, so unchanged fields notify no one.
fn set_if_changed<T: PartialEq + 'static>(signal: RwSignal<T>, value: T) {
    if signal.with_untracked(|current| *current != value) {
        signal.set(value);
    }
}

/// Returns true while signals are being updated from settings another tab stored.
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
            sync_unless_pending(|| parse_settings(&json));
        }
    }) as Box<dyn Fn(web_sys::StorageEvent)>);

    // Re-read localStorage on tab/window focus changes
    let sync_from_storage = Closure::wrap(Box::new(move || {
        sync_unless_pending(load_settings);
    }) as Box<dyn Fn()>);

    if let Some(window) = web_sys::window() {
        let _ =
            window.add_event_listener_with_callback("storage", storage_cb.as_ref().unchecked_ref());
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

/// Follow the stored settings, unless this tab has an unsaved change: then
/// that change is saved instead (other tabs follow it), so a sync never
/// overwrites the user's latest edit with older stored values.
fn sync_unless_pending(read: impl FnOnce() -> StoredSettings) {
    if !flush_pending_save() {
        apply_external(&read());
    }
}

/// Apply settings read from storage (written by another tab) into the shared
/// signals, suppressing the save effect to avoid a feedback loop.
fn apply_external(stored: &StoredSettings) {
    let settings = SharedSettings::instance();
    APPLYING_EXTERNAL.with(|c| c.set(true));
    settings.apply(stored);
    APPLYING_EXTERNAL.with(|c| c.set(false));
}

/// Save the current settings once `SAVE_DEBOUNCE_MS` pass without another call.
pub fn schedule_save() {
    schedule_debounced(&SAVE_DEBOUNCE, save_current_settings);
}

/// Save now if a scheduled save is pending; returns whether it saved.
pub fn flush_pending_save() -> bool {
    flush_debounced(&SAVE_DEBOUNCE, save_current_settings)
}

/// Record an input on `debounce` and arm a timer that runs `run` once the
/// debounce's quiet period has passed without another input.
fn schedule_debounced(debounce: &'static LocalKey<RefCell<Debounce>>, run: fn()) {
    let now = js_sys::Date::now();
    let deadline = debounce.with(|d| d.borrow_mut().on_input(now));
    arm_debounce_timer(debounce, run, deadline, (deadline - now).max(0.0));
}

/// Arm a timer for `deadline`. A timer can fire a little before its deadline
/// by `Date::now()`; if the input it was armed for is still the pending one,
/// it re-arms for the remaining time, or the save would be lost until the
/// next input. A later input has armed its own timer, so this one then stops.
fn arm_debounce_timer(
    debounce: &'static LocalKey<RefCell<Debounce>>,
    run: fn(),
    deadline: f64,
    delay_ms: f64,
) {
    set_timeout(
        move || {
            let now = js_sys::Date::now();
            if debounce.with(|d| d.borrow_mut().on_timer(now)) {
                run();
            } else if debounce.with(|d| d.borrow().deadline()) == Some(deadline) {
                // At least 1 ms, so a clock that has not yet ticked over to
                // the deadline does not spin zero-delay timers.
                arm_debounce_timer(debounce, run, deadline, (deadline - now).max(1.0));
            }
        },
        Duration::from_secs_f64(delay_ms / 1000.0),
    );
}

/// Run `run` now if `debounce` has a pending input; returns whether it ran.
fn flush_debounced(debounce: &'static LocalKey<RefCell<Debounce>>, run: fn()) -> bool {
    let pending = debounce.with(|d| d.borrow_mut().take_pending());
    if pending {
        run();
    }
    pending
}

/// Flush a pending save on `pagehide` (tab closed, navigation or reload, or
/// entering the back/forward cache), so the last change of a drag is never lost.
pub fn flush_pending_save_on_pagehide() {
    let closure = Closure::wrap(Box::new(move || {
        flush_pending_save();
    }) as Box<dyn Fn()>);
    if let Some(window) = web_sys::window() {
        let _ =
            window.add_event_listener_with_callback("pagehide", closure.as_ref().unchecked_ref());
    }
    closure.forget();
}

fn save_current_settings() {
    save_settings(&untrack(|| SharedSettings::instance().to_stored()));
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

/// Read one field from a parsed JSON object; a missing or wrongly-typed
/// field yields `None` so the caller falls back to that field's default alone.
fn field<T: serde::de::DeserializeOwned>(
    map: &serde_json::Map<String, serde_json::Value>,
    key: &str,
) -> Option<T> {
    map.get(key)
        .and_then(|v| serde_json::from_value(v.clone()).ok())
}

/// Parse stored settings JSON, recovering field by field.
///
/// Each known field is read independently, so one bad field defaults alone
/// instead of discarding the rest. Invalid JSON gives all defaults, and a
/// newer `version` is read best-effort. The result is validated.
pub fn parse_settings(json: &str) -> StoredSettings {
    let d = StoredSettings::default();
    let map = match serde_json::from_str::<serde_json::Value>(json) {
        Ok(serde_json::Value::Object(m)) => m,
        _ => return d,
    };
    let mut s = StoredSettings {
        version: field(&map, "version").unwrap_or(d.version),
        material: field(&map, "material").unwrap_or(d.material),
        volume_ml: field(&map, "volume_ml").unwrap_or(d.volume_ml),
        wall_thickness_mm: field(&map, "wall_thickness_mm").unwrap_or(d.wall_thickness_mm),
        coffee_temp_c: field(&map, "coffee_temp_c").unwrap_or(d.coffee_temp_c),
        room_temp_c: field(&map, "room_temp_c").unwrap_or(d.room_temp_c),
        min_drinkable_c: field(&map, "min_drinkable_c").unwrap_or(d.min_drinkable_c),
        use_fahrenheit: field(&map, "use_fahrenheit").unwrap_or(d.use_fahrenheit),
        duration_minutes: field(&map, "duration_minutes").unwrap_or(d.duration_minutes),
        preheat_temp_c: field(&map, "preheat_temp_c").unwrap_or(d.preheat_temp_c),
        h_conv: field(&map, "h_conv").unwrap_or(d.h_conv),
        exposure_pct: field(&map, "exposure_pct").unwrap_or(d.exposure_pct),
        emissivity: field(&map, "emissivity").unwrap_or(d.emissivity),
        active_tab: field(&map, "active_tab").unwrap_or(d.active_tab),
        use_analytical_lumped: field(&map, "use_analytical_lumped")
            .unwrap_or(d.use_analytical_lumped),
        use_analytical_heat_eq: field(&map, "use_analytical_heat_eq")
            .unwrap_or(d.use_analytical_heat_eq),
        use_near_wall: field(&map, "use_near_wall").unwrap_or(d.use_near_wall),
    };
    s.validate();
    s
}

/// Build the JSON to store: the existing stored object (so fields written by
/// a newer app version survive) with this app's fields overwritten by
/// `current`, and `version` kept at the higher of stored and current.
pub fn settings_json_for_save(existing: Option<&str>, current: &StoredSettings) -> String {
    let mut map = match existing.map(serde_json::from_str::<serde_json::Value>) {
        Some(Ok(serde_json::Value::Object(m))) => m,
        _ => serde_json::Map::new(),
    };
    let existing_version = map
        .get("version")
        .and_then(|v| v.as_u64())
        .and_then(|v| u32::try_from(v).ok())
        .unwrap_or(0);
    if let Ok(serde_json::Value::Object(fields)) = serde_json::to_value(current) {
        map.extend(fields);
    }
    map.insert(
        "version".to_string(),
        serde_json::Value::from(existing_version.max(SETTINGS_VERSION)),
    );
    serde_json::Value::Object(map).to_string()
}

/// Load user settings from localStorage.
///
/// Returns default settings if localStorage is unavailable (private browsing,
/// disabled, etc.) or nothing has been saved yet; stored data goes through
/// `parse_settings`, so corruption only costs the affected fields.
pub fn load_settings() -> StoredSettings {
    match get_storage() {
        Some(storage) => match storage.get_item(STORAGE_KEY) {
            Ok(Some(json)) => parse_settings(&json),
            Ok(None) => StoredSettings::default(),
            Err(_) => {
                log_warning("Prewarm: Could not read from localStorage (using defaults)");
                StoredSettings::default()
            }
        },
        // localStorage unavailable - likely private browsing mode (expected, no warning)
        None => StoredSettings::default(),
    }
}

/// Save user settings to localStorage, merging into whatever is stored.
///
/// Fails silently (with a console warning on write errors) if localStorage is
/// unavailable or full; the app keeps working without persistence.
pub fn save_settings(settings: &StoredSettings) {
    let storage = match get_storage() {
        Some(s) => s,
        None => return, // localStorage unavailable, silently skip (expected in private mode)
    };

    let existing = storage.get_item(STORAGE_KEY).ok().flatten();
    let json = settings_json_for_save(existing.as_deref(), settings);

    if storage.set_item(STORAGE_KEY, &json).is_err() {
        log_warning(
            "Prewarm: Could not save settings to localStorage. \
             You may be in private browsing mode or storage quota exceeded.",
        );
    }
}

/// Prefix of the per-path scroll key; the full key is [`scroll_key`].
const SCROLL_KEY_PREFIX: &str = "prewarm_scroll:";
/// Set just before the `popstate` reload workaround reloads the page: that
/// reload is a Back/forward, so the page restores its position after it.
const SCROLL_RESTORE_FLAG: &str = "prewarm_scroll_restore";
/// How long after load the saved position is applied, so the page has rendered.
const SCROLL_RESTORE_DELAY_MS: i32 = 50;

/// Session-storage key holding the scroll position of the page at `path`.
pub fn scroll_key(path: &str) -> String {
    format!("{SCROLL_KEY_PREFIX}{path}")
}

/// Whether a page load restores the saved position: only a Back/forward does
/// (a browser `back_forward` load, or the reload the `popstate` workaround
/// makes for one). A link, a typed URL or a plain reload opens at the top.
pub fn should_restore_scroll(navigation_type: Option<&str>, popstate_reload: bool) -> bool {
    popstate_reload || navigation_type == Some("back_forward")
}

/// Parse a stored scroll position; anything but a finite, non-negative number
/// is ignored.
pub fn parse_scroll_y(stored: &str) -> Option<f64> {
    stored
        .parse::<f64>()
        .ok()
        .filter(|y| y.is_finite() && *y >= 0.0)
}

/// Set up per-page scroll restoration - call once on app init, before mount.
///
/// Every page's position is kept in session storage under [`scroll_key`] of
/// its path. Scrolling saves it (debounced), leaving a page saves it (route
/// change, `popstate`, `pagehide`), and arriving by a link stores the top.
pub fn setup_scroll_restoration() {
    let window = window();
    // The browser restores one position per history entry, which the
    // popstate reload workaround defeats; the app restores per path itself.
    if let Ok(history) = window.history() {
        let _ = history.set_scroll_restoration(web_sys::ScrollRestoration::Manual);
    }

    let on_scroll = Closure::wrap(Box::new(record_scroll) as Box<dyn Fn()>);
    let _ = window.add_event_listener_with_callback("scroll", on_scroll.as_ref().unchecked_ref());
    on_scroll.forget();

    let on_pagehide = Closure::wrap(Box::new(|| {
        flush_scroll_save();
    }) as Box<dyn Fn()>);
    let _ =
        window.add_event_listener_with_callback("pagehide", on_pagehide.as_ref().unchecked_ref());
    on_pagehide.forget();
}

thread_local! {
    /// Pending scroll save: the debounce, and the path and position it saves.
    /// The path is taken when the scroll happens, so a save that lands after
    /// a navigation still goes to the page that was scrolled.
    static SCROLL_DEBOUNCE: RefCell<Debounce> =
        RefCell::new(Debounce::new(SCROLL_SAVE_DEBOUNCE_MS));
    static PENDING_SCROLL: RefCell<Option<(String, f64)>> = const { RefCell::new(None) };
    /// Set once the `popstate` reload is under way. Until the page unloads,
    /// the router may already show the destination (a shorter page clamps
    /// `scrollY`, firing `scroll`), and nothing seen then belongs to any page.
    static LEAVING: Cell<bool> = const { Cell::new(false) };
}

fn record_scroll() {
    if LEAVING.get() {
        return;
    }
    let window = window();
    let path = window.location().pathname().unwrap_or_default();
    let y = window.scroll_y().unwrap_or(0.0);
    PENDING_SCROLL.with_borrow_mut(|pending| *pending = Some((path, y)));
    schedule_debounced(&SCROLL_DEBOUNCE, save_pending_scroll);
}

fn save_pending_scroll() {
    if let Some((path, y)) = PENDING_SCROLL.with_borrow_mut(Option::take) {
        save_scroll(&path, y);
    }
}

/// Save a pending scroll position now (leaving the page it belongs to).
fn flush_scroll_save() {
    flush_debounced(&SCROLL_DEBOUNCE, save_pending_scroll);
}

fn save_scroll(path: &str, y: f64) {
    if let Some(storage) = session_storage() {
        let _ = storage.set_item(&scroll_key(path), &y.to_string());
    }
}

fn session_storage() -> Option<web_sys::Storage> {
    window().session_storage().ok().flatten()
}

/// Save the position of the page being left before the `popstate` reload
/// workaround reloads, stop recording until the page unloads, and mark that
/// reload as one that restores.
pub fn prepare_popstate_reload() {
    flush_scroll_save();
    LEAVING.set(true);
    if let Some(storage) = session_storage() {
        let _ = storage.set_item(SCROLL_RESTORE_FLAG, "1");
    }
}

/// A client-side route change (a link) to `path`: save the page being left,
/// and store the top for `path`, where the link opens it.
///
/// Skipped while a `popstate` reload is under way: the router also follows
/// that `popstate`, and the destination's saved position must survive the reload.
pub fn on_route_change(path: &str) {
    if LEAVING.get() {
        return;
    }
    flush_scroll_save();
    save_scroll(path, 0.0);
}

/// On load: restore the page's saved position if this load is a Back/forward,
/// else store the top for it (where it opens). Call once after mount.
pub fn restore_scroll_on_load() {
    let popstate_reload = session_storage()
        .and_then(|s| s.get_item(SCROLL_RESTORE_FLAG).ok().flatten())
        .is_some();
    if popstate_reload {
        if let Some(storage) = session_storage() {
            let _ = storage.remove_item(SCROLL_RESTORE_FLAG);
        }
    }
    let path = window().location().pathname().unwrap_or_default();
    if should_restore_scroll(navigation_type().as_deref(), popstate_reload) {
        set_timeout(
            move || restore_scroll(&path),
            Duration::from_millis(SCROLL_RESTORE_DELAY_MS as u64),
        );
    } else {
        save_scroll(&path, 0.0);
    }
}

fn restore_scroll(path: &str) {
    let saved = session_storage()
        .and_then(|s| s.get_item(&scroll_key(path)).ok().flatten())
        .and_then(|stored| parse_scroll_y(&stored));
    if let Some(y) = saved {
        window().scroll_to_with_x_and_y(0.0, y);
    }
}

/// `performance.getEntriesByType("navigation")[0].type` of this page load.
fn navigation_type() -> Option<String> {
    let entry = window()
        .performance()?
        .get_entries_by_type("navigation")
        .get(0);
    js_sys::Reflect::get(&entry, &"type".into())
        .ok()?
        .as_string()
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
            volume_ml: 1000.0,      // over max 500
            wall_thickness_mm: 0.1, // under min 1
            coffee_temp_c: 200.0,   // over max 100
            room_temp_c: -50.0,     // under min -10
            h_conv: 100.0,          // over max 50
            exposure_pct: 150.0,    // over max 100
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
            room_temp_c: 80.0, // equal to coffee — invalid
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
    fn parse_keeps_good_fields_when_one_is_bad() {
        let s = parse_settings(r#"{"volume_ml":"x","room_temp_c":18}"#);
        assert_eq!(s.volume_ml, StoredSettings::default().volume_ml);
        assert_eq!(s.room_temp_c, 18.0);
    }

    #[test]
    fn parse_invalid_json_gives_defaults() {
        let s = parse_settings("not json {");
        assert_eq!(s.volume_ml, StoredSettings::default().volume_ml);
        assert_eq!(s.version, SETTINGS_VERSION);
        let s = parse_settings("[1,2]");
        assert_eq!(s.active_tab, "lumped");
    }

    #[test]
    fn parse_validates_active_tab() {
        assert_eq!(
            parse_settings(r#"{"active_tab":"bogus"}"#).active_tab,
            "lumped"
        );
        assert_eq!(
            parse_settings(r#"{"active_tab":"compare"}"#).active_tab,
            "compare"
        );
    }

    #[test]
    fn parse_reads_newer_version_best_effort() {
        let s = parse_settings(r#"{"version":7,"volume_ml":250,"future_field":true}"#);
        assert_eq!(s.version, 7);
        assert_eq!(s.volume_ml, 250.0);
    }

    fn as_map(json: &str) -> serde_json::Map<String, serde_json::Value> {
        match serde_json::from_str(json).unwrap() {
            serde_json::Value::Object(m) => m,
            _ => panic!("not an object"),
        }
    }

    #[test]
    fn save_preserves_unknown_fields() {
        let current = StoredSettings {
            volume_ml: 200.0,
            ..Default::default()
        };
        let out = settings_json_for_save(
            Some(r#"{"version":7,"future_field":1,"volume_ml":100}"#),
            &current,
        );
        let m = as_map(&out);
        assert_eq!(m["future_field"], 1);
        assert_eq!(m["version"], 7);
        assert_eq!(m["volume_ml"], 200.0);
    }

    #[test]
    fn save_without_existing_writes_current_version() {
        for existing in [None, Some("garbage")] {
            let out = settings_json_for_save(existing, &StoredSettings::default());
            let m = as_map(&out);
            assert_eq!(m["version"], SETTINGS_VERSION);
            assert_eq!(m["material"], "ceramic");
        }
    }

    #[test]
    fn save_then_parse_roundtrips() {
        let current = StoredSettings {
            volume_ml: 321.0,
            active_tab: "compare".to_string(),
            ..Default::default()
        };
        let back = parse_settings(&settings_json_for_save(None, &current));
        assert_eq!(back.volume_ml, 321.0);
        assert_eq!(back.active_tab, "compare");
    }

    #[test]
    fn save_then_parse_roundtrips_every_field() {
        // Every field differs from its default, so a key that `parse_settings`
        // reads under a different name than serde writes falls back and fails.
        let current = StoredSettings {
            version: SETTINGS_VERSION + 1,
            material: "glass".to_string(),
            volume_ml: 321.0,
            wall_thickness_mm: 7.5,
            coffee_temp_c: 85.0,
            room_temp_c: 18.0,
            min_drinkable_c: 55.0,
            use_fahrenheit: false,
            duration_minutes: 45,
            preheat_temp_c: 95.0,
            h_conv: 12.0,
            exposure_pct: 35.0,
            emissivity: 0.5,
            active_tab: "heat_equation".to_string(),
            use_analytical_lumped: false,
            use_analytical_heat_eq: true,
            use_near_wall: true,
        };
        let defaults = serde_json::to_value(StoredSettings::default()).unwrap();
        let expected = serde_json::to_value(&current).unwrap();
        for (key, value) in expected.as_object().unwrap() {
            assert_ne!(&defaults[key], value, "{key} must differ from its default");
        }
        let back = parse_settings(&serde_json::to_string(&current).unwrap());
        assert_eq!(serde_json::to_value(&back).unwrap(), expected);
    }

    #[test]
    fn test_validate_clamps_preheat_and_duration() {
        let mut s = StoredSettings {
            preheat_temp_c: 200.0, // over max 100
            duration_minutes: 500, // over max 120
            ..Default::default()
        };
        s.validate();
        assert_eq!(s.preheat_temp_c, 110.0);
        assert_eq!(s.duration_minutes, 120);
    }

    #[test]
    fn test_validate_clamps_below_min() {
        let mut s = StoredSettings {
            preheat_temp_c: -20.0, // under min -10
            duration_minutes: 1,   // under min 10
            h_conv: 0.5,           // under min 2
            exposure_pct: -10.0,   // under min 0
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

    /// Settings signals plus an effect that reads all of them, counting its runs.
    fn settings_with_run_counter() -> (SharedSettings, std::rc::Rc<Cell<u32>>) {
        let settings = SharedSettings::from_stored(&StoredSettings::default());
        let runs = std::rc::Rc::new(Cell::new(0));
        let counter = runs.clone();
        create_effect(move |_| {
            settings.to_stored();
            counter.set(counter.get() + 1);
        });
        (settings, runs)
    }

    #[test]
    fn apply_skips_unchanged_fields() {
        let rt = create_runtime();
        let (settings, runs) = settings_with_run_counter();
        assert_eq!(runs.get(), 1);
        settings.apply(&StoredSettings::default());
        assert_eq!(runs.get(), 1);
        rt.dispose();
    }

    #[test]
    fn apply_changes_once_for_many_fields() {
        let rt = create_runtime();
        let (settings, runs) = settings_with_run_counter();
        settings.apply(&StoredSettings {
            material: "glass".to_string(),
            volume_ml: 250.0,
            room_temp_c: 18.0,
            use_fahrenheit: false,
            active_tab: "compare".to_string(),
            ..Default::default()
        });
        assert_eq!(runs.get(), 2);
        rt.dispose();
    }

    #[test]
    fn apply_sets_changed_values() {
        let rt = create_runtime();
        let settings = SharedSettings::from_stored(&StoredSettings::default());
        let target = StoredSettings {
            material: "stainless".to_string(),
            volume_ml: 250.0,
            wall_thickness_mm: 3.0,
            coffee_temp_c: 85.0,
            room_temp_c: 18.0,
            min_drinkable_c: 55.0,
            use_fahrenheit: false,
            duration_minutes: 45,
            preheat_temp_c: 70.0,
            h_conv: 10.0,
            exposure_pct: 50.0,
            emissivity: 0.5,
            active_tab: "compare".to_string(),
            use_analytical_lumped: false,
            use_analytical_heat_eq: false,
            use_near_wall: true,
            ..Default::default()
        };
        settings.apply(&target);
        let expected = serde_json::to_value(&target).unwrap();
        assert_eq!(
            serde_json::to_value(settings.to_stored()).unwrap(),
            expected
        );
        rt.dispose();
    }

    #[test]
    fn reset_restores_defaults_but_keeps_preferences() {
        let rt = create_runtime();
        let prefs = StoredSettings {
            volume_ml: 250.0,
            material: "glass".to_string(),
            use_fahrenheit: false,
            active_tab: "compare".to_string(),
            use_analytical_lumped: false,
            use_analytical_heat_eq: false,
            use_near_wall: true,
            ..Default::default()
        };
        let settings = SharedSettings::from_stored(&prefs);
        settings.reset_to_defaults();
        let after = settings.to_stored();
        assert_eq!(after.volume_ml, StoredSettings::default().volume_ml);
        assert_eq!(after.material, "ceramic");
        assert!(!after.use_fahrenheit);
        assert_eq!(after.active_tab, "compare");
        assert!(!after.use_analytical_lumped);
        assert!(!after.use_analytical_heat_eq);
        assert!(after.use_near_wall);
        rt.dispose();
    }

    #[test]
    fn set_material_updates_emissivity_in_one_run() {
        let rt = create_runtime();
        let (settings, runs) = settings_with_run_counter();
        settings.set_material(MugMaterial::StainlessSteel);
        assert_eq!(runs.get(), 2);
        assert_eq!(
            settings.emissivity.get_untracked(),
            MugMaterial::StainlessSteel.emissivity()
        );
        rt.dispose();
    }

    #[test]
    fn scroll_key_is_per_path() {
        assert_eq!(scroll_key("/"), "prewarm_scroll:/");
        assert_eq!(scroll_key("/theory"), "prewarm_scroll:/theory");
        assert_ne!(scroll_key("/theory"), scroll_key("/theory/lumped"));
    }

    #[test]
    fn restores_only_on_back_forward() {
        assert!(should_restore_scroll(Some("back_forward"), false));
        assert!(should_restore_scroll(Some("reload"), true));
        assert!(!should_restore_scroll(Some("reload"), false));
        assert!(!should_restore_scroll(Some("navigate"), false));
        assert!(!should_restore_scroll(None, false));
    }

    #[test]
    fn parse_scroll_y_accepts_only_finite_non_negative() {
        assert_eq!(parse_scroll_y("1200"), Some(1200.0));
        assert_eq!(parse_scroll_y("0"), Some(0.0));
        assert_eq!(parse_scroll_y("12.5"), Some(12.5));
        assert_eq!(parse_scroll_y("-3"), None);
        assert_eq!(parse_scroll_y("NaN"), None);
        assert_eq!(parse_scroll_y("inf"), None);
        assert_eq!(parse_scroll_y("abc"), None);
    }
}

use leptos::*;
use std::cell::Cell;

/// Settings of the heat-equation and Compare panels' own controls and views.
///
/// A singleton like `SharedSettings`, so the controls and the panel that
/// solves with them can live in different parts of the tree (the overlay
/// renders its controls outside the panel it swaps on a tab change). Not
/// persisted: a reload starts from the defaults.
#[derive(Clone, Copy)]
pub struct PanelSettings {
    /// Output frames per heat-equation solve.
    pub num_snapshots: RwSignal<usize>,
    /// Eigenfunction series terms (analytical heat equation).
    pub num_modes: RwSignal<usize>,
    /// Simulation domain extent as a multiple of the mug's outer radius.
    pub r_max_mult: RwSignal<f64>,
    /// Compare panel scenario: preheated (true) or cold mug.
    pub compare_use_preheated: RwSignal<bool>,
    /// Heat-equation profile chart: index of the snapshot shown.
    pub selected_time_idx: RwSignal<usize>,
}

const DEFAULT_NUM_SNAPSHOTS: usize = 100;
const DEFAULT_NUM_MODES: usize = 50;
const DEFAULT_R_MAX_MULT: f64 = 10.0;
const DEFAULT_COMPARE_USE_PREHEATED: bool = true;
const DEFAULT_SELECTED_TIME_IDX: usize = 2;

thread_local! {
    static PANEL_SETTINGS: Cell<Option<PanelSettings>> = const { Cell::new(None) };
}

impl PanelSettings {
    /// Get or create the singleton. The app creates it at its root (as it
    /// does `SharedSettings`), so the signals outlive every route.
    pub fn instance() -> Self {
        PANEL_SETTINGS.with(|cell| {
            if let Some(s) = cell.get() {
                s
            } else {
                let s = Self::new();
                cell.set(Some(s));
                s
            }
        })
    }

    fn new() -> Self {
        Self {
            num_snapshots: create_rw_signal(DEFAULT_NUM_SNAPSHOTS),
            num_modes: create_rw_signal(DEFAULT_NUM_MODES),
            r_max_mult: create_rw_signal(DEFAULT_R_MAX_MULT),
            compare_use_preheated: create_rw_signal(DEFAULT_COMPARE_USE_PREHEATED),
            selected_time_idx: create_rw_signal(DEFAULT_SELECTED_TIME_IDX),
        }
    }

    /// Restore the defaults as one update.
    pub fn reset(&self) {
        batch(|| {
            self.num_snapshots.set(DEFAULT_NUM_SNAPSHOTS);
            self.num_modes.set(DEFAULT_NUM_MODES);
            self.r_max_mult.set(DEFAULT_R_MAX_MULT);
            self.compare_use_preheated
                .set(DEFAULT_COMPARE_USE_PREHEATED);
            self.selected_time_idx.set(DEFAULT_SELECTED_TIME_IDX);
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::rc::Rc;

    #[test]
    fn reset_restores_defaults_in_one_update() {
        let rt = create_runtime();
        let s = PanelSettings::new();
        s.num_snapshots.set(20);
        s.num_modes.set(3);
        s.r_max_mult.set(1.5);
        s.compare_use_preheated.set(false);
        s.selected_time_idx.set(7);

        let runs = Rc::new(Cell::new(0));
        let counter = runs.clone();
        create_effect(move |_| {
            s.num_snapshots.get();
            s.num_modes.get();
            s.r_max_mult.get();
            s.compare_use_preheated.get();
            s.selected_time_idx.get();
            counter.set(counter.get() + 1);
        });
        assert_eq!(runs.get(), 1);

        s.reset();
        assert_eq!(runs.get(), 2);
        assert_eq!(s.num_snapshots.get_untracked(), 100);
        assert_eq!(s.num_modes.get_untracked(), 50);
        assert_eq!(s.r_max_mult.get_untracked(), 10.0);
        assert!(s.compare_use_preheated.get_untracked());
        assert_eq!(s.selected_time_idx.get_untracked(), 2);
        rt.dispose();
    }
}

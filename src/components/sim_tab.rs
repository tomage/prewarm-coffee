use crate::components::{
    CompareExtraControls, ComparePanel, ControlsLayout, HeatEqExtraControls, HeatEquationPanel,
    LumpedPanel, PanelSettings,
};
use crate::storage::SharedSettings;
use leptos::*;

/// The simulator's tabs, each a model panel with its own extra controls.
/// `SharedSettings::active_tab` stores the tab's key.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SimTab {
    Lumped,
    HeatEquation,
    Compare,
}

impl SimTab {
    /// The tab stored under `key`; an unknown key is the lumped tab.
    pub fn from_key(key: &str) -> Self {
        match key {
            "heat_equation" => Self::HeatEquation,
            "compare" => Self::Compare,
            _ => Self::Lumped,
        }
    }

    /// The tab's model panel. `with_controls = false` leaves the controls to
    /// the caller.
    pub fn panel(self, with_controls: bool) -> View {
        match self {
            Self::Lumped => view! { <LumpedPanel with_controls=with_controls/> }.into_view(),
            Self::HeatEquation => {
                view! { <HeatEquationPanel with_controls=with_controls/> }.into_view()
            }
            Self::Compare => view! { <ComparePanel with_controls=with_controls/> }.into_view(),
        }
    }

    /// The tab's own controls, shown inside `MugControls` after the shared ones.
    pub fn extra_controls(self, layout: ControlsLayout) -> View {
        match self {
            Self::Lumped => ().into_view(),
            Self::HeatEquation => view! { <HeatEqExtraControls layout=layout/> }.into_view(),
            Self::Compare => view! { <CompareExtraControls layout=layout/> }.into_view(),
        }
    }

    /// Restore the tab's own settings to their defaults: what Reset does
    /// beyond the shared reset, on desktop and in the overlay alike.
    pub fn reset_own_settings(self) {
        match self {
            Self::Lumped => {}
            Self::HeatEquation => batch(|| {
                PanelSettings::instance().reset();
                SharedSettings::instance().use_near_wall.set(false);
            }),
            Self::Compare => PanelSettings::instance().reset(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn from_key_reads_each_stored_key_and_defaults_to_lumped() {
        assert_eq!(SimTab::from_key("lumped"), SimTab::Lumped);
        assert_eq!(SimTab::from_key("heat_equation"), SimTab::HeatEquation);
        assert_eq!(SimTab::from_key("compare"), SimTab::Compare);
        assert_eq!(SimTab::from_key("bogus"), SimTab::Lumped);
    }
}

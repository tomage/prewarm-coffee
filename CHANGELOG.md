# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Version 5] - 2026-10-29

### Added

- **Quality gates.** The toolchain is pinned (`rust-toolchain.toml`). `make lint` (`cargo fmt --check` and clippy with `-D warnings`, host and wasm32) and `make test-ci` (lint plus unit tests over an 85% line-coverage floor) are what CI gates on, and a pre-commit hook (`make install-git-hooks`) runs `make lint`.
- **Opt-in checks outside the default test run.** `make bench` runs the timing-sensitive perf checks against baselines, now including the analytical heat-equation solver in the browser (`?solver=analytical`). `make grid` reports how far the four solvers disagree over 576 mugs and conditions.
- **Tests that pin the physics:** energy conservation, the home page's stated numbers against both models, and analytical-vs-numerical agreement.

### Fixed

- **Physics corrections** in all four solvers (from the 2026-10 project review):
  - The heat-equation numerical solver conserves energy: interface nodes take heat capacity from both materials (it over-counted the mug's heat capacity), and grid segments carry their own material's conductivity.
  - No more step-to-step ringing at the mug wall: the solver takes 4 backward-Euler steps before switching to Crank–Nicolson.
  - The heat-equation analytical solver's lid-opening loss had the wrong sign; it now cools the coffee, and the post-hoc exposure correction is gone.
  - Exposure means lid openness everywhere (0% sealed lid, 100% open mug with no lid; the mug always exists) and results vary smoothly across the slider. The 100% "no mug" special case is gone.
  - One evaporation model, shared by all solvers, with no cutoff: coffee keeps evaporating at room temperature and can settle a few degrees below it.
  - The lumped coffee-to-mug conductance comes from the wall geometry (a water-side film of 500 W/(m²·K) in series with half the wall) instead of `2·k`. Steel mugs couple less strongly and plastic ones more.
  - Each physical constant and material property is defined once. Stainless steel's properties come from the one material table in every solver (conductivity 16 W/(m·K), density 8000 kg/m³).
  - When the coffee never drops below the drinkable threshold, the app says so instead of reporting the last time step, and the extra drinkable time reads "at least +N min".
  - The lumped tab's crossings and averages cover the displayed duration, as the chart does, not the full 120-minute run behind it.
- **Theory pages** describe the code as it is (boundary conditions, flux balance at the interfaces, time-step rule, single-region coefficients, evaporation and radiation magnitudes, coffee-to-mug conductance), and each lists its model's assumptions, including the measured gap between the two heat-equation solvers. The work-in-progress notice is gone. The Compare tab notes that the lumped and heat-equation models are different idealisations.
- **Sliders stay live without stalling** (from the 2026-10 project review). Dragging a slider re-solves the heat-equation and Compare charts at most every 100 ms and always once more with the final value, so the chart follows the drag and rests on the latest inputs. Each tab re-solves only the solves whose inputs changed; the lumped tab re-solves on every input.
- **Settings sync between tabs without storms.** A drag is saved once, 250 ms after it settles (and on leaving the page, so the last change is kept). Another tab applies the stored settings as one update and only where they differ, and a tab with an unsaved change keeps it rather than taking older stored values.
- **Settings survive bad data.** Stored settings carry a version and are read field by field: a bad field falls back to its default alone, settings from a newer version are read as far as they are understood and never wiped, and saving keeps fields this version does not know. Nothing is written until a setting changes (a `?tab=` link selecting a different tab than the stored one counts as a change).
- **One simulator on mobile.** Small screens mount only the entry card and its overlay, not a hidden copy of the desktop simulator, so each change solves once and element ids are unique. Switching tabs in the overlay swaps only the chart panel and the tab's own controls; open accordion groups and the controls' scroll position stay. Desktop, overlay and theory pages share one set of mug controls.
- **The mobile simulator overlay is an accessible dialog.** It has `role="dialog"` with a label, takes focus when it opens, keeps Tab inside, and closes on Escape, Done or the Back button. Back closes it without reloading the page, and focus returns to "Open Simulator".
- **Scroll position per page.** Back and forward restore the position of the page they land on; following a link opens the page at the top.
- **Math renders on theory deep links.** KaTeX loads from absolute paths, so a theory page opened directly (or reloaded) renders its equations.

### Changed

- **The heat-equation analytical solver follows the coffee as it cools.** It re-linearises the open-top loss and the mug's radiation over ten time segments (midpoint rule), solving each exactly and starting each from the previous one's temperature profile. It used to linearise them once and overstate what preheating buys by about half a minute of drinkable time; it now agrees with the numerical solver to within 0.1 minute and 0.2 °C for ceramic mugs at the default lid opening, with crossings 0.2–0.4 minutes early. The note about the gap is gone. Each analytical curve takes about 1.7× as long to compute as before at the default settings (about 2× at 100 series terms and 120 minutes).

- **The heat-equation tab defaults to the numerical solver**, which is both faithful and the cheaper of the two. Visitors who already chose a solver keep their choice.

- **Home states the conclusion.** A new "What It Adds Up To" section gives what preheating buys (about 3½ more minutes above drinking temperature and 5–8 °C warmer for a typical ceramic mug) and how wall thickness changes it, and the intro gives the equalisation drop (about 13 °C cold vs 2 °C preheated). A test checks these numbers against both models.

- **Cooling numbers changed.** Default settings (180 ml ceramic, 5 mm wall, 90 °C coffee, 80 °C preheat, 22 °C room, 60 °C drinkable, h_conv = 20 W/(m²·K), lid opening 20%, 30 min), numerical solvers:

  | | Cold-mug penalty (before → after) | Cold mug reaches 60 °C | Preheated mug reaches 60 °C |
  |---|---|---|---|
  | Lumped | 6.27 → 6.94 °C | 10.0 → 9.4 min | 13.4 → 12.8 min |
  | Heat equation | 8.01 → 7.43 °C | 8.9 → 9.1 min | 12.8 → 12.6 min |

  The penalty is how much warmer, on average, coffee in a preheated mug is while coffee in a cold mug stays above the drinkable threshold.

## [Version 4] - 2026-03-01

### Added

- **Mobile simulator overlay**: Fullscreen interactive mode on mobile (≤700px) with chart pinned at top and scrollable accordion-grouped controls below
- New `AccordionGroup` component for collapsible control sections with `aria-expanded` accessibility support
- New `SimulatorOverlay` component with tab switching, body scroll lock, and independent chart/controls scrolling
- "Open Simulator" entry card on mobile replaces the inline simulator layout
- **Offline KaTeX**: Bundled KaTeX JS, CSS, and 9 woff2 fonts locally, removing CDN dependency for full offline support

### Changed

- Mobile layout: controls reorganized into 4 accordion groups (Mug Properties, Temperatures, Simulation, Advanced) instead of flat list
- Desktop layout unchanged — all mobile changes gated behind ≤700px media query
- Overlay chart area scrolls independently for panels with multiple charts (Heat Equation profile, Compare table)

### Dependencies

- Added `web-sys` features: `DomTokenList`, `HtmlElement` (for body scroll lock)

## [Version 3] - 2026-02-28

### Added

- ~900 lines of cross-solver tests validating agreement between analytical, numerical, and lumped approaches
- Expanded performance benchmarks in `perf_tests.rs`
- Minor E2E smoke test additions

### Changed

- Expanded and refactored lumped analytical solver (`lumped_analytical.rs`)
- Expanded and refactored lumped numerical solver (`lumped_numerical.rs`)
- Tuned cooling rates to better match real-world observations
- Renamed terminology throughout UI and theory pages: "liquid" to "coffee", "shell" to "mug" for clarity
- Refined heat equation analytical solver
- Updated comparison metrics display
- Added material properties for new solver configurations

## [Version 2] - 2026-02-04

### Added

- **Heat equation solver**: Full three-region eigenfunction expansion analytical solver (`heat_equation_analytical.rs`, 1,571 lines)
- **Numerical solver**: Finite-difference numerical heat equation solver (`heat_equation_numerical.rs`, 931 lines)
- **Client-side routing** via `leptos_router` with multi-page navigation
- **Pages**: About, Theory (overview, lumped model, heat equation), and 404 page
- **Comparison mode**: Side-by-side preheated vs non-preheated panel (`compare_panel.rs`)
- **New components**: cooling chart, profile chart, solver toggle, comparison metrics, simulation results, and navigation bar
- **Cross-solver tests**: 693 lines validating analytical vs numerical consistency
- **Performance tests**: 183 lines of benchmarks with baseline tracking
- **E2E test suite**: Playwright-based tests for navigation, smoke testing, and performance
- Lumped analytical model (`lumped_analytical.rs`)
- Temperature metrics calculations module (`metrics.rs`)

### Changed

- Restructured `App` component to use router and new panel-based layout
- Renamed `cooling.rs` to `lumped_numerical.rs` to clarify its role
- Removed `temperature_display.rs`, replaced by new chart and panel components
- Expanded `storage.rs` with cross-tab sync via `StorageEvent`
- Bumped dev `opt-level` from 1 to 2 for faster development builds
- Overhauled `styles.css` with 1,700+ lines of new styles

### Dependencies

- Added `leptos_router` 0.6 for client-side routing
- Added `web-sys` features: `StorageEvent`, `Element`, `MouseEvent`, `DomRect`

## [Version 1] - 2026-01-26

### Added

- Initial release of the Prewarm Coffee interactive website
- **Physics engine**: Lumped capacitance heat transfer model (`cooling.rs`, 534 lines) with mug material properties
- **Interactive UI**: Mug controls with sliders for material, thickness, volume, and other properties
- **Temperature visualization**: Chart showing cooling curves over time
- **Mobile support**: Responsive layout with mobile summary component
- **State persistence**: LocalStorage-based settings and scroll position restoration
- **Analytics**: Privacy-focused analytics integration
- **Haptics**: Haptic feedback support for mobile devices
- **Deployment**: GitHub Actions CI/CD pipeline with custom domain support
- Built with Rust + Leptos 0.6 (CSR) compiled to WASM via Trunk
- AGPL-3.0 license

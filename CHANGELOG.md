# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

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

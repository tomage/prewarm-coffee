//! Exposure is lid openness in every solver: coffee temperature after a fixed
//! time falls monotonically and smoothly as the lid opens from sealed (0) to
//! no lid (1). A model switch at either end shows up as one outsized step.

use super::heat_equation_analytical::relinearised_liquid_average;
use super::*;

const VOLUME_ML: f64 = 350.0;
const WALL_MM: f64 = 5.0;
const COFFEE_C: f64 = 90.0;
const ROOM_C: f64 = 22.0;
const H_CONV: f64 = 20.0;
const DURATION_MIN: usize = 30;
/// Sphere domain extent as a multiple of the outer radius, as the panels use.
const R_MAX_MULT: f64 = 10.0;
/// Eigenmodes for the analytical heat equation, as the panels use.
const NUM_MODES: usize = 50;

/// Coffee temperature at 30 min for exposure 0.00, 0.01, …, 1.00.
fn sweep(coffee_at_30_min: impl Fn(f64) -> f64) -> Vec<f64> {
    (0..=100)
        .map(|k| coffee_at_30_min(k as f64 / 100.0))
        .collect()
}

/// Each step is non-increasing, and no single step's drop exceeds three
/// times the mean step drop.
fn assert_monotone_and_smooth(solver: &str, temps: &[f64]) {
    let n = temps.len() - 1;
    let mean_drop = (temps[0] - temps[n]) / n as f64;
    assert!(
        mean_drop > 0.0,
        "{solver}: opening the lid should cool the coffee (sealed {:.3} °C, open {:.3} °C)",
        temps[0],
        temps[n]
    );
    for k in 0..n {
        let drop = temps[k] - temps[k + 1];
        assert!(
            drop >= -1e-9,
            "{solver}: warmer at exposure {:.2} ({:.4} °C) than at {:.2} ({:.4} °C)",
            (k + 1) as f64 / 100.0,
            temps[k + 1],
            k as f64 / 100.0,
            temps[k]
        );
        assert!(
            drop <= 3.0 * mean_drop,
            "{solver}: step {:.2}→{:.2} drops {drop:.4} K, more than 3× the mean step {mean_drop:.4} K",
            k as f64 / 100.0,
            (k + 1) as f64 / 100.0
        );
    }
}

fn mug(exposure: f64) -> MugParameters {
    MugParameters {
        material: MugMaterial::Ceramic,
        volume_ml: VOLUME_ML,
        wall_thickness_mm: WALL_MM,
        coffee_temp_c: COFFEE_C,
        room_temp_c: ROOM_C,
        preheated: false,
        preheat_temp_c: ROOM_C,
        h_conv: H_CONV,
        emissivity: MugMaterial::Ceramic.emissivity(),
        exposure,
    }
}

/// Inner and outer sphere radii in m for the same volume and wall.
fn sphere_radii() -> (f64, f64) {
    let inner = (3.0 * VOLUME_ML * 1e-6 / (4.0 * std::f64::consts::PI)).cbrt();
    (inner, inner + WALL_MM * 1e-3)
}

fn shell_material() -> HeatMaterial {
    HeatMaterial::from_mug(MugMaterial::Ceramic, MugMaterial::Ceramic.emissivity())
}

#[test]
fn exposure_response_is_monotone_and_smooth_lumped_numerical() {
    let temps = sweep(|e| {
        let data = simulate_lumped_numerical(&mug(e), DURATION_MIN);
        data.last().unwrap().coffee_temp_c
    });
    assert_monotone_and_smooth("lumped numerical", &temps);
}

#[test]
fn exposure_response_is_monotone_and_smooth_lumped_analytical() {
    let temps = sweep(|e| {
        let data = simulate_lumped_analytical(&mug(e), DURATION_MIN);
        data.last().unwrap().coffee_temp_c
    });
    assert_monotone_and_smooth("lumped analytical", &temps);
}

#[test]
fn exposure_response_is_monotone_and_smooth_heat_eq_numerical() {
    let (inner, outer) = sphere_radii();
    let temps = sweep(|e| {
        let params = SphereHeatParams {
            outer_radius_m: outer,
            inner_radius_m: inner,
            shell_temp_c: ROOM_C,
            liquid_temp_c: COFFEE_C,
            room_temp_c: ROOM_C,
            shell_material: shell_material(),
            liquid_material: HeatMaterial::water(),
            medium_material: HeatMaterial::air(),
            exposure: e,
            convection_coeff: H_CONV,
            r_max_multiplier: R_MAX_MULT,
            outer_boundary: OuterBoundary::Ambient,
        };
        let snaps = simulate_sphere_heat_numerical(&params, DURATION_MIN as f64 * 60.0, 2);
        liquid_average_temperature(&snaps).last().unwrap().1
    });
    assert_monotone_and_smooth("heat-eq numerical", &temps);
}

#[test]
fn exposure_response_is_monotone_and_smooth_heat_eq_analytical() {
    let (inner, outer) = sphere_radii();
    let shell = shell_material();
    let temps = sweep(|e| {
        let mug = ThreeRegionMug {
            r_inner: inner,
            r_outer: outer,
            r_max: outer * R_MAX_MULT,
            shell,
            liquid_temp_c: COFFEE_C,
            shell_temp_c: ROOM_C,
            room_temp_c: ROOM_C,
            h_conv: H_CONV,
            exposure: e,
            num_modes: NUM_MODES,
        };
        let avg =
            relinearised_liquid_average(&mug, ANALYTICAL_SEGMENTS, DURATION_MIN as f64 * 60.0, 2);
        avg.last().unwrap().1
    });
    assert_monotone_and_smooth("heat-eq analytical", &temps);
}

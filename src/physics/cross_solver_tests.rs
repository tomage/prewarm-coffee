//! Cross-solver comparison tests.
//!
//! Verifies that analytical and numerical solvers produce consistent
//! temperature curves when given identical physical parameters.
//! Covers both the lumped capacitance model and the three-region
//! eigenfunction expansion / finite-difference heat equation.

use super::*;
use heat_equation_analytical::three_region_liquid_average_closed_form;

// ─── Lumped capacitance model ──────────────────────────────────────

#[test]
fn test_lumped_analytical_vs_numerical_cold() {
    let params = MugParameters {
        material: MugMaterial::Ceramic,
        volume_ml: 350.0,
        wall_thickness_mm: 5.0,
        coffee_temp_c: 95.0,
        room_temp_c: 22.0,
        preheated: false,
        preheat_temp_c: 90.0,
        h_conv: 10.0,
        emissivity: MugMaterial::Ceramic.emissivity(),
        exposure_pct: 0.2,
    };

    let numerical = simulate_lumped_numerical(&params, 60);
    let analytical = simulate_lumped_analytical(&params, 60);

    // Both should produce similar number of data points
    assert!(
        (numerical.len() as i64 - analytical.len() as i64).unsigned_abs() <= 2,
        "Point counts should be similar: numerical={}, analytical={}",
        numerical.len(),
        analytical.len(),
    );

    // With piecewise linearization, the analytical solver closely tracks the numerical one.
    let count = numerical.len().min(analytical.len());
    let mut max_coffee_diff = 0.0_f64;
    let mut max_mug_diff = 0.0_f64;
    for i in 0..count {
        let dt = (numerical[i].time_minutes - analytical[i].time_minutes).abs();
        assert!(dt < 0.01, "Time mismatch at index {}: {:.3} vs {:.3}", i, numerical[i].time_minutes, analytical[i].time_minutes);

        let coffee_diff = (numerical[i].coffee_temp_c - analytical[i].coffee_temp_c).abs();
        max_coffee_diff = max_coffee_diff.max(coffee_diff);
        assert!(
            coffee_diff < 0.5,
            "Coffee temp mismatch at t={:.1}: numerical={:.2}, analytical={:.2}, diff={:.3}",
            numerical[i].time_minutes, numerical[i].coffee_temp_c, analytical[i].coffee_temp_c, coffee_diff
        );

        let mug_diff = (numerical[i].mug_temp_c - analytical[i].mug_temp_c).abs();
        max_mug_diff = max_mug_diff.max(mug_diff);
        assert!(
            mug_diff < 0.5,
            "Mug temp mismatch at t={:.1}: numerical={:.2}, analytical={:.2}, diff={:.3}",
            numerical[i].time_minutes, numerical[i].mug_temp_c, analytical[i].mug_temp_c, mug_diff
        );
    }
    println!("Cold mug: max coffee diff={:.4}°C, max mug diff={:.4}°C", max_coffee_diff, max_mug_diff);
}

#[test]
fn test_lumped_analytical_vs_numerical_preheated() {
    let params = MugParameters {
        material: MugMaterial::Ceramic,
        volume_ml: 350.0,
        wall_thickness_mm: 5.0,
        coffee_temp_c: 95.0,
        room_temp_c: 22.0,
        preheated: true,
        preheat_temp_c: 90.0,
        h_conv: 10.0,
        emissivity: MugMaterial::Ceramic.emissivity(),
        exposure_pct: 0.2,
    };

    let numerical = simulate_lumped_numerical(&params, 30);
    let analytical = simulate_lumped_analytical(&params, 30);

    // Piecewise linearization keeps divergence under 1°C.
    let count = numerical.len().min(analytical.len());
    let mut max_diff = 0.0_f64;
    for i in 0..count {
        let coffee_diff = (numerical[i].coffee_temp_c - analytical[i].coffee_temp_c).abs();
        max_diff = max_diff.max(coffee_diff);
        assert!(
            coffee_diff < 0.5,
            "Preheated coffee mismatch at t={:.1}: diff={:.3}",
            numerical[i].time_minutes, coffee_diff
        );
    }
    println!("Preheated mug: max coffee diff={:.4}°C", max_diff);
}

// ─── Three-region heat equation ────────────────────────────────────

/// Build matching params for both solvers from a common physical setup.
struct CommonSetup {
    analytical: ThreeRegionParams,
    numerical: SphereHeatParams,
}

fn setup(liquid_temp: f64, shell_temp: f64, room_temp: f64, h_conv: f64) -> CommonSetup {
    let water = HeatMaterial::water();
    let ceramic = HeatMaterial::ceramic();
    let air = HeatMaterial::air();

    let volume_ml = 180.0;
    let inner_r_cm = (3.0 * volume_ml / (4.0 * std::f64::consts::PI)).cbrt();
    let wall_mm = 5.0;
    let outer_r_cm = inner_r_cm + wall_mm / 10.0;
    let r_max_mult = 10.0;

    let r_inner = inner_r_cm / 100.0;
    let r_outer = outer_r_cm / 100.0;
    let r_max = outer_r_cm * r_max_mult / 100.0;

    CommonSetup {
        analytical: ThreeRegionParams {
            r_inner,
            r_outer,
            r_max,
            alpha_1: water.diffusivity(),
            alpha_2: ceramic.diffusivity(),
            alpha_3: air.diffusivity(),
            k_1: water.conductivity,
            k_2: ceramic.conductivity,
            k_3: air.conductivity,
            t_liquid_init: liquid_temp,
            t_shell_init: shell_temp,
            t_ambient: room_temp,
            num_modes: 50,
            h_conv,
        },
        numerical: SphereHeatParams {
            outer_radius_m: r_outer,
            inner_radius_m: r_inner,
            shell_temp_c: shell_temp,
            liquid_temp_c: liquid_temp,
            room_temp_c: room_temp,
            shell_material: ceramic,
            liquid_material: water,
            medium_material: air,
            evaporation_factor: 0.0,
            convection_coeff: h_conv,
            r_max_multiplier: r_max_mult,
        },
    }
}

/// Compare liquid average temperature curves from both solvers.
/// Returns (max_diff, rms_diff) in °C.
fn compare_liquid_avg(
    setup: &CommonSetup,
    duration_s: f64,
    num_snapshots: usize,
) -> (f64, f64) {
    // Analytical — closed-form volume average eliminates all radial discretization
    let modes = compute_eigenmodes(&setup.analytical);
    let analytical_avg = three_region_liquid_average_closed_form(
        &setup.analytical,
        &modes,
        duration_s,
        num_snapshots,
    );

    // Numerical
    let numerical_snaps =
        simulate_sphere_heat_numerical(&setup.numerical, duration_s, num_snapshots);
    let numerical_avg = liquid_average_temperature(&numerical_snaps);

    assert_eq!(
        analytical_avg.len(),
        numerical_avg.len(),
        "Snapshot count mismatch"
    );

    let mut max_diff = 0.0_f64;
    let mut sum_sq = 0.0_f64;
    let n = analytical_avg.len();

    for i in 0..n {
        let (_t_a, temp_a) = analytical_avg[i];
        let (_t_n, temp_n) = numerical_avg[i];
        let diff = (temp_a - temp_n).abs();
        max_diff = max_diff.max(diff);
        sum_sq += diff * diff;

        // Print every 10th point for diagnosis
        if i % (n / 10).max(1) == 0 || diff > 3.0 {
            let t_min = _t_a / 60.0;
            println!(
                "  t={:5.1}min: analytical={:.2}°C, numerical={:.2}°C, diff={:.2}°C",
                t_min, temp_a, temp_n, diff,
            );
        }
    }

    let rms = (sum_sq / n as f64).sqrt();
    println!("  max_diff={:.2}°C, rms={:.2}°C", max_diff, rms);
    (max_diff, rms)
}

#[test]
fn test_heat_eq_analytical_vs_numerical_preheated() {
    // Preheated ceramic mug: shell starts at 80°C
    // Note: the numerical solver includes radiation (Stefan-Boltzmann T⁴)
    // while the analytical eigenfunction expansion does not, so we expect
    // ~15°C divergence over 60 minutes.
    let s = setup(95.0, 80.0, 22.0, 10.0);
    let (max_diff, rms) = compare_liquid_avg(&s, 3600.0, 100);

    assert!(
        max_diff < 16.0,
        "Preheated mug: max diff {:.2}°C exceeds 16°C tolerance",
        max_diff,
    );
    assert!(
        rms < 12.0,
        "Preheated mug: RMS diff {:.2}°C exceeds 12°C tolerance",
        rms,
    );
}

#[test]
fn test_heat_eq_analytical_vs_numerical_cold() {
    // Cold mug: shell starts at room temperature.
    // Numerical solver has radiation; analytical does not.
    let s = setup(95.0, 22.0, 22.0, 10.0);
    let (max_diff, rms) = compare_liquid_avg(&s, 3600.0, 100);

    assert!(
        max_diff < 16.0,
        "Cold mug: max diff {:.2}°C exceeds 16°C tolerance",
        max_diff,
    );
    assert!(
        rms < 12.0,
        "Cold mug: RMS diff {:.2}°C exceeds 12°C tolerance",
        rms,
    );
}

#[test]
fn test_heat_eq_analytical_vs_numerical_isothermal() {
    // Isothermal start: liquid and shell both at 95°C
    // Numerical solver has radiation; analytical does not.
    let s = setup(95.0, 95.0, 22.0, 10.0);
    let (max_diff, rms) = compare_liquid_avg(&s, 3600.0, 100);

    assert!(
        max_diff < 16.0,
        "Isothermal: max diff {:.2}°C exceeds 16°C tolerance",
        max_diff,
    );
    assert!(
        rms < 13.0,
        "Isothermal: RMS diff {:.2}°C exceeds 13°C tolerance",
        rms,
    );
}

/// Grid convergence study: run the cold-mug case at varying resolutions
/// to see how the error scales with grid spacing.
/// Not a pass/fail test — just prints diagnostics.
#[test]
fn convergence_study_cold_mug() {
    use std::time::Instant;

    let water = HeatMaterial::water();
    let ceramic = HeatMaterial::ceramic();
    let air = HeatMaterial::air();

    let volume_ml = 180.0;
    let inner_r_cm = (3.0 * volume_ml / (4.0 * std::f64::consts::PI)).cbrt();
    let wall_mm = 5.0;
    let outer_r_cm = inner_r_cm + wall_mm / 10.0;
    let r_max_mult = 10.0;

    let r_inner = inner_r_cm / 100.0;
    let r_outer = outer_r_cm / 100.0;
    let r_max = outer_r_cm * r_max_mult / 100.0;

    // Analytical reference — closed-form volume average (no radial discretization)
    let analytical_params = ThreeRegionParams {
        r_inner,
        r_outer,
        r_max,
        alpha_1: water.diffusivity(),
        alpha_2: ceramic.diffusivity(),
        alpha_3: air.diffusivity(),
        k_1: water.conductivity,
        k_2: ceramic.conductivity,
        k_3: air.conductivity,
        t_liquid_init: 95.0,
        t_shell_init: 22.0,
        t_ambient: 22.0,
        num_modes: 50,
        h_conv: 10.0,
    };
    let modes = compute_eigenmodes(&analytical_params);
    let analytical_avg = three_region_liquid_average_closed_form(
        &analytical_params, &modes, 3600.0, 50,
    );

    let numerical_params = SphereHeatParams {
        outer_radius_m: r_outer,
        inner_radius_m: r_inner,
        shell_temp_c: 22.0,
        liquid_temp_c: 95.0,
        room_temp_c: 22.0,
        shell_material: ceramic,
        liquid_material: water,
        medium_material: air,
        evaporation_factor: 0.0,
        convection_coeff: 10.0,
        r_max_multiplier: r_max_mult,
    };

    // Test at the current default resolution
    println!("\n=== Grid convergence study (cold mug 95/22°C) ===");
    println!("Analytical reference: 50 modes, closed-form volume average\n");

    let t0 = Instant::now();
    let num_snaps = simulate_sphere_heat_numerical(&numerical_params, 3600.0, 50);
    let elapsed = t0.elapsed();
    let num_avg = liquid_average_temperature(&num_snaps);

    let n_points = num_snaps[0].radii.len();
    let inner_idx = num_snaps[0].inner_interface_idx;
    let dr_at_interface = num_snaps[0].radii[inner_idx] - num_snaps[0].radii[inner_idx - 1];

    println!("Grid: {} points, inner_idx={}, dr at interface = {:.4}mm",
             n_points, inner_idx, dr_at_interface * 1000.0);
    println!("Time: {:.1}ms", elapsed.as_secs_f64() * 1000.0);

    let mut max_diff = 0.0_f64;
    let mut sum_sq = 0.0_f64;
    for i in 0..analytical_avg.len().min(num_avg.len()) {
        let diff = (analytical_avg[i].1 - num_avg[i].1).abs();
        max_diff = max_diff.max(diff);
        sum_sq += diff * diff;
    }
    let rms = (sum_sq / analytical_avg.len() as f64).sqrt();
    println!("max_diff = {:.3}°C, rms = {:.3}°C", max_diff, rms);

    // Print first 10 snapshots for detail
    println!("\nTime-resolved (first 10 snapshots):");
    for i in 0..10.min(analytical_avg.len()) {
        let (t_a, temp_a) = analytical_avg[i];
        let (_, temp_n) = num_avg[i];
        let diff = temp_a - temp_n;
        println!("  t={:5.1}min: analytical={:.2}°C  numerical={:.2}°C  diff={:+.2}°C",
                 t_a / 60.0, temp_a, temp_n, diff);
    }

    // Print grid cells around the inner interface
    println!("\nGrid around r_inner ({:.4}mm):", r_inner * 1000.0);
    let snap = &num_snaps[0];
    let lo = inner_idx.saturating_sub(5);
    let hi = (inner_idx + 5).min(snap.radii.len() - 1);
    for i in lo..=hi {
        let dr = if i > 0 { snap.radii[i] - snap.radii[i - 1] } else { 0.0 };
        let marker = if i == inner_idx { " <-- inner_interface_idx" } else { "" };
        println!("  [{}] r={:.4}mm  T={:.1}°C  dr={:.4}mm{}",
                 i, snap.radii[i] * 1000.0, snap.temperatures[i], dr * 1000.0, marker);
    }
}

// ─── Lumped vs Heat Equation cross-model comparison ───────────────
//
// The lumped and heat equation models differ in two fundamental ways:
//
// 1. **Geometry**: lumped = cylinder (h = 2r); heat equation = sphere.
//    For the same volume, the cylinder has ~14% more surface area → faster cooling.
//
// 2. **Spatial resolution**: lumped assumes uniform temperature within
//    each body; the heat equation resolves the radial T(r) profile.
//
// Despite these structural differences, both models should predict
// qualitatively identical behavior and agree to within ~15°C over
// 60 minutes under convection-only conditions (no radiation/evaporation).

/// Build `MugParameters` for the lumped model from common physical inputs.
fn lumped_params(
    volume_ml: f64,
    wall_mm: f64,
    liquid_temp: f64,
    shell_temp: f64,
    room_temp: f64,
    h_conv: f64,
    preheated: bool,
) -> MugParameters {
    MugParameters {
        material: MugMaterial::Ceramic,
        volume_ml,
        wall_thickness_mm: wall_mm,
        coffee_temp_c: liquid_temp,
        room_temp_c: room_temp,
        preheated,
        preheat_temp_c: shell_temp,
        h_conv,
        emissivity: 0.0,       // disable radiation for cross-model fairness
        exposure_pct: 0.0,     // disable evaporation for cross-model fairness
    }
}

/// Build `CommonSetup` (analytical + numerical heat equation params) from
/// the same physical inputs, also with radiation disabled.
fn heat_eq_setup(
    volume_ml: f64,
    wall_mm: f64,
    liquid_temp: f64,
    shell_temp: f64,
    room_temp: f64,
    h_conv: f64,
) -> CommonSetup {
    let water = HeatMaterial::water();
    let ceramic = HeatMaterial::ceramic();
    let air = HeatMaterial::air();

    let inner_r_cm = (3.0 * volume_ml / (4.0 * std::f64::consts::PI)).cbrt();
    let outer_r_cm = inner_r_cm + wall_mm / 10.0;
    let r_max_mult = 10.0;

    let r_inner = inner_r_cm / 100.0;
    let r_outer = outer_r_cm / 100.0;
    let r_max = outer_r_cm * r_max_mult / 100.0;

    CommonSetup {
        analytical: ThreeRegionParams {
            r_inner,
            r_outer,
            r_max,
            alpha_1: water.diffusivity(),
            alpha_2: ceramic.diffusivity(),
            alpha_3: air.diffusivity(),
            k_1: water.conductivity,
            k_2: ceramic.conductivity,
            k_3: air.conductivity,
            t_liquid_init: liquid_temp,
            t_shell_init: shell_temp,
            t_ambient: room_temp,
            num_modes: 50,
            h_conv,  // convection only, no radiation enhancement
        },
        numerical: SphereHeatParams {
            outer_radius_m: r_outer,
            inner_radius_m: r_inner,
            shell_temp_c: shell_temp,
            liquid_temp_c: liquid_temp,
            room_temp_c: room_temp,
            shell_material: ceramic,
            liquid_material: water,
            medium_material: air,
            evaporation_factor: 0.0,
            convection_coeff: h_conv,
            r_max_multiplier: r_max_mult,
        },
    }
}

/// Linearly interpolate a time series at a given time.
fn interpolate(series: &[(f64, f64)], t: f64) -> f64 {
    if series.is_empty() {
        return 0.0;
    }
    if t <= series[0].0 {
        return series[0].1;
    }
    for i in 1..series.len() {
        if series[i].0 >= t {
            let frac = (t - series[i - 1].0) / (series[i].0 - series[i - 1].0);
            return series[i - 1].1 + frac * (series[i].1 - series[i - 1].1);
        }
    }
    series.last().unwrap().1
}

/// Compare lumped vs heat equation coffee temperature curves.
/// Returns (max_diff, rms_diff) in °C.
fn compare_lumped_vs_heat_eq(
    lumped_series: &[(f64, f64)],         // (time_minutes, temp_c)
    heat_eq_series: &[(f64, f64)],        // (time_seconds, temp_c)
    duration_min: f64,
    label: &str,
) -> (f64, f64) {
    // Sample at 1-minute intervals
    let n_samples = duration_min as usize;
    let mut max_diff = 0.0_f64;
    let mut sum_sq = 0.0_f64;
    let mut count = 0;

    println!("\n  {label}:");
    for min in 0..=n_samples {
        let t_min = min as f64;
        let t_sec = t_min * 60.0;
        let t_lumped = interpolate(lumped_series, t_min);
        let t_heat_eq = interpolate(heat_eq_series, t_sec);
        let diff = (t_lumped - t_heat_eq).abs();
        max_diff = max_diff.max(diff);
        sum_sq += diff * diff;
        count += 1;

        if min % 10 == 0 || diff > 10.0 {
            println!(
                "    t={:3}min: lumped={:.1}°C, heat_eq={:.1}°C, diff={:.1}°C",
                min, t_lumped, t_heat_eq, diff,
            );
        }
    }

    let rms = (sum_sq / count as f64).sqrt();
    println!("    max_diff={:.1}°C, rms={:.1}°C", max_diff, rms);
    (max_diff, rms)
}

#[test]
fn test_lumped_analytical_vs_heat_eq_analytical_cold() {
    let volume = 180.0;
    let wall = 5.0;
    let liquid = 95.0;
    let room = 22.0;
    let h = 10.0;
    let duration_min = 60;

    // Lumped analytical — cold mug (shell starts at room temp)
    let lp = lumped_params(volume, wall, liquid, room, room, h, false);
    let lumped_data = simulate_lumped_analytical(&lp, duration_min);
    let lumped_coffee: Vec<(f64, f64)> = CoolingDataPoint::coffee_series(&lumped_data);

    // Heat equation analytical — cold mug
    let he = heat_eq_setup(volume, wall, liquid, room, room, h);
    let modes = compute_eigenmodes(&he.analytical);
    let heat_eq_avg = three_region_liquid_average_closed_form(
        &he.analytical, &modes, duration_min as f64 * 60.0, duration_min * 5,
    );

    let (max_diff, rms) = compare_lumped_vs_heat_eq(
        &lumped_coffee, &heat_eq_avg, duration_min as f64,
        "Analytical lumped vs analytical heat equation (cold mug)",
    );

    // Both analytical solvers omit radiation, so the only difference is
    // geometry (cylinder vs sphere) and lumped-vs-resolved spatial detail.
    // Observed: ~5°C max, ~4°C RMS.
    assert!(
        max_diff < 8.0,
        "Analytical cold: max diff {:.1}°C exceeds 8°C",
        max_diff,
    );
    assert!(
        rms < 6.0,
        "Analytical cold: RMS {:.1}°C exceeds 6°C",
        rms,
    );

    // Both should start at the same temperature
    assert!(
        (lumped_coffee[0].1 - heat_eq_avg[0].1).abs() < 0.5,
        "Both should start at ~{}°C", liquid,
    );
}

#[test]
fn test_lumped_analytical_vs_heat_eq_analytical_preheated() {
    let volume = 180.0;
    let wall = 5.0;
    let liquid = 95.0;
    let shell = 80.0;
    let room = 22.0;
    let h = 10.0;
    let duration_min = 60;

    let lp = lumped_params(volume, wall, liquid, shell, room, h, true);
    let lumped_data = simulate_lumped_analytical(&lp, duration_min);
    let lumped_coffee: Vec<(f64, f64)> = CoolingDataPoint::coffee_series(&lumped_data);

    let he = heat_eq_setup(volume, wall, liquid, shell, room, h);
    let modes = compute_eigenmodes(&he.analytical);
    let heat_eq_avg = three_region_liquid_average_closed_form(
        &he.analytical, &modes, duration_min as f64 * 60.0, duration_min * 5,
    );

    let (max_diff, rms) = compare_lumped_vs_heat_eq(
        &lumped_coffee, &heat_eq_avg, duration_min as f64,
        "Analytical lumped vs analytical heat equation (preheated mug)",
    );

    // Both analytical, no radiation — geometry-only difference.
    // Observed: ~5°C max, ~4°C RMS.
    assert!(
        max_diff < 8.0,
        "Analytical preheated: max diff {:.1}°C exceeds 8°C",
        max_diff,
    );
    assert!(
        rms < 6.0,
        "Analytical preheated: RMS {:.1}°C exceeds 6°C",
        rms,
    );
}

#[test]
fn test_lumped_numerical_vs_heat_eq_numerical_cold() {
    let volume = 180.0;
    let wall = 5.0;
    let liquid = 95.0;
    let room = 22.0;
    let h = 10.0;
    let duration_min = 60;

    // Lumped numerical — cold mug
    let lp = lumped_params(volume, wall, liquid, room, room, h, false);
    let lumped_data = simulate_lumped_numerical(&lp, duration_min);
    let lumped_coffee: Vec<(f64, f64)> = CoolingDataPoint::coffee_series(&lumped_data);

    // Heat equation numerical — cold mug
    let he = heat_eq_setup(volume, wall, liquid, room, room, h);
    let numerical_snaps = simulate_sphere_heat_numerical(
        &he.numerical, duration_min as f64 * 60.0, duration_min * 5,
    );
    let heat_eq_avg = liquid_average_temperature(&numerical_snaps);

    let (max_diff, rms) = compare_lumped_vs_heat_eq(
        &lumped_coffee, &heat_eq_avg, duration_min as f64,
        "Numerical lumped vs numerical heat equation (cold mug)",
    );

    // The numerical heat equation includes Stefan-Boltzmann radiation from the
    // ceramic shell (emissivity ~0.9) while the lumped model has emissivity=0,
    // so the heat equation cools faster. Combined with the geometry difference,
    // this gives a larger gap than the analytical comparison.
    // Observed: ~13°C max, ~11°C RMS.
    assert!(
        max_diff < 18.0,
        "Numerical cold: max diff {:.1}°C exceeds 18°C",
        max_diff,
    );
    assert!(
        rms < 14.0,
        "Numerical cold: RMS {:.1}°C exceeds 14°C",
        rms,
    );

    assert!(
        (lumped_coffee[0].1 - heat_eq_avg[0].1).abs() < 0.5,
        "Both should start at ~{}°C", liquid,
    );
}

#[test]
fn test_lumped_numerical_vs_heat_eq_numerical_preheated() {
    let volume = 180.0;
    let wall = 5.0;
    let liquid = 95.0;
    let shell = 80.0;
    let room = 22.0;
    let h = 10.0;
    let duration_min = 60;

    let lp = lumped_params(volume, wall, liquid, shell, room, h, true);
    let lumped_data = simulate_lumped_numerical(&lp, duration_min);
    let lumped_coffee: Vec<(f64, f64)> = CoolingDataPoint::coffee_series(&lumped_data);

    let he = heat_eq_setup(volume, wall, liquid, shell, room, h);
    let numerical_snaps = simulate_sphere_heat_numerical(
        &he.numerical, duration_min as f64 * 60.0, duration_min * 5,
    );
    let heat_eq_avg = liquid_average_temperature(&numerical_snaps);

    let (max_diff, rms) = compare_lumped_vs_heat_eq(
        &lumped_coffee, &heat_eq_avg, duration_min as f64,
        "Numerical lumped vs numerical heat equation (preheated mug)",
    );

    // Radiation asymmetry + geometry → larger gap for preheated case.
    // Observed: ~15.5°C max, ~12°C RMS.
    assert!(
        max_diff < 18.0,
        "Numerical preheated: max diff {:.1}°C exceeds 18°C",
        max_diff,
    );
    assert!(
        rms < 14.0,
        "Numerical preheated: RMS {:.1}°C exceeds 14°C",
        rms,
    );
}

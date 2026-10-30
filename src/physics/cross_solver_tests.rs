//! Cross-solver comparison tests.
//!
//! Verifies that analytical and numerical solvers produce consistent
//! temperature curves when given identical physical parameters.
//! Covers both the lumped capacitance model and the three-region
//! eigenfunction expansion / finite-difference heat equation.

use super::*;
use heat_equation_analytical::{
    relinearised_liquid_average, three_region_liquid_average_closed_form,
};

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
        exposure: 0.2,
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
        assert!(
            dt < 0.01,
            "Time mismatch at index {}: {:.3} vs {:.3}",
            i,
            numerical[i].time_minutes,
            analytical[i].time_minutes
        );

        let coffee_diff = (numerical[i].coffee_temp_c - analytical[i].coffee_temp_c).abs();
        max_coffee_diff = max_coffee_diff.max(coffee_diff);
        assert!(
            coffee_diff < 0.5,
            "Coffee temp mismatch at t={:.1}: numerical={:.2}, analytical={:.2}, diff={:.3}",
            numerical[i].time_minutes,
            numerical[i].coffee_temp_c,
            analytical[i].coffee_temp_c,
            coffee_diff
        );

        let mug_diff = (numerical[i].mug_temp_c - analytical[i].mug_temp_c).abs();
        max_mug_diff = max_mug_diff.max(mug_diff);
        assert!(
            mug_diff < 0.5,
            "Mug temp mismatch at t={:.1}: numerical={:.2}, analytical={:.2}, diff={:.3}",
            numerical[i].time_minutes,
            numerical[i].mug_temp_c,
            analytical[i].mug_temp_c,
            mug_diff
        );
    }
    println!(
        "Cold mug: max coffee diff={:.4}°C, max mug diff={:.4}°C",
        max_coffee_diff, max_mug_diff
    );
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
        exposure: 0.2,
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
            numerical[i].time_minutes,
            coffee_diff
        );
    }
    println!("Preheated mug: max coffee diff={:.4}°C", max_diff);
}

// ─── Three-region heat equation ────────────────────────────────────

/// Physical inputs for one heat-equation case, shared by both solvers.
#[derive(Debug, Clone, Copy)]
struct HeatEqCase {
    volume_ml: f64,
    wall_mm: f64,
    liquid_temp: f64,
    shell_temp: f64,
    room_temp: f64,
    h_conv: f64,
    /// Ceramic shell emissivity; 0 disables outer radiation in both solvers.
    shell_emissivity: f64,
    /// Lid openness 0–1.
    exposure: f64,
}

impl HeatEqCase {
    /// The panel's default mug: 180 ml ceramic, 5 mm wall, ceramic emissivity.
    fn ceramic_180ml(liquid_temp: f64, shell_temp: f64, room_temp: f64, h_conv: f64) -> Self {
        Self {
            volume_ml: 180.0,
            wall_mm: 5.0,
            liquid_temp,
            shell_temp,
            room_temp,
            h_conv,
            shell_emissivity: MugMaterial::Ceramic.emissivity(),
            exposure: 0.0,
        }
    }

    /// Both solvers' params, built the way `heat_equation_panel.rs` builds
    /// them, so the two see identical physics: the same materials, outer
    /// radiation (numerical: T⁴ via the shell emissivity; analytical: the
    /// segment's linearised `h_rad` added to `h_conv`), and the same lid
    /// openness (numerical: `exposure`; analytical: `three_region_h_inner`).
    fn build(&self) -> CommonSetup {
        let water = HeatMaterial::water();
        let shell = HeatMaterial::from_mug(MugMaterial::Ceramic, self.shell_emissivity);
        let air = HeatMaterial::air();

        let inner_r_cm = (3.0 * self.volume_ml / (4.0 * std::f64::consts::PI)).cbrt();
        let outer_r_cm = inner_r_cm + self.wall_mm / 10.0;
        let r_max_mult = 10.0;

        let r_inner = inner_r_cm / 100.0;
        let r_outer = outer_r_cm / 100.0;
        let r_max = outer_r_cm * r_max_mult / 100.0;

        CommonSetup {
            analytical: ThreeRegionMug {
                r_inner,
                r_outer,
                r_max,
                shell,
                liquid_temp_c: self.liquid_temp,
                shell_temp_c: self.shell_temp,
                room_temp_c: self.room_temp,
                h_conv: self.h_conv,
                exposure: self.exposure,
                num_modes: 50,
            },
            numerical: SphereHeatParams {
                outer_radius_m: r_outer,
                inner_radius_m: r_inner,
                shell_temp_c: self.shell_temp,
                liquid_temp_c: self.liquid_temp,
                room_temp_c: self.room_temp,
                shell_material: shell,
                liquid_material: water,
                medium_material: air,
                exposure: self.exposure,
                convection_coeff: self.h_conv,
                r_max_multiplier: r_max_mult,
                outer_boundary: OuterBoundary::Ambient,
            },
        }
    }
}

/// Matching params for both heat-equation solvers.
struct CommonSetup {
    analytical: ThreeRegionMug,
    numerical: SphereHeatParams,
}

/// Compare liquid average temperature curves from both solvers.
/// Returns (max_diff, rms_diff) in °C.
fn compare_liquid_avg(setup: &CommonSetup, duration_s: f64, num_snapshots: usize) -> (f64, f64) {
    // Analytical — closed-form volume average eliminates all radial discretization
    let analytical_avg = relinearised_liquid_average(
        &setup.analytical,
        ANALYTICAL_SEGMENTS,
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

/// Within-family agreement bounds at lid openness `exposure`, as
/// (max, RMS) difference of the liquid average in °C. Both heat-equation
/// solvers integrate the same PDE with the same physics; the modelled
/// difference is that the analytical one linearises radiation and the
/// lid-opening loss per segment, at the segment's liquid-average and
/// outer-surface temperatures, while the numerical one applies the open-top
/// loss at the open surface's own temperature every step. The open surface
/// is cooler than the liquid average, so the analytical coffee runs a little
/// cold, more so as the lid opens. Measured for the 180 ml mug over an hour
/// (max / RMS): sealed preheated 0.05 / 0.03, sealed cold 0.22 / 0.14, 20%
/// open 0.20 / 0.14, no lid 1.74 / 0.77 (preheated) and 1.82 / 0.95 °C
/// (cold). The bounds follow that rise with openness plus a ~10% margin.
fn heat_eq_bounds(exposure: f64) -> (f64, f64) {
    (0.4 + 1.6 * exposure, 0.25 + 0.8 * exposure)
}

fn assert_heat_eq_solvers_agree(label: &str, case: HeatEqCase) {
    let (max_bound, rms_bound) = heat_eq_bounds(case.exposure);
    let (max_diff, rms) = compare_liquid_avg(&case.build(), 3600.0, 100);
    assert!(
        max_diff < max_bound,
        "{label}: max diff {max_diff:.2}°C exceeds {max_bound}°C"
    );
    assert!(
        rms < rms_bound,
        "{label}: RMS diff {rms:.2}°C exceeds {rms_bound}°C"
    );
}

#[test]
fn test_heat_eq_analytical_vs_numerical_preheated() {
    assert_heat_eq_solvers_agree(
        "Preheated",
        HeatEqCase::ceramic_180ml(95.0, 80.0, 22.0, 10.0),
    );
}

#[test]
fn test_heat_eq_analytical_vs_numerical_cold() {
    assert_heat_eq_solvers_agree("Cold", HeatEqCase::ceramic_180ml(95.0, 22.0, 22.0, 10.0));
}

#[test]
fn test_heat_eq_analytical_vs_numerical_isothermal() {
    assert_heat_eq_solvers_agree(
        "Isothermal",
        HeatEqCase::ceramic_180ml(95.0, 95.0, 22.0, 10.0),
    );
}

#[test]
fn test_heat_eq_analytical_vs_numerical_exposed() {
    // The site's default lid openness: the open-top loss acts in both
    // solvers (numerically exact, analytically as h_inner).
    let case = HeatEqCase {
        exposure: 0.2,
        ..HeatEqCase::ceramic_180ml(95.0, 80.0, 22.0, 10.0)
    };
    assert_heat_eq_solvers_agree("Exposed 20%", case);
}

#[test]
fn test_heat_eq_analytical_vs_numerical_no_lid_preheated() {
    // No lid at the site's default convection: the open-top loss is at its
    // largest, so this is where linearising at the liquid average differs most.
    let case = HeatEqCase {
        exposure: 1.0,
        ..HeatEqCase::ceramic_180ml(95.0, 80.0, 22.0, 20.0)
    };
    assert_heat_eq_solvers_agree("No lid, preheated", case);
}

#[test]
fn test_heat_eq_analytical_vs_numerical_no_lid_cold() {
    let case = HeatEqCase {
        exposure: 1.0,
        ..HeatEqCase::ceramic_180ml(95.0, 22.0, 22.0, 20.0)
    };
    assert_heat_eq_solvers_agree("No lid, cold", case);
}

/// Systematic factor isolation: test each physics feature independently to
/// identify which causes the largest analytical/numerical divergence.
///
/// Factors tested:
///   1. Pure conduction (h_conv=0, emissivity=0, evaporation=0)
///   2. Convection only (h_conv=10, emissivity=0, evaporation=0)
///   3. Convection + linearized radiation (h_conv=10, emissivity=0.9)
///   4. Panel-like radiation (analytical gets h_conv_eff = h_conv + h_rad)
#[test]
fn test_factor_isolation() {
    let water = HeatMaterial::water();
    let ceramic_no_rad = HeatMaterial {
        emissivity: 0.0,
        ..HeatMaterial::ceramic()
    };
    let ceramic_with_rad = HeatMaterial::ceramic(); // emissivity = 0.90
    let air = HeatMaterial::air();

    let volume_ml = 180.0;
    let inner_r_cm = (3.0 * volume_ml / (4.0 * std::f64::consts::PI)).cbrt();
    let wall_mm = 5.0;
    let outer_r_cm = inner_r_cm + wall_mm / 10.0;
    let r_max_mult = 10.0;

    let r_inner = inner_r_cm / 100.0;
    let r_outer = outer_r_cm / 100.0;
    let r_max = outer_r_cm * r_max_mult / 100.0;

    let liquid_temp = 95.0;
    let shell_temp = 22.0; // cold mug — worst case for thermal gradients
    let room_temp = 22.0;

    // Helper to build analytical params
    let analytical_params = |h_conv: f64| -> ThreeRegionParams {
        ThreeRegionParams {
            r_inner,
            r_outer,
            r_max,
            alpha_1: water.diffusivity(),
            alpha_2: ceramic_no_rad.diffusivity(),
            alpha_3: air.diffusivity(),
            k_1: water.conductivity,
            k_2: ceramic_no_rad.conductivity,
            k_3: air.conductivity,
            t_liquid_init: liquid_temp,
            t_shell_init: shell_temp,
            t_ambient: room_temp,
            num_modes: 50,
            h_conv,
            h_inner: 0.0,
        }
    };

    // Helper to build numerical params
    let numerical_params =
        |h_conv: f64, shell_mat: HeatMaterial, exposure: f64| -> SphereHeatParams {
            SphereHeatParams {
                outer_radius_m: r_outer,
                inner_radius_m: r_inner,
                shell_temp_c: shell_temp,
                liquid_temp_c: liquid_temp,
                room_temp_c: room_temp,
                shell_material: shell_mat,
                liquid_material: water,
                medium_material: air,
                exposure,
                convection_coeff: h_conv,
                r_max_multiplier: r_max_mult,
                outer_boundary: OuterBoundary::Ambient,
            }
        };

    let duration = 3600.0;
    let snaps = 100;

    println!("\n======================================================================");
    println!("FACTOR ISOLATION: cold mug (liquid=95, shell=22, room=22)");
    println!("======================================================================");

    // ── 1. Pure conduction ──
    {
        let ap = analytical_params(0.0);
        let np = numerical_params(0.0, ceramic_no_rad, 0.0);
        let modes = compute_eigenmodes(&ap);
        let a_avg = three_region_liquid_average_closed_form(&ap, &modes, duration, snaps);
        let n_snaps = simulate_sphere_heat_numerical(&np, duration, snaps);
        let n_avg = liquid_average_temperature(&n_snaps);

        let mut max_diff = 0.0_f64;
        let mut sum_sq = 0.0;
        for i in 0..a_avg.len().min(n_avg.len()) {
            let diff = (a_avg[i].1 - n_avg[i].1).abs();
            max_diff = max_diff.max(diff);
            sum_sq += diff * diff;
        }
        let rms = (sum_sq / a_avg.len() as f64).sqrt();
        println!("\n1. PURE CONDUCTION (h=0, ε=0, evap=0):");
        println!("   max_diff={:.2}°C, rms={:.2}°C", max_diff, rms);
        println!(
            "   t=0:  analytical={:.1}°C, numerical={:.1}°C",
            a_avg[0].1, n_avg[0].1
        );
        println!(
            "   t=30: analytical={:.1}°C, numerical={:.1}°C",
            a_avg[snaps / 2].1,
            n_avg[snaps / 2].1
        );
        println!(
            "   t=60: analytical={:.1}°C, numerical={:.1}°C",
            a_avg[snaps - 1].1,
            n_avg[snaps - 1].1
        );
    }

    // ── 2. Convection only ──
    {
        let h_conv = 10.0;
        let ap = analytical_params(h_conv);
        let np = numerical_params(h_conv, ceramic_no_rad, 0.0);
        let modes = compute_eigenmodes(&ap);
        let a_avg = three_region_liquid_average_closed_form(&ap, &modes, duration, snaps);
        let n_snaps = simulate_sphere_heat_numerical(&np, duration, snaps);
        let n_avg = liquid_average_temperature(&n_snaps);

        let mut max_diff = 0.0_f64;
        let mut sum_sq = 0.0;
        for i in 0..a_avg.len().min(n_avg.len()) {
            let diff = (a_avg[i].1 - n_avg[i].1).abs();
            max_diff = max_diff.max(diff);
            sum_sq += diff * diff;
        }
        let rms = (sum_sq / a_avg.len() as f64).sqrt();
        println!("\n2. CONVECTION ONLY (h=10, ε=0, evap=0):");
        println!("   max_diff={:.2}°C, rms={:.2}°C", max_diff, rms);
        println!(
            "   t=0:  analytical={:.1}°C, numerical={:.1}°C",
            a_avg[0].1, n_avg[0].1
        );
        println!(
            "   t=30: analytical={:.1}°C, numerical={:.1}°C",
            a_avg[snaps / 2].1,
            n_avg[snaps / 2].1
        );
        println!(
            "   t=60: analytical={:.1}°C, numerical={:.1}°C",
            a_avg[snaps - 1].1,
            n_avg[snaps - 1].1
        );
    }

    // ── 3. Convection + radiation (numerical has radiation, analytical does NOT) ──
    {
        let h_conv = 10.0;
        let ap = analytical_params(h_conv); // no radiation in eigenvalue problem
        let np = numerical_params(h_conv, ceramic_with_rad, 0.0);
        let modes = compute_eigenmodes(&ap);
        let a_avg = three_region_liquid_average_closed_form(&ap, &modes, duration, snaps);
        let n_snaps = simulate_sphere_heat_numerical(&np, duration, snaps);
        let n_avg = liquid_average_temperature(&n_snaps);

        let mut max_diff = 0.0_f64;
        let mut sum_sq = 0.0;
        for i in 0..a_avg.len().min(n_avg.len()) {
            let diff = (a_avg[i].1 - n_avg[i].1).abs();
            max_diff = max_diff.max(diff);
            sum_sq += diff * diff;
        }
        let rms = (sum_sq / a_avg.len() as f64).sqrt();
        println!("\n3. CONV + RAD (numerical only) (h=10, ε=0.9 numerical / 0 analytical):");
        println!("   max_diff={:.2}°C, rms={:.2}°C", max_diff, rms);
        println!(
            "   t=0:  analytical={:.1}°C, numerical={:.1}°C",
            a_avg[0].1, n_avg[0].1
        );
        println!(
            "   t=30: analytical={:.1}°C, numerical={:.1}°C",
            a_avg[snaps / 2].1,
            n_avg[snaps / 2].1
        );
        println!(
            "   t=60: analytical={:.1}°C, numerical={:.1}°C",
            a_avg[snaps - 1].1,
            n_avg[snaps - 1].1
        );
    }

    // ── 4. Panel-like: conv + linearized radiation in BOTH ──
    //       analytical gets h_eff = h_conv + h_rad_outer (fixed linearization)
    //       numerical gets time-varying h_rad (re-linearized each step)
    {
        let h_conv = 10.0;
        let h_conv_eff = heat_equation_analytical::outer_surface_h(
            h_conv,
            ceramic_with_rad.emissivity,
            liquid_temp,
            room_temp,
        );
        let h_rad = h_conv_eff - h_conv;
        println!(
            "\n4. PANEL-LIKE: linearized radiation in both (h_rad={:.1} → h_eff={:.1}):",
            h_rad, h_conv_eff
        );

        let ap = analytical_params(h_conv_eff);
        let np = numerical_params(h_conv, ceramic_with_rad, 0.0);
        let modes = compute_eigenmodes(&ap);
        let a_avg = three_region_liquid_average_closed_form(&ap, &modes, duration, snaps);
        let n_snaps = simulate_sphere_heat_numerical(&np, duration, snaps);
        let n_avg = liquid_average_temperature(&n_snaps);

        let mut max_diff = 0.0_f64;
        let mut sum_sq = 0.0;
        for i in 0..a_avg.len().min(n_avg.len()) {
            let diff = (a_avg[i].1 - n_avg[i].1).abs();
            max_diff = max_diff.max(diff);
            sum_sq += diff * diff;
        }
        let rms = (sum_sq / a_avg.len() as f64).sqrt();
        println!("   max_diff={:.2}°C, rms={:.2}°C", max_diff, rms);
        println!(
            "   t=0:  analytical={:.1}°C, numerical={:.1}°C",
            a_avg[0].1, n_avg[0].1
        );
        println!(
            "   t=30: analytical={:.1}°C, numerical={:.1}°C",
            a_avg[snaps / 2].1,
            n_avg[snaps / 2].1
        );
        println!(
            "   t=60: analytical={:.1}°C, numerical={:.1}°C",
            a_avg[snaps - 1].1,
            n_avg[snaps - 1].1
        );
    }

    // ── 5. High convection, no radiation ──
    {
        let h_conv = 25.0;
        let ap = analytical_params(h_conv);
        let np = numerical_params(h_conv, ceramic_no_rad, 0.0);
        let modes = compute_eigenmodes(&ap);
        let a_avg = three_region_liquid_average_closed_form(&ap, &modes, duration, snaps);
        let n_snaps = simulate_sphere_heat_numerical(&np, duration, snaps);
        let n_avg = liquid_average_temperature(&n_snaps);

        let mut max_diff = 0.0_f64;
        let mut sum_sq = 0.0;
        for i in 0..a_avg.len().min(n_avg.len()) {
            let diff = (a_avg[i].1 - n_avg[i].1).abs();
            max_diff = max_diff.max(diff);
            sum_sq += diff * diff;
        }
        let rms = (sum_sq / a_avg.len() as f64).sqrt();
        println!("\n5. HIGH CONVECTION (h=25, ε=0, evap=0):");
        println!("   max_diff={:.2}°C, rms={:.2}°C", max_diff, rms);
        println!(
            "   t=0:  analytical={:.1}°C, numerical={:.1}°C",
            a_avg[0].1, n_avg[0].1
        );
        println!(
            "   t=30: analytical={:.1}°C, numerical={:.1}°C",
            a_avg[snaps / 2].1,
            n_avg[snaps / 2].1
        );
        println!(
            "   t=60: analytical={:.1}°C, numerical={:.1}°C",
            a_avg[snaps - 1].1,
            n_avg[snaps - 1].1
        );
    }
}

/// Exposure isolation: both heat-equation solvers with the lid opening on.
///
///   Analytical: 3-region eigenfunction with `h_inner` from
///               `surface::three_region_h_inner`, re-linearised per segment
///   Numerical:  open-top loss at the inner interface, re-linearised each step
///
/// Sub-tests:
///   A. Vary exposure (0%, 10%, 20%, 50%, 100%) — panel flow, cold mug
///   B. Time-resolved detail at 20% exposure
///   C. Exposure with no outer radiation, to isolate the inner-surface loss
///   D. Preheated mug with exposure
#[test]
fn test_exposure_isolation() {
    use surface::three_region_h_inner;

    let water = HeatMaterial::water();
    let ceramic = HeatMaterial::ceramic();
    let ceramic_no_rad = HeatMaterial {
        emissivity: 0.0,
        ..ceramic
    };
    let air = HeatMaterial::air();

    let volume_ml = 180.0;
    let inner_r_cm = (3.0 * volume_ml / (4.0 * std::f64::consts::PI)).cbrt();
    let wall_mm = 5.0;
    let outer_r_cm = inner_r_cm + wall_mm / 10.0;
    let r_max_mult = 10.0;
    let h_conv = 10.0;

    let r_inner = inner_r_cm / 100.0;
    let r_outer = outer_r_cm / 100.0;
    let r_max = outer_r_cm * r_max_mult / 100.0;

    let duration = 3600.0;
    let snaps = 100;
    let liquid_temp = 95.0;
    let room_temp = 22.0;

    println!("\n======================================================================");
    println!("EXPOSURE ISOLATION: liquid=95, room=22");
    println!("======================================================================");

    // The first segment's single-linearisation h_conv_eff and h_inner, for reference
    let h_conv_eff = heat_equation_analytical::outer_surface_h(
        h_conv,
        ceramic.emissivity,
        liquid_temp,
        room_temp,
    );
    let h_rad_outer = h_conv_eff - h_conv;
    println!(
        "Analytical h_conv_eff = {:.1} (h_conv={:.1} + h_rad_outer={:.1}); h_inner at 100% = {:.1}",
        h_conv_eff,
        h_conv,
        h_rad_outer,
        three_region_h_inner(1.0, 0.5 * (liquid_temp + room_temp), room_temp, h_conv)
    );

    let analytical_avg = |shell: HeatMaterial, shell_temp: f64, exposure: f64| -> Vec<(f64, f64)> {
        let mug = ThreeRegionMug {
            r_inner,
            r_outer,
            r_max,
            shell,
            liquid_temp_c: liquid_temp,
            shell_temp_c: shell_temp,
            room_temp_c: room_temp,
            h_conv,
            exposure,
            num_modes: 50,
        };
        relinearised_liquid_average(&mug, ANALYTICAL_SEGMENTS, duration, snaps)
    };
    let numerical_avg = |shell: HeatMaterial, shell_temp: f64, exposure: f64| -> Vec<(f64, f64)> {
        let np = SphereHeatParams {
            outer_radius_m: r_outer,
            inner_radius_m: r_inner,
            shell_temp_c: shell_temp,
            liquid_temp_c: liquid_temp,
            room_temp_c: room_temp,
            shell_material: shell,
            liquid_material: water,
            medium_material: air,
            exposure,
            convection_coeff: h_conv,
            r_max_multiplier: r_max_mult,
            outer_boundary: OuterBoundary::Ambient,
        };
        liquid_average_temperature(&simulate_sphere_heat_numerical(&np, duration, snaps))
    };
    let print_row = |exposure: f64, a_avg: &[(f64, f64)], n_avg: &[(f64, f64)]| {
        let mut max_diff = 0.0_f64;
        let mut sum_sq = 0.0;
        for i in 0..a_avg.len().min(n_avg.len()) {
            let diff = (a_avg[i].1 - n_avg[i].1).abs();
            max_diff = max_diff.max(diff);
            sum_sq += diff * diff;
        }
        let rms = (sum_sq / a_avg.len() as f64).sqrt();
        let last = a_avg.len() - 1;
        println!(
            "{:<12} {:>10.2} {:>10.2} {:>10.1} {:>10.1}",
            format!("{:.0}%", exposure * 100.0),
            max_diff,
            rms,
            a_avg[last].1,
            n_avg[last].1
        );
    };
    let header = || {
        println!(
            "{:<12} {:>10} {:>10} {:>10} {:>10}",
            "exposure", "max_diff", "rms", "a_t=60", "n_t=60"
        )
    };

    // ── A. Vary exposure level, cold mug ──
    println!("\n--- A. Varying exposure level (cold mug) ---");
    header();
    for &exposure in &[0.0, 0.10, 0.20, 0.50, 1.0] {
        let a_avg = analytical_avg(ceramic, room_temp, exposure);
        let n_avg = numerical_avg(ceramic, room_temp, exposure);
        print_row(exposure, &a_avg, &n_avg);
    }

    // ── B. Time-resolved detail at 20% exposure ──
    println!("\n--- B. Time-resolved at 20% exposure (sealed vs open analytical vs numerical) ---");
    {
        let exposure = 0.20;
        let a_sealed = analytical_avg(ceramic, room_temp, 0.0);
        let a_open = analytical_avg(ceramic, room_temp, exposure);
        let n_avg = numerical_avg(ceramic, room_temp, exposure);
        println!(
            "{:<8} {:>10} {:>10} {:>10} {:>10}",
            "t(min)", "sealed", "open", "numerical", "open-num"
        );
        for i in (0..snaps).step_by(10) {
            println!(
                "{:<8.1} {:>10.1} {:>10.1} {:>10.1} {:>+10.1}",
                a_open[i].0 / 60.0,
                a_sealed[i].1,
                a_open[i].1,
                n_avg[i].1,
                a_open[i].1 - n_avg[i].1
            );
        }
    }

    // ── C. Exposure with no outer radiation (isolate inner-surface effects) ──
    println!("\n--- C. Exposure isolation (no outer radiation, ε=0 on shell) ---");
    header();
    for &exposure in &[0.0, 0.10, 0.20, 0.50, 1.0] {
        let a_avg = analytical_avg(ceramic_no_rad, room_temp, exposure);
        let n_avg = numerical_avg(ceramic_no_rad, room_temp, exposure);
        print_row(exposure, &a_avg, &n_avg);
    }

    // ── D. Preheated mug with exposure ──
    println!("\n--- D. Preheated mug (shell=80°C) with exposure ---");
    header();
    for &exposure in &[0.0, 0.20, 0.50, 1.0] {
        let a_avg = analytical_avg(ceramic, 80.0, exposure);
        let n_avg = numerical_avg(ceramic, 80.0, exposure);
        print_row(exposure, &a_avg, &n_avg);
    }
}

/// End-to-end panel-like comparison: analytical (grid-based avg, the
/// re-linearised solver) vs numerical, exactly as the UI computes them.
#[test]
fn test_panel_endtoend_comparison() {
    use super::convert_snapshots;

    let water = HeatMaterial::water();
    let ceramic = HeatMaterial::ceramic(); // emissivity 0.9
    let air = HeatMaterial::air();

    let volume_ml = 180.0;
    let inner_r_cm = (3.0 * volume_ml / (4.0 * std::f64::consts::PI)).cbrt();
    let wall_mm = 5.0;
    let outer_r_cm = inner_r_cm + wall_mm / 10.0;
    let r_max_mult = 10.0;
    let h_conv = 10.0;
    let display_r_max = 0.10;

    let r_inner = inner_r_cm / 100.0;
    let r_outer = outer_r_cm / 100.0;
    let r_max = outer_r_cm * r_max_mult / 100.0;

    let duration = 3600.0;
    let snaps = 100;
    let liquid_temp = 95.0;
    let room_temp = 22.0;

    for &(label, shell_temp) in &[("COLD MUG", 22.0), ("PREHEATED", 80.0)] {
        let mug = ThreeRegionMug {
            r_inner,
            r_outer,
            r_max,
            shell: ceramic,
            liquid_temp_c: liquid_temp,
            shell_temp_c: shell_temp,
            room_temp_c: room_temp,
            h_conv,
            exposure: 0.0,
            num_modes: 50,
        };

        // Grid-based analytical (same as panel with 80 pts)
        let a_snaps = simulate_three_region_relinearised(
            &mug,
            ANALYTICAL_SEGMENTS,
            duration,
            snaps,
            80,
            display_r_max,
        );
        let a_heat = convert_snapshots(a_snaps);
        let a_avg = liquid_average_temperature(&a_heat);

        // Closed-form analytical (ground truth for analytical solver)
        let a_closed = relinearised_liquid_average(&mug, ANALYTICAL_SEGMENTS, duration, snaps);

        // Numerical
        let np = SphereHeatParams {
            outer_radius_m: r_outer,
            inner_radius_m: r_inner,
            shell_temp_c: shell_temp,
            liquid_temp_c: liquid_temp,
            room_temp_c: room_temp,
            shell_material: ceramic,
            liquid_material: water,
            medium_material: air,
            exposure: 0.0,
            convection_coeff: h_conv,
            r_max_multiplier: r_max_mult,
            outer_boundary: OuterBoundary::Ambient,
        };
        let n_snaps_data = simulate_sphere_heat_numerical(&np, duration, snaps);
        let n_avg = liquid_average_temperature(&n_snaps_data);

        println!("\n{} (shell={}°C):", label, shell_temp);
        let mut max_grid_vs_num = 0.0_f64;
        let mut max_closed_vs_num = 0.0_f64;
        let mut max_grid_vs_closed = 0.0_f64;
        for i in 0..snaps {
            let d1 = (a_avg[i].1 - n_avg[i].1).abs();
            let d2 = (a_closed[i].1 - n_avg[i].1).abs();
            let d3 = (a_avg[i].1 - a_closed[i].1).abs();
            max_grid_vs_num = max_grid_vs_num.max(d1);
            max_closed_vs_num = max_closed_vs_num.max(d2);
            max_grid_vs_closed = max_grid_vs_closed.max(d3);
            if i % 10 == 0 {
                println!("  t={:5.1}min: grid={:.1} closed={:.1} numerical={:.1} | g-n={:.1} c-n={:.1} g-c={:.1}",
                    a_avg[i].0 / 60.0, a_avg[i].1, a_closed[i].1, n_avg[i].1,
                    a_avg[i].1 - n_avg[i].1, a_closed[i].1 - n_avg[i].1, a_avg[i].1 - a_closed[i].1);
            }
        }
        println!(
            "  MAX DIFFS: grid_vs_num={:.2} closed_vs_num={:.2} grid_vs_closed={:.2}",
            max_grid_vs_num, max_closed_vs_num, max_grid_vs_closed
        );
    }
}

/// Test the grid-based liquid average (as the panel uses it) vs closed-form.
/// The panel uses convert_snapshots + liquid_average_temperature on a 80-point
/// uniform grid, which can introduce artifacts at the interface.
#[test]
fn test_panel_grid_vs_closedform() {
    use super::convert_snapshots;
    use heat_equation_analytical::simulate_three_region_analytical;

    let water = HeatMaterial::water();
    let ceramic = HeatMaterial {
        emissivity: 0.0,
        ..HeatMaterial::ceramic()
    };
    let air = HeatMaterial::air();

    let volume_ml = 180.0;
    let inner_r_cm = (3.0 * volume_ml / (4.0 * std::f64::consts::PI)).cbrt();
    let wall_mm = 5.0;
    let outer_r_cm = inner_r_cm + wall_mm / 10.0;
    let r_max_mult = 10.0;

    let r_inner = inner_r_cm / 100.0;
    let r_outer = outer_r_cm / 100.0;
    let r_max = outer_r_cm * r_max_mult / 100.0;
    let display_r_max = 0.10; // 10cm, same as panel

    let params = ThreeRegionParams {
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
        t_shell_init: 22.0, // cold mug
        t_ambient: 22.0,
        num_modes: 50,
        h_conv: 10.0,
        h_inner: 0.0,
    };

    let modes = compute_eigenmodes(&params);
    let duration = 3600.0;
    let snaps = 100;

    // Closed-form liquid average (ground truth)
    let closed_avg = three_region_liquid_average_closed_form(&params, &modes, duration, snaps);

    // Grid-based liquid average (what the panel uses)
    let num_radial = 80;
    let analytical_snaps = simulate_three_region_analytical(
        &params,
        &modes,
        duration,
        snaps,
        num_radial,
        display_r_max,
    );
    let heat_snaps = convert_snapshots(analytical_snaps);
    let grid_avg = liquid_average_temperature(&heat_snaps);

    println!("\nPanel grid (80 pts) vs closed-form liquid average:");
    println!("  r_inner = {:.4}m", r_inner);
    let dr = display_r_max / (num_radial - 1) as f64;
    let inner_idx_grid = (r_inner / dr).ceil() as usize;
    println!(
        "  grid dr = {:.4}mm, inner_idx = {}, r[inner_idx] = {:.4}m (overshoot = {:.3}mm)",
        dr * 1000.0,
        inner_idx_grid,
        inner_idx_grid as f64 * dr,
        (inner_idx_grid as f64 * dr - r_inner) * 1000.0
    );

    let mut max_diff = 0.0_f64;
    for i in 0..closed_avg.len().min(grid_avg.len()) {
        let diff = (closed_avg[i].1 - grid_avg[i].1).abs();
        max_diff = max_diff.max(diff);
        if i % 10 == 0 {
            println!(
                "  t={:5.1}min: closed={:.2}°C, grid={:.2}°C, diff={:.2}°C",
                closed_avg[i].0 / 60.0,
                closed_avg[i].1,
                grid_avg[i].1,
                diff
            );
        }
    }
    println!("  max_diff = {:.2}°C", max_diff);

    // Also test with more radial points
    let num_radial_fine = 300;
    let analytical_snaps_fine = simulate_three_region_analytical(
        &params,
        &modes,
        duration,
        snaps,
        num_radial_fine,
        display_r_max,
    );
    let heat_snaps_fine = convert_snapshots(analytical_snaps_fine);
    let grid_avg_fine = liquid_average_temperature(&heat_snaps_fine);

    println!("\nPanel grid (300 pts) vs closed-form:");
    let dr_fine = display_r_max / (num_radial_fine - 1) as f64;
    let inner_idx_fine = (r_inner / dr_fine).ceil() as usize;
    println!(
        "  grid dr = {:.4}mm, inner_idx = {}, r[inner_idx] = {:.4}m (overshoot = {:.3}mm)",
        dr_fine * 1000.0,
        inner_idx_fine,
        inner_idx_fine as f64 * dr_fine,
        (inner_idx_fine as f64 * dr_fine - r_inner) * 1000.0
    );

    let mut max_diff_fine = 0.0_f64;
    for i in 0..closed_avg.len().min(grid_avg_fine.len()) {
        let diff = (closed_avg[i].1 - grid_avg_fine[i].1).abs();
        max_diff_fine = max_diff_fine.max(diff);
        if i % 10 == 0 {
            println!(
                "  t={:5.1}min: closed={:.2}°C, grid={:.2}°C, diff={:.2}°C",
                closed_avg[i].0 / 60.0,
                closed_avg[i].1,
                grid_avg_fine[i].1,
                diff
            );
        }
    }
    println!("  max_diff = {:.2}°C", max_diff_fine);
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
        h_inner: 0.0,
    };
    let modes = compute_eigenmodes(&analytical_params);
    let analytical_avg =
        three_region_liquid_average_closed_form(&analytical_params, &modes, 3600.0, 50);

    let numerical_params = SphereHeatParams {
        outer_radius_m: r_outer,
        inner_radius_m: r_inner,
        shell_temp_c: 22.0,
        liquid_temp_c: 95.0,
        room_temp_c: 22.0,
        shell_material: ceramic,
        liquid_material: water,
        medium_material: air,
        exposure: 0.0,
        convection_coeff: 10.0,
        r_max_multiplier: r_max_mult,
        outer_boundary: OuterBoundary::Ambient,
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

    println!(
        "Grid: {} points, inner_idx={}, dr at interface = {:.4}mm",
        n_points,
        inner_idx,
        dr_at_interface * 1000.0
    );
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
        println!(
            "  t={:5.1}min: analytical={:.2}°C  numerical={:.2}°C  diff={:+.2}°C",
            t_a / 60.0,
            temp_a,
            temp_n,
            diff
        );
    }

    // Print grid cells around the inner interface
    println!("\nGrid around r_inner ({:.4}mm):", r_inner * 1000.0);
    let snap = &num_snaps[0];
    let lo = inner_idx.saturating_sub(5);
    let hi = (inner_idx + 5).min(snap.radii.len() - 1);
    for i in lo..=hi {
        let dr = if i > 0 {
            snap.radii[i] - snap.radii[i - 1]
        } else {
            0.0
        };
        let marker = if i == inner_idx {
            " <-- inner_interface_idx"
        } else {
            ""
        };
        println!(
            "  [{}] r={:.4}mm  T={:.1}°C  dr={:.4}mm{}",
            i,
            snap.radii[i] * 1000.0,
            snap.temperatures[i],
            dr * 1000.0,
            marker
        );
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
// They are different idealisations of a mug, so their agreement is a
// sanity band, not a correctness check: both run convection only (no
// radiation, sealed lid), and each tolerance is the observed max + 2 °C,
// rounded up. Correctness is held by the within-family comparisons.

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
        emissivity: 0.0, // disable radiation for cross-model fairness
        exposure: 0.0,   // disable evaporation for cross-model fairness
    }
}

/// Heat-equation setup matching `lumped_params`: no radiation, sealed lid.
fn cross_family_heat_eq(
    volume_ml: f64,
    wall_mm: f64,
    liquid_temp: f64,
    shell_temp: f64,
    room_temp: f64,
    h_conv: f64,
) -> CommonSetup {
    HeatEqCase {
        volume_ml,
        wall_mm,
        liquid_temp,
        shell_temp,
        room_temp,
        h_conv,
        shell_emissivity: 0.0,
        exposure: 0.0,
    }
    .build()
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
    lumped_series: &[(f64, f64)],  // (time_minutes, temp_c)
    heat_eq_series: &[(f64, f64)], // (time_seconds, temp_c)
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
    println!("    max_diff={:.2}°C, rms={:.2}°C", max_diff, rms);
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
    let he = cross_family_heat_eq(volume, wall, liquid, room, room, h);
    let heat_eq_avg = relinearised_liquid_average(
        &he.analytical,
        ANALYTICAL_SEGMENTS,
        duration_min as f64 * 60.0,
        duration_min * 5,
    );

    let (max_diff, rms) = compare_lumped_vs_heat_eq(
        &lumped_coffee,
        &heat_eq_avg,
        duration_min as f64,
        "Analytical lumped vs analytical heat equation (cold mug)",
    );

    // Sanity band (cylinder vs sphere). Observed: 3.84 °C max, 2.96 °C RMS.
    assert!(
        max_diff < 6.0,
        "Analytical cold: max diff {:.1}°C exceeds 6°C",
        max_diff,
    );
    assert!(rms < 5.0, "Analytical cold: RMS {:.1}°C exceeds 5°C", rms,);

    // Both should start at the same temperature
    assert!(
        (lumped_coffee[0].1 - heat_eq_avg[0].1).abs() < 0.5,
        "Both should start at ~{}°C",
        liquid,
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

    let he = cross_family_heat_eq(volume, wall, liquid, shell, room, h);
    let heat_eq_avg = relinearised_liquid_average(
        &he.analytical,
        ANALYTICAL_SEGMENTS,
        duration_min as f64 * 60.0,
        duration_min * 5,
    );

    let (max_diff, rms) = compare_lumped_vs_heat_eq(
        &lumped_coffee,
        &heat_eq_avg,
        duration_min as f64,
        "Analytical lumped vs analytical heat equation (preheated mug)",
    );

    // Sanity band (cylinder vs sphere). Observed: 4.21 °C max, 3.06 °C RMS.
    assert!(
        max_diff < 7.0,
        "Analytical preheated: max diff {:.1}°C exceeds 7°C",
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
    let he = cross_family_heat_eq(volume, wall, liquid, room, room, h);
    let numerical_snaps =
        simulate_sphere_heat_numerical(&he.numerical, duration_min as f64 * 60.0, duration_min * 5);
    let heat_eq_avg = liquid_average_temperature(&numerical_snaps);

    let (max_diff, rms) = compare_lumped_vs_heat_eq(
        &lumped_coffee,
        &heat_eq_avg,
        duration_min as f64,
        "Numerical lumped vs numerical heat equation (cold mug)",
    );

    // Sanity band (cylinder vs sphere). Observed: 3.85 °C max, 2.96 °C RMS.
    assert!(
        max_diff < 6.0,
        "Numerical cold: max diff {:.1}°C exceeds 6°C",
        max_diff,
    );
    assert!(rms < 5.0, "Numerical cold: RMS {:.1}°C exceeds 5°C", rms,);

    assert!(
        (lumped_coffee[0].1 - heat_eq_avg[0].1).abs() < 0.5,
        "Both should start at ~{}°C",
        liquid,
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

    let he = cross_family_heat_eq(volume, wall, liquid, shell, room, h);
    let numerical_snaps =
        simulate_sphere_heat_numerical(&he.numerical, duration_min as f64 * 60.0, duration_min * 5);
    let heat_eq_avg = liquid_average_temperature(&numerical_snaps);

    let (max_diff, rms) = compare_lumped_vs_heat_eq(
        &lumped_coffee,
        &heat_eq_avg,
        duration_min as f64,
        "Numerical lumped vs numerical heat equation (preheated mug)",
    );

    // Sanity band (cylinder vs sphere). Observed: 4.22 °C max, 3.07 °C RMS.
    assert!(
        max_diff < 7.0,
        "Numerical preheated: max diff {:.1}°C exceeds 7°C",
        max_diff,
    );
    assert!(
        rms < 6.0,
        "Numerical preheated: RMS {:.1}°C exceeds 6°C",
        rms,
    );
}

/// Diagnostic: eigenvalue structure with h_inner > 0 (sink at the inner surface).
/// Prints the eigenvalues, amplitudes and volume integrals for several sink
/// strengths, as a fraction of the whole-sphere open-top secant coefficient.
#[test]
fn test_robin_bc_eigenvalue_diagnostic() {
    use heat_equation_analytical::{compute_eigenmodes, three_region_liquid_average_closed_form};

    let water = HeatMaterial::water();
    let ceramic = HeatMaterial::ceramic();

    let volume_ml = 180.0;
    let inner_r_cm = (3.0 * volume_ml / (4.0 * std::f64::consts::PI)).cbrt();
    let wall_mm = 5.0;
    let outer_r_cm = inner_r_cm + wall_mm / 10.0;
    let r_max_mult = 10.0;
    let h_conv = 10.0;
    let liquid_temp = 95.0;
    let room_temp = 22.0;

    let r_inner = inner_r_cm / 100.0;
    let r_outer = outer_r_cm / 100.0;
    let r_max = outer_r_cm * r_max_mult / 100.0;

    let h_open_top_secant =
        surface::open_top_secant_h(0.5 * (liquid_temp + room_temp), room_temp, h_conv);

    let base_params = ThreeRegionParams::for_mug(&ThreeRegionMug {
        r_inner,
        r_outer,
        r_max,
        shell: ceramic,
        liquid_temp_c: liquid_temp,
        shell_temp_c: room_temp,
        room_temp_c: room_temp,
        h_conv,
        exposure: 0.0,
        num_modes: 50,
    });

    println!("\n======================================================================");
    println!("ROBIN BC EIGENVALUE DIAGNOSTIC");
    println!("======================================================================");
    println!(
        "r_inner = {:.4}m, h_open_top_secant = {:.1} W/m²K",
        r_inner, h_open_top_secant
    );
    println!(
        "Bi_lumped = h_eff × r / k = {:.2}",
        h_open_top_secant * r_inner / water.conductivity
    );
    println!(
        "Lumped τ = ρcV/(hA) = {:.0}s",
        water.density * water.specific_heat * r_inner / (3.0 * h_open_top_secant)
    );

    for fraction in [0.0, 0.1, 0.2, 0.5] {
        let h_inner_val = fraction * h_open_top_secant;
        let label = format!("h_inner={h_inner_val:.1} ({:.0}%)", fraction * 100.0);
        let params = ThreeRegionParams {
            h_inner: h_inner_val,
            ..base_params.clone()
        };
        let modes = compute_eigenmodes(&params);

        println!("\n--- {} ---", label);
        println!("Found {} modes", modes.len());

        // Print first 10 eigenvalues with decay times
        for (i, m) in modes.iter().take(10).enumerate() {
            let tau = 1.0 / (m.lambda * m.lambda);
            let beta1 = m.lambda / params.alpha_1.sqrt();
            let vol_integral =
                m.coeff_a * super::heat_equation_analytical::integral_r_sin(beta1, 0.0, r_inner);
            println!(
                "  mode {:2}: λ={:.6}, τ={:8.1}s, amp={:10.3}, A={:8.4}, vol_int={:10.6}",
                i + 1,
                m.lambda,
                tau,
                m.amplitude,
                m.coeff_a,
                vol_integral
            );
        }

        // Check Parseval at t=0+
        let r_i3 = r_inner * r_inner * r_inner;
        let coeff = 3.0 / r_i3;
        let sum_at_t0: f64 = modes
            .iter()
            .map(|m| {
                let beta1 = m.lambda / params.alpha_1.sqrt();
                let vol_int = m.coeff_a
                    * super::heat_equation_analytical::integral_r_sin(beta1, 0.0, r_inner);
                m.amplitude * vol_int
            })
            .sum();
        let t_avg_t0 = params.t_ambient + coeff * sum_at_t0;
        println!(
            "  Sum at t=0+: {:.2}, T_avg(0+) = {:.1}°C (should be {:.0}°C)",
            sum_at_t0, t_avg_t0, liquid_temp
        );

        // Check at a few time points
        let avg = three_region_liquid_average_closed_form(&params, &modes, 3600.0, 100);
        println!(
            "  T_avg curve: t=0={:.1}, t=6min={:.1}, t=30min={:.1}, t=60min={:.1}",
            avg[0].1, avg[10].1, avg[50].1, avg[99].1
        );

        // Dense determinant scan to find ALL sign changes
        if h_inner_val > 0.0 {
            // Scan from very small λ up to 0.25
            for &(lo, hi, label) in &[
                (0.0001, 0.001, "[0.0001, 0.001]"),
                (0.001, 0.02, "[0.001, 0.02]"),
                (0.02, 0.05, "[0.02, 0.05]"),
            ] {
                let n_scan = 10000;
                let dl = (hi - lo) / n_scan as f64;
                let mut prev_det =
                    super::heat_equation_analytical::interface_determinant(lo, &params);
                let mut prev_lam = lo;
                let mut found = 0;
                for i in 1..=n_scan {
                    let lam = lo + i as f64 * dl;
                    let det = super::heat_equation_analytical::interface_determinant(lam, &params);
                    if prev_det.is_finite() && det.is_finite() && prev_det * det < 0.0 {
                        found += 1;
                        let tau = 1.0 / ((prev_lam + lam) / 2.0).powi(2);
                        println!(
                            "    {} sign change at λ≈{:.6} τ≈{:.0}s (det: {:.3e} → {:.3e})",
                            label,
                            (prev_lam + lam) / 2.0,
                            tau,
                            prev_det,
                            det
                        );
                    }
                    prev_det = det;
                    prev_lam = lam;
                }
                if found == 0 {
                    println!("    {} no sign changes", label);
                }
            }

            // Print determinant values at the old eigenvalue location
            println!("  Determinant at old λ=0.019:");
            for &lam in &[0.015, 0.017, 0.018, 0.0188, 0.019, 0.020, 0.022, 0.025] {
                let det = super::heat_equation_analytical::interface_determinant(lam, &params);
                println!("    det({:.4}) = {:.6e}", lam, det);
            }
        }
    }
}

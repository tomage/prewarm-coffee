//! Energy-conservation tests: a closed coffee + mug system with every loss
//! switched off must settle at the heat-capacity-weighted mean of its initial
//! temperatures, `T_eq = (C_coffee·T_coffee + C_mug·T_mug)/(C_coffee + C_mug)`.
//!
//! Each model is checked against its own heat capacities, so the tests catch
//! any solver that creates or destroys heat, independent of geometry choices.

use super::{
    liquid_average_temperature, shell_average_temperature, simulate_lumped_analytical,
    simulate_lumped_numerical, simulate_sphere_heat_numerical, HeatMaterial, MugMaterial,
    MugParameters, OuterBoundary, SphereHeatParams,
};

const VOLUME_ML: f64 = 350.0;
const WALL_MM: f64 = 5.0;
const COFFEE_C: f64 = 95.0;
const MUG_C: f64 = 22.0;
const TOLERANCE_C: f64 = 0.1;
/// The heat-equation solver conserves Σ C_i·T_i exactly, so its
/// capacity-weighted mean only drifts by round-off.
const MEAN_TOLERANCE_C: f64 = 1e-3;

fn equilibrium(c_coffee: f64, c_mug: f64) -> f64 {
    (c_coffee * COFFEE_C + c_mug * MUG_C) / (c_coffee + c_mug)
}

/// Ceramic mug with no convection, radiation or open-top losses.
fn closed_mug() -> MugParameters {
    MugParameters {
        material: MugMaterial::Ceramic,
        volume_ml: VOLUME_ML,
        wall_thickness_mm: WALL_MM,
        coffee_temp_c: COFFEE_C,
        room_temp_c: MUG_C,
        preheated: true,
        preheat_temp_c: MUG_C,
        h_conv: 0.0,
        emissivity: 0.0,
        exposure: 0.0,
    }
}

fn assert_lumped_settles(name: &str, final_coffee: f64, final_mug: f64) {
    let params = closed_mug();
    let t_eq = equilibrium(params.coffee_thermal_mass(), params.mug_thermal_mass());
    assert!(
        (final_coffee - t_eq).abs() < TOLERANCE_C,
        "{name}: coffee {final_coffee:.3}°C, expected {t_eq:.3}°C"
    );
    assert!(
        (final_mug - t_eq).abs() < TOLERANCE_C,
        "{name}: mug {final_mug:.3}°C, expected {t_eq:.3}°C"
    );
}

#[test]
fn lumped_numerical_conserves_energy() {
    let data = simulate_lumped_numerical(&closed_mug(), 120);
    let last = data.last().unwrap();
    assert_lumped_settles("lumped numerical", last.coffee_temp_c, last.mug_temp_c);
}

#[test]
fn lumped_analytical_conserves_energy() {
    let data = simulate_lumped_analytical(&closed_mug(), 120);
    let last = data.last().unwrap();
    assert_lumped_settles("lumped analytical", last.coffee_temp_c, last.mug_temp_c);
}

/// Closed sphere of coffee in a shell of `material`, insulated outside.
fn assert_heat_eq_settles(material: MugMaterial) {
    let volume_m3 = VOLUME_ML * 1e-6;
    let r_inner = (3.0 * volume_m3 / (4.0 * std::f64::consts::PI)).cbrt();
    let r_outer = r_inner + WALL_MM * 1e-3;
    let shell = HeatMaterial::from_mug(material, 0.0);
    let liquid = HeatMaterial::water();

    let params = SphereHeatParams {
        outer_radius_m: r_outer,
        inner_radius_m: r_inner,
        shell_temp_c: MUG_C,
        liquid_temp_c: COFFEE_C,
        room_temp_c: MUG_C,
        shell_material: shell,
        liquid_material: liquid,
        medium_material: HeatMaterial::air(),
        exposure: 0.0,
        convection_coeff: 0.0,
        r_max_multiplier: 1.0,
        outer_boundary: OuterBoundary::Insulated,
    };
    let snapshots = simulate_sphere_heat_numerical(&params, 3600.0, 10);

    // Heat capacities of the continuous spheres (the 4π/3 factor cancels).
    let c_liquid = liquid.density * liquid.specific_heat * r_inner.powi(3);
    let c_shell = shell.density * shell.specific_heat * (r_outer.powi(3) - r_inner.powi(3));
    let t_eq = equilibrium(c_liquid, c_shell);

    let t_liquid = liquid_average_temperature(&snapshots).last().unwrap().1;
    let t_shell = shell_average_temperature(&snapshots).last().unwrap().1;
    let t_mean = (c_liquid * t_liquid + c_shell * t_shell) / (c_liquid + c_shell);

    assert!(
        (t_mean - t_eq).abs() < MEAN_TOLERANCE_C,
        "{material:?}: capacity-weighted mean {t_mean:.6}°C, expected {t_eq:.6}°C"
    );
    assert!(
        (t_liquid - t_eq).abs() < TOLERANCE_C,
        "{material:?}: liquid average {t_liquid:.3}°C, expected {t_eq:.3}°C"
    );
    assert!(
        (t_shell - t_eq).abs() < TOLERANCE_C,
        "{material:?}: shell average {t_shell:.3}°C, expected {t_eq:.3}°C"
    );
}

#[test]
fn heat_eq_numerical_conserves_energy() {
    assert_heat_eq_settles(MugMaterial::Ceramic);
}

/// Steel's conductivity is ~16× ceramic's, so the interface faces carry
/// very different conductances on each side of the shared node.
#[test]
fn heat_eq_numerical_steel_conserves_energy() {
    assert_heat_eq_settles(MugMaterial::StainlessSteel);
}

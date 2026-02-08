mod heat_equation_analytical;
mod heat_equation_numerical;
mod lumped_analytical;
mod lumped_numerical;
mod materials;
mod metrics;
#[cfg(test)]
mod cross_solver_tests;
#[cfg(test)]
mod perf_tests;

pub use heat_equation_analytical::{
    compute_eigenmodes, compute_exposed_h_eff, exposed_liquid_average,
    simulate_three_region_analytical, Eigenmode, ExposedLiquidParams, ThreeRegionParams,
    ThreeRegionSnapshot,
};
pub use heat_equation_numerical::{
    ball_average_temperature, ball_center_temperature, liquid_average_temperature,
    liquid_center_temperature, liquid_near_wall_temperature, shell_average_temperature,
    simulate_sphere_heat_numerical, HeatMaterial, HeatSnapshot, SphereHeatParams,
    K_AIR, EXPOSED_CONVECTION_MULTIPLIER,
};
pub use lumped_analytical::simulate_lumped_analytical;
pub use lumped_numerical::{
    simulate_lumped_numerical, CoolingDataPoint, MAX_SIMULATION_MINUTES,
};
pub use materials::{MugMaterial, MugParameters};
pub use metrics::{
    average_during_window, average_while_above, compute_comparison,
    find_crossing_time_series, seconds_to_minutes, ComparisonData,
};

/// Geometric scaling factor for mapping lid-opening exposure from a cylinder
/// (lumped model) to a sphere (heat equation model).
///
/// The lumped model uses a cylinder with h = 2r, where the lid opening is a
/// fraction of the top circle (area = πr_cyl²).  The heat equation uses a
/// sphere, where inner-surface heat loss is distributed over 4πr_sphere².
///
/// For the same volume V:
///   r_cyl = (V/(2π))^(1/3),  r_sphere = (3V/(4π))^(1/3)
///   ratio = πr_cyl² / (4πr_sphere²) = (2/3)^(2/3) / 4 ≈ 0.191
///
/// Without this factor, the sphere model applies evaporation to ~5× more
/// surface area than the cylinder, causing the numerical Heat Eq to cool
/// significantly faster than the lumped models.
pub const LID_TO_SPHERE_AREA_RATIO: f64 = 0.191;

/// Apply exposure correction to liquid temperatures via multiplicative decay.
///
/// Models heat loss through the lid opening as a decoupled exponential correction:
///   T(r,t) = T_amb + [T_3region(r,t) - T_amb] × exp(-t/τ)
///
/// Uses only convection + radiation (no evaporation) for h_eff.  The full
/// h_eff_exposed (~340 W/m²K) is dominated by evaporative cooling which
/// cannot be accurately captured by a uniform multiplicative decay on the
/// entire liquid.  Convection + radiation alone (~18 W/m²K) gives a mild,
/// physically reasonable correction.
///
/// Time constant: τ = τ_surface + τ_diffusion = ρcr/(3h) + r²/(π²α)
///
/// Returns the snapshots unchanged if `exposure <= 0`.
pub fn apply_exposure_correction(
    snapshots: Vec<HeatSnapshot>,
    exposure: f64,
    room_temp: f64,
    h_conv: f64,
    inner_radius_m: f64,
    liquid_temp_init: f64,
) -> Vec<HeatSnapshot> {
    if exposure <= 0.0 {
        return snapshots;
    }
    let water = HeatMaterial::water();

    // Convection + linearized radiation from the water surface (no evaporation).
    let sigma: f64 = 5.670_374e-8;
    let water_emissivity: f64 = 0.96;
    let t_avg_k = 0.5 * (liquid_temp_init + room_temp) + 273.15;
    let h_rad_water = 4.0 * water_emissivity * sigma * t_avg_k.powi(3);
    let h_eff = h_conv + h_rad_water;
    let h_inner = exposure * h_eff;

    // Time constant combining surface resistance and internal diffusion:
    //   τ = τ_surface + τ_diffusion  =  ρcr/(3h) + r²/(π²α)
    let rho_c = water.density * water.specific_heat;
    let tau_surface = rho_c * inner_radius_m / (3.0 * h_inner);
    let tau_diffusion = inner_radius_m * inner_radius_m
        / (std::f64::consts::PI.powi(2) * water.diffusivity());
    let tau = tau_surface + tau_diffusion;

    snapshots
        .into_iter()
        .map(|mut snap| {
            let decay = (-snap.time_s / tau).exp();
            let inner_idx = snap.inner_interface_idx;
            for i in 0..=inner_idx {
                snap.temperatures[i] =
                    room_temp + (snap.temperatures[i] - room_temp) * decay;
            }
            snap
        })
        .collect()
}

/// Convert analytical ThreeRegionSnapshot to HeatSnapshot so downstream code works unchanged.
pub fn convert_snapshots(snaps: Vec<ThreeRegionSnapshot>) -> Vec<HeatSnapshot> {
    snaps
        .into_iter()
        .map(|s| {
            let inner_idx = s
                .radii
                .iter()
                .position(|&r| r >= s.r_inner)
                .unwrap_or(0);
            let outer_idx = s
                .radii
                .iter()
                .position(|&r| r >= s.r_outer)
                .unwrap_or(s.radii.len() - 1);
            HeatSnapshot {
                time_s: s.time_s,
                radii: s.radii,
                temperatures: s.temperatures,
                inner_interface_idx: inner_idx,
                outer_interface_idx: outer_idx,
            }
        })
        .collect()
}

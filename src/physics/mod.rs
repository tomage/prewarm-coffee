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

/// Blend liquid-region temperatures with the exposed-sphere model.
///
/// For each snapshot, each liquid grid point (index <= inner_interface_idx) is blended:
/// `T = (1 - exposure) * T_3region(r,t) + exposure * T_exposed_avg(t)`
///
/// Returns the snapshots unchanged if `exposure <= 0`.
pub fn blend_exposure(
    snapshots: Vec<HeatSnapshot>,
    exposure: f64,
    liquid_temp_init: f64,
    room_temp: f64,
    h_conv: f64,
    inner_radius_m: f64,
    num_modes: usize,
    duration_s: f64,
) -> Vec<HeatSnapshot> {
    if exposure <= 0.0 {
        return snapshots;
    }
    let water = HeatMaterial::water();
    let h_eff = compute_exposed_h_eff(liquid_temp_init, room_temp, h_conv);
    let exposed_params = ExposedLiquidParams {
        radius: inner_radius_m,
        alpha: water.diffusivity(),
        k: water.conductivity,
        h_eff,
        t_init: liquid_temp_init,
        t_ambient: room_temp,
    };
    let exposed_avg = exposed_liquid_average(&exposed_params, num_modes, duration_s, snapshots.len());

    snapshots
        .into_iter()
        .zip(exposed_avg.iter())
        .map(|(mut snap, &(_, t_exp))| {
            let inner_idx = snap.inner_interface_idx;
            for i in 0..=inner_idx {
                snap.temperatures[i] =
                    (1.0 - exposure) * snap.temperatures[i] + exposure * t_exp;
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

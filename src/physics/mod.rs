#[cfg(test)]
mod claims_tests;
#[cfg(test)]
mod conservation_tests;
pub mod constants;
#[cfg(test)]
mod cross_solver_tests;
#[cfg(test)]
mod exposure_tests;
mod heat_equation_analytical;
mod heat_equation_numerical;
mod lumped_analytical;
mod lumped_numerical;
mod materials;
mod metrics;
#[cfg(test)]
mod perf_tests;
#[cfg(test)]
mod solver_grid_tests;
pub mod surface;

pub use heat_equation_analytical::{
    compute_eigenmodes, simulate_three_region_relinearised, Eigenmode, ThreeRegionMug,
    ThreeRegionParams, ThreeRegionSnapshot, ANALYTICAL_SEGMENTS,
};
pub use heat_equation_numerical::{
    liquid_average_temperature, liquid_center_temperature, liquid_near_wall_temperature,
    shell_average_temperature, simulate_sphere_heat_numerical, HeatMaterial, HeatSnapshot,
    OuterBoundary, SphereHeatParams,
};
pub use lumped_analytical::simulate_lumped_analytical;
pub use lumped_numerical::{simulate_lumped_numerical, CoolingDataPoint, MAX_SIMULATION_MINUTES};
pub use materials::{MugMaterial, MugParameters};
pub use metrics::{
    average_during_window, average_while_above, compute_comparison, extra_time,
    find_crossing_time_series, seconds_to_minutes, series_up_to, ComparisonData, ExtraTime,
};

/// Convert analytical ThreeRegionSnapshot to HeatSnapshot so downstream code works unchanged.
pub fn convert_snapshots(snaps: Vec<ThreeRegionSnapshot>) -> Vec<HeatSnapshot> {
    snaps
        .into_iter()
        .map(|s| {
            let inner_idx = s.radii.iter().position(|&r| r >= s.r_inner).unwrap_or(0);
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

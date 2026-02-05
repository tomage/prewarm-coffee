//! Simulation of the heat equation for a sphere in an infinite medium
//!
//! Models temperature distribution T(r, t) for a shell (ceramic, steel, etc.)
//! containing liquid, surrounded by air, solving the spherical heat equation:
//!
//! ∂T/∂t = α (∂²T/∂r² + (2/r)(∂T/∂r))
//!
//! where α = k/(ρc) is the thermal diffusivity.
//!
//! Features:
//! - Three-region model: liquid core, shell, surrounding air
//! - Interface conditions using harmonic mean of diffusivities
//! - Optional evaporative cooling at the liquid surface
//! - Convective heat transfer at the outer boundary

/// Stefan-Boltzmann constant W/(m²·K⁴)
const STEFAN_BOLTZMANN: f64 = 5.670_374e-8;

// Evaporation constants (same as lumped model)
const LATENT_HEAT_VAPORIZATION: f64 = 2.26e6; // J/kg
const R_VAPOR: f64 = 461.5; // J/(kg·K)
const RHO_AIR: f64 = 1.18; // kg/m³
const CP_AIR: f64 = 1005.0; // J/(kg·K)
const RELATIVE_HUMIDITY: f64 = 0.50;
const EVAP_ENHANCEMENT: f64 = 2.5;

fn saturation_vapor_pressure(t_celsius: f64) -> f64 {
    610.78 * (17.27 * t_celsius / (t_celsius + 237.3)).exp()
}

/// Material properties for the heat equation
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HeatMaterial {
    /// Thermal conductivity W/(m·K)
    pub conductivity: f64,
    /// Density kg/m³
    pub density: f64,
    /// Specific heat J/(kg·K)
    pub specific_heat: f64,
    /// Emissivity (0–1)
    pub emissivity: f64,
}

/// Physical thermal conductivity of air at ~25°C in W/(m·K).
///
/// Used for boundary-layer profile rendering in the profile chart and for
/// convection calculations at the mug outer surface.
pub const K_AIR: f64 = 0.026;

/// Convection coefficient multiplier for exposed (no shell) simulations.
pub const EXPOSED_CONVECTION_MULTIPLIER: f64 = 1.5;

impl HeatMaterial {
    /// Thermal diffusivity α = k/(ρc) in m²/s
    pub fn diffusivity(&self) -> f64 {
        self.conductivity / (self.density * self.specific_heat)
    }

    /// Steel properties
    pub fn steel() -> Self {
        Self {
            conductivity: 50.0,      // W/(m·K)
            density: 7800.0,         // kg/m³
            specific_heat: 500.0,    // J/(kg·K)
            emissivity: 0.10,
        }
    }

    /// Ceramic properties
    pub fn ceramic() -> Self {
        Self {
            conductivity: 1.0,       // W/(m·K)
            density: 2400.0,         // kg/m³
            specific_heat: 880.0,    // J/(kg·K)
            emissivity: 0.90,
        }
    }

    /// Air properties (at ~25°C)
    pub fn air() -> Self {
        Self {
            conductivity: 0.026,     // W/(m·K)
            density: 1.18,           // kg/m³
            specific_heat: 1005.0,   // J/(kg·K)
            emissivity: 1.0,         // not physically meaningful for air medium
        }
    }

    /// Water properties (effective conductivity includes natural convection)
    ///
    /// Molecular conductivity of water is 0.6 W/(m·K), but hot coffee has
    /// vigorous buoyancy-driven convection (Ra ~ 10⁹) that keeps the bulk
    /// nearly well-mixed. Using an effective conductivity of ~10× molecular
    /// approximates this mixing without a full Navier-Stokes solve.
    pub fn water() -> Self {
        Self {
            conductivity: 6.0,       // W/(m·K) effective (≈10× molecular)
            density: 1000.0,         // kg/m³
            specific_heat: 4186.0,   // J/(kg·K)
            emissivity: 0.96,
        }
    }

    /// Coffee (similar to water)
    pub fn coffee() -> Self {
        Self::water()
    }
}

/// Parameters for the sphere heat equation simulation
#[derive(Debug, Clone, PartialEq)]
pub struct SphereHeatParams {
    /// Outer radius of steel shell in meters
    pub outer_radius_m: f64,
    /// Inner radius of steel shell (liquid cavity) in meters
    pub inner_radius_m: f64,
    /// Initial steel shell temperature in °C
    pub shell_temp_c: f64,
    /// Initial liquid temperature in °C
    pub liquid_temp_c: f64,
    /// Room/air temperature in °C
    pub room_temp_c: f64,
    /// Shell material
    pub shell_material: HeatMaterial,
    /// Liquid material (inside shell)
    pub liquid_material: HeatMaterial,
    /// Surrounding medium material
    pub medium_material: HeatMaterial,
    /// Evaporation factor (0.0 = none, 1.0 = typical for water at ~80°C)
    /// Adds extra heat loss at liquid surface due to evaporative cooling
    pub evaporation_factor: f64,
    /// Convective heat transfer coefficient at outer boundary W/(m²·K)
    /// Typical values: 5-25 for natural convection in air, 10-200 for forced convection
    /// Set to 0.0 to use pure conduction (Dirichlet boundary at r_max)
    pub convection_coeff: f64,
    /// r_max = outer_radius_m * r_max_multiplier. Default 3.0.
    pub r_max_multiplier: f64,
}

impl Default for SphereHeatParams {
    fn default() -> Self {
        Self {
            outer_radius_m: 0.05,  // 5 cm outer radius
            inner_radius_m: 0.04,  // 4 cm inner radius (1 cm thick shell)
            shell_temp_c: 100.0,   // 100°C shell
            liquid_temp_c: 200.0,  // 200°C liquid (hotter than shell)
            room_temp_c: 20.0,     // 20°C room
            shell_material: HeatMaterial::steel(),
            liquid_material: HeatMaterial::water(),
            medium_material: HeatMaterial::air(),
            evaporation_factor: 0.0, // No evaporation by default
            convection_coeff: 10.0,  // Natural convection in air ~10 W/(m²·K)
            r_max_multiplier: 3.0,
        }
    }
}

/// Result of the heat equation simulation at a single time
#[derive(Debug, Clone, PartialEq)]
pub struct HeatSnapshot {
    /// Time in seconds
    pub time_s: f64,
    /// Radial positions in meters
    pub radii: Vec<f64>,
    /// Temperature at each radial position in °C
    pub temperatures: Vec<f64>,
    /// Index where liquid/shell interface is
    pub inner_interface_idx: usize,
    /// Index where shell/air interface is
    pub outer_interface_idx: usize,
}

/// Solve tridiagonal system Ax = rhs in-place using the Thomas algorithm.
///
/// `sub[i]`, `dia[i]`, `sup[i]` are the sub/main/super diagonal coefficients.
/// `rhs` is overwritten with the solution on output.
/// `dia` is modified during forward elimination.
fn thomas_solve(sub: &[f64], dia: &mut [f64], sup: &[f64], rhs: &mut [f64]) {
    let n = dia.len();
    // Forward elimination
    for i in 1..n {
        let w = sub[i] / dia[i - 1];
        dia[i] -= w * sup[i - 1];
        rhs[i] -= w * rhs[i - 1];
    }
    // Back substitution
    rhs[n - 1] /= dia[n - 1];
    for i in (0..n - 1).rev() {
        rhs[i] = (rhs[i] - sup[i] * rhs[i + 1]) / dia[i];
    }
}

/// Simulate the heat equation for a steel shell with liquid inside, surrounded by air
///
/// Uses finite difference method with Crank-Nicolson time stepping.
/// Three regions: liquid (0 to r_inner), shell (r_inner to r_outer), air (r_outer to infinity)
/// Returns snapshots of temperature distribution at specified times.
pub fn simulate_sphere_heat_numerical(
    params: &SphereHeatParams,
    duration_s: f64,
    num_snapshots: usize,
) -> Vec<HeatSnapshot> {
    // Spatial discretization with adaptive grid
    // Uses piecewise-uniform grids: coarse in bulk, fine near interfaces
    // Scale resolution inversely with duration to keep total work bounded
    // Domain size from multiplier. Convection modeled via Robin BC at r_max.
    let r_max = params.outer_radius_m * params.r_max_multiplier;
    let r_inner = params.inner_radius_m;
    let r_outer = params.outer_radius_m;
    let shell_thickness = r_outer - r_inner;

    // Grid layout (4 zones):
    //   [0, liquid_coarse_end]  — bulk liquid, uniform coarse
    //   (liquid_coarse_end, r_inner] — geometric refinement toward r_inner
    //   [r_inner, r_outer]     — shell, uniform fine
    //   (r_outer, r_max]       — air, geometric coarse
    //
    // The geometric refinement concentrates resolution at the liquid-shell
    // interface where temperature discontinuities cause the largest error.
    let scale = (120.0_f64 / duration_s.max(120.0)).cbrt().max(0.7);
    let n_bulk_liquid = (40.0 * scale).round().max(25.0) as usize;
    let n_shell = (80.0 * scale).round().max(50.0) as usize;

    // dr_shell: uniform spacing in the shell
    let dr_shell = shell_thickness / n_shell as f64;

    // Geometric refinement zone on the liquid side of r_inner.
    // Finest cell at r_inner matches shell spacing (dr_shell) for a smooth
    // transition across the interface. Cells grow by `growth` moving away
    // from the interface into the liquid bulk.
    let n_refine = (30.0 * scale).round().max(20.0) as usize;
    let growth = 1.12_f64;
    let dr_min = dr_shell;
    let refine_len = dr_min * (growth.powi(n_refine as i32) - 1.0) / (growth - 1.0);
    let liquid_coarse_end = (r_inner - refine_len).max(0.0);

    let mut radii = Vec::with_capacity(n_bulk_liquid + n_refine + n_shell + 40);

    // Zone 1: Bulk liquid [0, liquid_coarse_end] — uniform coarse
    if liquid_coarse_end > 1e-10 {
        for i in 0..n_bulk_liquid {
            radii.push(liquid_coarse_end * i as f64 / (n_bulk_liquid - 1) as f64);
        }
    } else {
        radii.push(0.0);
    }

    // Zone 2: Geometric refinement (liquid_coarse_end, r_inner] — finest cells
    // at r_inner, growing toward bulk liquid. Build from r_inner backward.
    {
        let mut refine_pts = Vec::with_capacity(n_refine);
        let mut r = r_inner;
        let mut dr = dr_min;
        for _ in 0..n_refine {
            refine_pts.push(r);
            r -= dr;
            dr *= growth;
        }
        refine_pts.reverse();
        let last_bulk = *radii.last().unwrap_or(&0.0);
        for &rp in &refine_pts {
            if rp > last_bulk + 1e-12 {
                radii.push(rp);
            }
        }
    }

    // Zone 3: Shell [r_inner, r_outer] — uniform fine
    for i in 1..=n_shell {
        radii.push(r_inner + shell_thickness * i as f64 / n_shell as f64);
    }

    // Zone 4: Air (r_outer, r_max] — geometric coarse
    let air_len = r_max - r_outer;
    if air_len > 1e-10 {
        let dr_start = dr_shell * 10.0;
        let growth = 1.4_f64;
        let mut r = r_outer + dr_start;
        while r < r_max {
            radii.push(r);
            let dr = (r - radii[radii.len() - 2]) * growth;
            r += dr;
        }
        if let Some(&last) = radii.last() {
            if (last - r_max).abs() > 1e-10 {
                radii.push(r_max);
            }
        }
    } else if let Some(&last) = radii.last() {
        if (last - r_max).abs() > 1e-10 {
            radii.push(r_max);
        }
    }

    let n_points = radii.len();

    // Find interface indices (closest grid point to each interface)
    let inner_interface_idx = radii.iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| (*a - r_inner).abs().partial_cmp(&(*b - r_inner).abs()).unwrap())
        .unwrap().0;
    let outer_interface_idx = radii.iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| (*a - r_outer).abs().partial_cmp(&(*b - r_outer).abs()).unwrap())
        .unwrap().0;

    // Initialize temperature field - three regions
    let mut temp: Vec<f64> = radii.iter().map(|&r| {
        if r <= params.inner_radius_m {
            params.liquid_temp_c
        } else if r <= params.outer_radius_m {
            params.shell_temp_c
        } else {
            params.room_temp_c
        }
    }).collect();

    // Pre-compute material properties at each grid point
    // Using conservative form: ρc ∂T/∂t = (1/r²) ∂/∂r(r² k ∂T/∂r)
    // This properly enforces heat flux continuity k₁∂T/∂r = k₂∂T/∂r at interfaces
    let k_arr: Vec<f64> = radii.iter().map(|&r| {
        if r <= params.inner_radius_m {
            params.liquid_material.conductivity
        } else if r <= params.outer_radius_m {
            params.shell_material.conductivity
        } else {
            params.medium_material.conductivity
        }
    }).collect();

    let rho_c_arr: Vec<f64> = radii.iter().map(|&r| {
        if r <= params.inner_radius_m {
            params.liquid_material.density * params.liquid_material.specific_heat
        } else if r <= params.outer_radius_m {
            params.shell_material.density * params.shell_material.specific_heat
        } else {
            params.medium_material.density * params.medium_material.specific_heat
        }
    }).collect();

    // Pre-compute spacing between adjacent points
    let dr_arr: Vec<f64> = (0..n_points - 1).map(|i| radii[i + 1] - radii[i]).collect();

    // Half-point conductivities: harmonic mean of k for flux continuity
    let k_half: Vec<f64> = (0..n_points - 1).map(|i| {
        2.0 * k_arr[i] * k_arr[i + 1] / (k_arr[i] + k_arr[i + 1])
    }).collect();

    // Half-point radii
    let r_half: Vec<f64> = (0..n_points - 1).map(|i| {
        0.5 * (radii[i] + radii[i + 1])
    }).collect();

    // Pre-compute tridiagonal operator coefficients for the diffusion operator
    // L(T)_i = lower[i]*T_{i-1} + diag[i]*T_i + upper[i]*T_{i+1} + source[i]
    let mut lower = vec![0.0; n_points];
    let mut upper = vec![0.0; n_points];
    let mut diag_op = vec![0.0; n_points];
    let mut source = vec![0.0; n_points];

    // Center (i=0): L'Hôpital's rule for spherical Laplacian at r=0
    // lim_{r→0} (1/r²)∂/∂r(r²k∂T/∂r) = 3k ∂²T/∂r² ≈ 6k(T₁-T₀)/dr₁²
    let dr1 = dr_arr[0];
    upper[0] = 6.0 * k_arr[0] / (rho_c_arr[0] * dr1 * dr1);
    diag_op[0] = -upper[0];

    // Interior points (1..n_points-1)
    for i in 1..n_points - 1 {
        let r = radii[i];
        let dr_right = dr_arr[i];
        let dr_left = dr_arr[i - 1];
        let dr_avg = 0.5 * (dr_left + dr_right);
        let r2 = r * r;
        let denom = rho_c_arr[i] * r2 * dr_avg;

        lower[i] = r_half[i - 1] * r_half[i - 1] * k_half[i - 1] / (dr_left * denom);
        upper[i] = r_half[i] * r_half[i] * k_half[i] / (dr_right * denom);
        diag_op[i] = -(lower[i] + upper[i]);

        // Convection at outer interface: implicit treatment (linear in T)
        if i == outer_interface_idx && params.convection_coeff > 0.0 {
            let conv_coeff = params.convection_coeff / (rho_c_arr[i] * dr_avg);
            diag_op[i] -= conv_coeff;
            source[i] += conv_coeff * params.room_temp_c;
        }
    }

    // Precompute dr_avg at the outer interface for per-step radiation updates
    let dr_avg_outer = {
        let dr_left = dr_arr[outer_interface_idx - 1];
        let dr_right = dr_arr[outer_interface_idx.min(dr_arr.len() - 1)];
        0.5 * (dr_left + dr_right)
    };
    // Store base operator values at outer interface (without radiation)
    let base_diag_outer = diag_op[outer_interface_idx];
    let base_source_outer = source[outer_interface_idx];

    // Precompute dr_avg at the inner interface for per-step evaporation updates
    let dr_avg_inner = {
        let dr_left = dr_arr[inner_interface_idx.saturating_sub(1)];
        let dr_right = dr_arr[inner_interface_idx.min(dr_arr.len() - 1)];
        0.5 * (dr_left + dr_right)
    };
    let base_diag_inner = diag_op[inner_interface_idx];
    let base_source_inner = source[inner_interface_idx];

    // Time stepping: CN is unconditionally stable, so dt is chosen for accuracy.
    // Use ~500 steps for the simulation, capped at snapshot_interval and 5s.
    let snapshot_interval = duration_s / (num_snapshots - 1) as f64;
    let dt = (duration_s / 500.0).min(0.5).min(snapshot_interval);

    let mut snapshots = Vec::with_capacity(num_snapshots);

    // Record initial state
    snapshots.push(HeatSnapshot {
        time_s: 0.0,
        radii: radii.clone(),
        temperatures: temp.clone(),
        inner_interface_idx,
        outer_interface_idx,
    });

    let mut time = 0.0;
    let mut next_snapshot = snapshot_interval;

    // Scratch buffers for the Thomas solver
    let mut rhs = vec![0.0; n_points];
    let mut dia_lhs = vec![0.0; n_points];
    let mut sub_lhs = vec![0.0; n_points];
    let mut sup_lhs = vec![0.0; n_points];

    while snapshots.len() < num_snapshots {
        let half_dt = 0.5 * dt;

        // Linearized radiation at mug outer surface: h_rad = 4·ε·σ·T_m³
        // Re-linearized each step for accuracy, treated implicitly like convection.
        {
            let t_k = temp[outer_interface_idx] + 273.15;
            let t_room_k = params.room_temp_c + 273.15;
            let t_m = 0.5 * (t_k + t_room_k);
            let h_rad = 4.0 * params.shell_material.emissivity * STEFAN_BOLTZMANN * t_m.powi(3);
            let rad_coeff = h_rad / (rho_c_arr[outer_interface_idx] * dr_avg_outer);
            diag_op[outer_interface_idx] = base_diag_outer - rad_coeff;
            source[outer_interface_idx] = base_source_outer + rad_coeff * params.room_temp_c;
        }

        // Coffee surface heat loss through the open top (scaled by evaporation_factor
        // = lid opening fraction). Includes evaporation, convection, and radiation —
        // all treated implicitly for stability.
        if params.evaporation_factor > 0.0 {
            let t_surf = temp[inner_interface_idx];
            let t_room = params.room_temp_c;
            let t_surf_k = t_surf + 273.15;
            let t_room_k = t_room + 273.15;
            let exposure = params.evaporation_factor;
            let rho_c_inner = rho_c_arr[inner_interface_idx];
            let mut total_coeff = 0.0; // accumulated h_eff / (ρc · dr) [1/s]

            // 1. Evaporation (Lewis analogy, same physics as lumped model)
            let h_mass = params.convection_coeff / (RHO_AIR * CP_AIR);
            let rho_v_surf = saturation_vapor_pressure(t_surf) / (R_VAPOR * t_surf_k);
            let rho_v_ambient = RELATIVE_HUMIDITY * saturation_vapor_pressure(t_room) / (R_VAPOR * t_room_k);
            let delta_t = t_surf - t_room;
            if rho_v_surf > rho_v_ambient && delta_t > 1.0 {
                let q_evap = EVAP_ENHANCEMENT * exposure
                    * h_mass * LATENT_HEAT_VAPORIZATION * (rho_v_surf - rho_v_ambient);
                total_coeff += q_evap / (delta_t * rho_c_inner * dr_avg_inner);
            }

            // 2. Direct convection from coffee surface through opening
            let h_conv_inner = exposure * params.convection_coeff;
            total_coeff += h_conv_inner / (rho_c_inner * dr_avg_inner);

            // 3. Radiation from coffee surface through opening
            let t_m = 0.5 * (t_surf_k + t_room_k);
            let h_rad_inner = exposure * 4.0 * 0.96 * STEFAN_BOLTZMANN * t_m.powi(3);
            total_coeff += h_rad_inner / (rho_c_inner * dr_avg_inner);

            diag_op[inner_interface_idx] = base_diag_inner - total_coeff;
            source[inner_interface_idx] = base_source_inner + total_coeff * t_room;
        }

        // Build RHS: (I + dt/2 * L) T^n + dt * source
        // Center (i=0)
        rhs[0] = (1.0 + half_dt * diag_op[0]) * temp[0]
                + half_dt * upper[0] * temp[1]
                + dt * source[0];

        // Interior points
        for i in 1..n_points - 1 {
            rhs[i] = half_dt * lower[i] * temp[i - 1]
                   + (1.0 + half_dt * diag_op[i]) * temp[i]
                   + half_dt * upper[i] * temp[i + 1]
                   + dt * source[i];
        }

        // Boundary at r_max: Dirichlet BC
        rhs[n_points - 1] = params.room_temp_c;

        // Build LHS tridiagonal: (I - dt/2 * L)
        // Center
        sub_lhs[0] = 0.0;
        dia_lhs[0] = 1.0 - half_dt * diag_op[0];
        sup_lhs[0] = -half_dt * upper[0];

        for i in 1..n_points - 1 {
            sub_lhs[i] = -half_dt * lower[i];
            dia_lhs[i] = 1.0 - half_dt * diag_op[i];
            sup_lhs[i] = -half_dt * upper[i];
        }

        // Boundary at r_max: identity row
        sub_lhs[n_points - 1] = 0.0;
        dia_lhs[n_points - 1] = 1.0;
        sup_lhs[n_points - 1] = 0.0;

        // Solve tridiagonal system
        thomas_solve(&sub_lhs, &mut dia_lhs, &sup_lhs, &mut rhs);

        // Copy solution back to temp
        temp.copy_from_slice(&rhs);
        time += dt;

        // Record snapshot if it's time
        if time >= next_snapshot - 1e-10 {
            snapshots.push(HeatSnapshot {
                time_s: time,
                radii: radii.clone(),
                temperatures: temp.clone(),
                inner_interface_idx,
                outer_interface_idx,
            });
            next_snapshot += snapshot_interval;
        }
    }

    snapshots
}

/// Get temperature at center of liquid over time
pub fn liquid_center_temperature(snapshots: &[HeatSnapshot]) -> Vec<(f64, f64)> {
    snapshots.iter().map(|s| (s.time_s, s.temperatures[0])).collect()
}

/// Get average liquid temperature over time
pub fn liquid_average_temperature(snapshots: &[HeatSnapshot]) -> Vec<(f64, f64)> {
    snapshots.iter().map(|s| {
        let interface = s.inner_interface_idx;
        // Volume-weighted average (4πr²dr weighting)
        let mut sum = 0.0;
        let mut volume = 0.0;
        for i in 0..=interface {
            let r = s.radii[i];
            let dr = if i == 0 {
                s.radii[1] - s.radii[0]
            } else {
                s.radii[i] - s.radii[i - 1]
            };
            // Use spherical shell volume: V = 4π/3 * (r_outer³ - r_inner³)
            // For point i, the shell extends from r - dr/2 to r + dr/2
            let r_inner = (r - dr / 2.0).max(0.0);
            let r_outer = r + dr / 2.0;
            let shell_volume = 4.0 / 3.0 * std::f64::consts::PI * (r_outer.powi(3) - r_inner.powi(3));
            sum += s.temperatures[i] * shell_volume;
            volume += shell_volume;
        }
        (s.time_s, if volume > 0.0 { sum / volume } else { s.temperatures[0] })
    }).collect()
}

/// Get volume-weighted liquid near-wall average temperature over time.
/// `wall_fraction` is the fraction of liquid radius from the wall inward
/// (e.g. 0.10 = outer 10% of the liquid radius).
pub fn liquid_near_wall_temperature(snapshots: &[HeatSnapshot], wall_fraction: f64) -> Vec<(f64, f64)> {
    snapshots.iter().map(|s| {
        let interface = s.inner_interface_idx;
        let r_inner = s.radii[interface];
        let r_start = r_inner * (1.0 - wall_fraction);
        let mut sum = 0.0;
        let mut volume = 0.0;
        for i in 0..=interface {
            let r = s.radii[i];
            if r < r_start {
                continue;
            }
            let dr = if i == 0 {
                s.radii[1] - s.radii[0]
            } else {
                s.radii[i] - s.radii[i - 1]
            };
            let r_lo = (r - dr / 2.0).max(r_start);
            let r_hi = (r + dr / 2.0).min(r_inner);
            let shell_volume = 4.0 / 3.0 * std::f64::consts::PI * (r_hi.powi(3) - r_lo.powi(3));
            sum += s.temperatures[i] * shell_volume;
            volume += shell_volume;
        }
        (s.time_s, if volume > 0.0 { sum / volume } else { s.temperatures[interface] })
    }).collect()
}

/// Get average shell temperature over time
pub fn shell_average_temperature(snapshots: &[HeatSnapshot]) -> Vec<(f64, f64)> {
    snapshots.iter().map(|s| {
        let inner = s.inner_interface_idx;
        let outer = s.outer_interface_idx;
        // Volume-weighted average
        let mut sum = 0.0;
        let mut volume = 0.0;
        for i in inner..=outer {
            let r = s.radii[i];
            let dr = if i == 0 {
                s.radii[1] - s.radii[0]
            } else {
                s.radii[i] - s.radii[i - 1]
            };
            // Use spherical shell volume: V = 4π/3 * (r_outer³ - r_inner³)
            let r_inner = (r - dr / 2.0).max(0.0);
            let r_outer = r + dr / 2.0;
            let shell_volume = 4.0 / 3.0 * std::f64::consts::PI * (r_outer.powi(3) - r_inner.powi(3));
            sum += s.temperatures[i] * shell_volume;
            volume += shell_volume;
        }
        (s.time_s, if volume > 0.0 { sum / volume } else { s.temperatures[inner] })
    }).collect()
}

// Legacy aliases for backward compatibility
pub fn ball_center_temperature(snapshots: &[HeatSnapshot]) -> Vec<(f64, f64)> {
    liquid_center_temperature(snapshots)
}

pub fn ball_average_temperature(snapshots: &[HeatSnapshot]) -> Vec<(f64, f64)> {
    liquid_average_temperature(snapshots)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simulation_runs() {
        let params = SphereHeatParams::default();
        let snapshots = simulate_sphere_heat_numerical(&params, 60.0, 5);
        assert_eq!(snapshots.len(), 5);
    }

    #[test]
    fn test_temperature_decreases() {
        let params = SphereHeatParams::default();
        let snapshots = simulate_sphere_heat_numerical(&params, 300.0, 10);

        // Ball center should cool down
        let initial_center = snapshots[0].temperatures[0];
        let final_center = snapshots.last().unwrap().temperatures[0];
        assert!(final_center < initial_center, "Ball should cool down");
    }

    #[test]
    fn test_approaches_equilibrium() {
        let params = SphereHeatParams::default();
        let snapshots = simulate_sphere_heat_numerical(&params, 100000.0, 20); // ~28 hours

        // Should approach room temperature (Robin BC with larger domain is slower to equilibrate)
        let final_center = snapshots.last().unwrap().temperatures[0];
        assert!(
            (final_center - params.room_temp_c).abs() < 50.0,
            "Should approach room temp, got {}",
            final_center
        );
    }

    #[test]
    fn test_ceramic_vs_steel_shell() {
        // With radiation included, ceramic (ε=0.90) radiates much more heat
        // than steel (ε=0.10). This dominates over steel's higher thermal
        // conductivity (k=50 vs k=1), making ceramic cool the liquid faster.
        let base = SphereHeatParams {
            liquid_temp_c: 90.0,
            shell_temp_c: 90.0,
            room_temp_c: 20.0,
            ..SphereHeatParams::default()
        };
        let ceramic = SphereHeatParams {
            shell_material: HeatMaterial::ceramic(),
            ..base.clone()
        };
        let steel = SphereHeatParams {
            shell_material: HeatMaterial::steel(),
            ..base
        };

        let snaps_ceramic = simulate_sphere_heat_numerical(&ceramic, 600.0, 20);
        let snaps_steel = simulate_sphere_heat_numerical(&steel, 600.0, 20);

        let center_ceramic = liquid_center_temperature(&snaps_ceramic);
        let center_steel = liquid_center_temperature(&snaps_steel);

        // At 10 minutes, ceramic's high emissivity drives faster cooling
        let final_ceramic = center_ceramic.last().unwrap().1;
        let final_steel = center_steel.last().unwrap().1;
        assert!(
            final_ceramic < final_steel,
            "Ceramic shell (ε=0.90) should cool faster than steel (ε=0.10): ceramic={:.1} steel={:.1}",
            final_ceramic, final_steel
        );
    }

    #[test]
    fn test_temperature_continuity_at_interface() {
        // Temperature should be approximately continuous across the
        // liquid-shell interface (no sharp jumps) after some time.
        let params = SphereHeatParams {
            liquid_temp_c: 90.0,
            shell_temp_c: 22.0,  // Cold shell, will show interface effects
            room_temp_c: 22.0,
            shell_material: HeatMaterial::ceramic(),
            ..SphereHeatParams::default()
        };
        let snapshots = simulate_sphere_heat_numerical(&params, 60.0, 10);

        // After 60 seconds, check the interface region
        let final_snap = snapshots.last().unwrap();
        let idx = final_snap.inner_interface_idx;

        // Temperature jump at interface should be small (< 10°C) after 60s
        let t_liquid_side = final_snap.temperatures[idx.saturating_sub(1)];
        let t_shell_side = final_snap.temperatures[idx + 1];
        let jump = (t_liquid_side - t_shell_side).abs();
        assert!(
            jump < 10.0,
            "Interface temperature jump should be small after 60s, got {:.1}°C ({:.1} vs {:.1})",
            jump, t_liquid_side, t_shell_side
        );
    }

    #[test]
    fn test_temperatures_stay_bounded() {
        // No temperature should exceed initial max or drop below room temp.
        let params = SphereHeatParams {
            liquid_temp_c: 95.0,
            shell_temp_c: 80.0,
            room_temp_c: 20.0,
            shell_material: HeatMaterial::ceramic(),
            ..SphereHeatParams::default()
        };
        let snapshots = simulate_sphere_heat_numerical(&params, 300.0, 20);

        let t_max_init = 95.0_f64;
        let t_min_init = 20.0_f64;

        for snap in &snapshots {
            for (i, &t) in snap.temperatures.iter().enumerate() {
                assert!(
                    t <= t_max_init + 1.0,
                    "Temperature at r={:.4} t={:.0}s is {:.1}°C, exceeds initial max {:.1}°C",
                    snap.radii[i], snap.time_s, t, t_max_init
                );
                assert!(
                    t >= t_min_init - 1.0,
                    "Temperature at r={:.4} t={:.0}s is {:.1}°C, below room temp {:.1}°C",
                    snap.radii[i], snap.time_s, t, t_min_init
                );
            }
        }
    }

    #[test]
    fn test_symmetry_bc_at_center() {
        // The center (r=0) should have zero radial gradient (symmetry BC).
        // After some evolution, T[0] and T[1] should be very close.
        let params = SphereHeatParams {
            liquid_temp_c: 90.0,
            shell_temp_c: 90.0,
            room_temp_c: 20.0,
            ..SphereHeatParams::default()
        };
        let snapshots = simulate_sphere_heat_numerical(&params, 300.0, 10);

        for snap in &snapshots[1..] {
            let gradient = (snap.temperatures[1] - snap.temperatures[0])
                / (snap.radii[1] - snap.radii[0]);
            // Gradient at center should be small relative to overall temp range
            let temp_range = snap.temperatures[0] - params.room_temp_c;
            if temp_range.abs() > 1.0 {
                let normalized_gradient = gradient.abs() * snap.radii[1] / temp_range.abs();
                assert!(
                    normalized_gradient < 0.5,
                    "Center gradient too large at t={:.0}s: normalized={:.3}",
                    snap.time_s, normalized_gradient
                );
            }
        }
    }

    #[test]
    fn test_preheated_vs_cold_shell() {
        // Coffee in a preheated shell should stay warmer than in a cold shell
        let preheated = SphereHeatParams {
            liquid_temp_c: 90.0,
            shell_temp_c: 90.0,
            room_temp_c: 20.0,
            shell_material: HeatMaterial::ceramic(),
            ..SphereHeatParams::default()
        };
        let cold = SphereHeatParams {
            shell_temp_c: 20.0,
            ..preheated.clone()
        };

        let snaps_hot = simulate_sphere_heat_numerical(&preheated, 600.0, 20);
        let snaps_cold = simulate_sphere_heat_numerical(&cold, 600.0, 20);

        let avg_hot = liquid_average_temperature(&snaps_hot);
        let avg_cold = liquid_average_temperature(&snaps_cold);

        // At every time step, preheated should be >= cold
        for (h, c) in avg_hot.iter().zip(avg_cold.iter()) {
            assert!(
                h.1 >= c.1 - 0.5,  // small tolerance for numerical noise
                "Preheated should be warmer: t={:.0}s hot={:.1} cold={:.1}",
                h.0, h.1, c.1
            );
        }
    }

    #[test]
    fn test_evaporative_cooling_effect() {
        // With evaporation enabled, liquid should cool faster
        let base = SphereHeatParams {
            liquid_temp_c: 90.0,
            shell_temp_c: 90.0,
            room_temp_c: 20.0,
            shell_material: HeatMaterial::ceramic(),
            evaporation_factor: 0.0,
            ..SphereHeatParams::default()
        };
        let with_evap = SphereHeatParams {
            evaporation_factor: 1.0,
            ..base.clone()
        };

        let snaps_no_evap = simulate_sphere_heat_numerical(&base, 300.0, 10);
        let snaps_evap = simulate_sphere_heat_numerical(&with_evap, 300.0, 10);

        let avg_no_evap = liquid_average_temperature(&snaps_no_evap);
        let avg_evap = liquid_average_temperature(&snaps_evap);

        let final_no_evap = avg_no_evap.last().unwrap().1;
        let final_evap = avg_evap.last().unwrap().1;
        assert!(
            final_evap < final_no_evap,
            "Evaporative cooling should make liquid cooler: evap={:.1} no_evap={:.1}",
            final_evap, final_no_evap
        );
    }

    #[test]
    fn test_convection_coefficient_effect() {
        // Higher convection coefficient should lead to faster cooling
        let low_h = SphereHeatParams {
            liquid_temp_c: 90.0,
            shell_temp_c: 90.0,
            room_temp_c: 20.0,
            convection_coeff: 5.0,
            shell_material: HeatMaterial::ceramic(),
            ..SphereHeatParams::default()
        };
        let high_h = SphereHeatParams {
            convection_coeff: 25.0,
            ..low_h.clone()
        };

        let snaps_low = simulate_sphere_heat_numerical(&low_h, 600.0, 10);
        let snaps_high = simulate_sphere_heat_numerical(&high_h, 600.0, 10);

        let avg_low = liquid_average_temperature(&snaps_low);
        let avg_high = liquid_average_temperature(&snaps_high);

        let final_low = avg_low.last().unwrap().1;
        let final_high = avg_high.last().unwrap().1;
        assert!(
            final_high < final_low,
            "Higher convection should cool faster: h=25 → {:.1}°C, h=5 → {:.1}°C",
            final_high, final_low
        );
    }

    /// Benchmark: all 4 panel simulations (hot/cold × enclosed/exposed)
    /// must complete in < 2 seconds to avoid WASM UI hangs.
    #[test]
    fn test_panel_timing() {
        let inner_r_cm = (3.0 * 180.0_f64 / (4.0 * std::f64::consts::PI)).cbrt();
        let outer_r_cm = inner_r_cm + 5.0 / 10.0;
        let chart_max_r_cm = 10.0_f64;
        let dur = 60.0 * 60.0;
        let snaps = 100;

        let base = SphereHeatParams {
            outer_radius_m: outer_r_cm / 100.0,
            inner_radius_m: inner_r_cm / 100.0,
            shell_temp_c: 90.0,
            liquid_temp_c: 90.0,
            room_temp_c: 22.0,
            shell_material: HeatMaterial::ceramic(),
            liquid_material: HeatMaterial::water(),
            medium_material: HeatMaterial::air(),
            evaporation_factor: 0.0,
            convection_coeff: 10.0,
            r_max_multiplier: chart_max_r_cm / outer_r_cm,
        };

        let t0 = std::time::Instant::now();
        simulate_sphere_heat_numerical(&base, dur, snaps);
        simulate_sphere_heat_numerical(&SphereHeatParams { shell_temp_c: 22.0, ..base.clone() }, dur, snaps);
        simulate_sphere_heat_numerical(&SphereHeatParams {
            shell_material: HeatMaterial::air(),
            evaporation_factor: 1.0,
            convection_coeff: 15.0,
            ..base.clone()
        }, dur, snaps);
        simulate_sphere_heat_numerical(&SphereHeatParams {
            shell_material: HeatMaterial::air(),
            evaporation_factor: 1.0,
            convection_coeff: 15.0,
            ..base
        }, dur, snaps);
        let total = t0.elapsed();

        assert!(total.as_secs() < 2, "4 sims took {:?}, must be < 2s", total);
    }
}

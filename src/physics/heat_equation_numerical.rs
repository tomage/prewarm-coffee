//! Simulation of the heat equation for a sphere in a finite air domain
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
//! - Finite-volume interfaces: each grid segment carries one material's
//!   conductivity, and interface nodes own capacity from both sides
//! - Open-top losses (convection, radiation, evaporation) at the liquid
//!   surface, scaled by lid openness
//! - Convection and radiation at the shell's outer surface; the air at r_max
//!   is held at room temperature (`OuterBoundary::Ambient`), or the domain
//!   ends insulated for closed-system tests (`OuterBoundary::Insulated`)

use super::constants::{
    CELSIUS_TO_KELVIN, CP_AIR, K_AIR, RHO_AIR, STEFAN_BOLTZMANN, WATER_DENSITY,
    WATER_EFFECTIVE_CONDUCTIVITY, WATER_EMISSIVITY, WATER_SPECIFIC_HEAT,
};
use super::materials::MugMaterial;
use super::surface::{open_top_loss, sphere_open_fraction};

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

impl HeatMaterial {
    /// Thermal diffusivity α = k/(ρc) in m²/s
    pub fn diffusivity(&self) -> f64 {
        self.conductivity / (self.density * self.specific_heat)
    }

    /// Shell properties for a mug material, with the surface emissivity
    /// supplied separately (the user can override it).
    pub fn from_mug(material: MugMaterial, emissivity: f64) -> Self {
        Self {
            conductivity: material.thermal_conductivity(),
            density: material.density(),
            specific_heat: material.specific_heat(),
            emissivity,
        }
    }

    /// Stainless steel shell
    pub fn steel() -> Self {
        let material = MugMaterial::StainlessSteel;
        Self::from_mug(material, material.emissivity())
    }

    /// Ceramic shell
    pub fn ceramic() -> Self {
        let material = MugMaterial::Ceramic;
        Self::from_mug(material, material.emissivity())
    }

    /// Air properties (at ~25°C)
    pub fn air() -> Self {
        Self {
            conductivity: K_AIR,
            density: RHO_AIR,
            specific_heat: CP_AIR,
            emissivity: 1.0, // not physically meaningful for air medium
        }
    }

    /// Water properties (effective conductivity includes natural convection;
    /// see `WATER_EFFECTIVE_CONDUCTIVITY`)
    pub fn water() -> Self {
        Self {
            conductivity: WATER_EFFECTIVE_CONDUCTIVITY,
            density: WATER_DENSITY,
            specific_heat: WATER_SPECIFIC_HEAT,
            emissivity: WATER_EMISSIVITY,
        }
    }

    /// Coffee (similar to water)
    pub fn coffee() -> Self {
        Self::water()
    }
}

/// Boundary condition at the last grid point, r_max.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum OuterBoundary {
    /// Fixed at room temperature (Dirichlet): the far field of the air region.
    #[default]
    Ambient,
    /// Zero heat flux. With `r_max_multiplier = 1.0` the domain ends at the
    /// shell's outer surface, so coffee and shell form a closed system.
    Insulated,
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
    /// Lid openness (0.0–1.0): 0 = sealed lid, 1 = open mug with no lid.
    /// The open share of the liquid surface, `surface::sphere_open_fraction`,
    /// loses heat by convection, radiation and evaporation.
    pub exposure: f64,
    /// Convective heat transfer coefficient of the room air in W/(m²·K), at
    /// the shell's outer surface and the open liquid surface.
    /// Typical values: 5-25 for natural convection in air, 10-200 for forced convection.
    /// 0.0 leaves radiation as the only loss to the room.
    pub convection_coeff: f64,
    /// r_max = outer_radius_m * r_max_multiplier. Default 3.0.
    pub r_max_multiplier: f64,
    /// Boundary condition at r_max.
    pub outer_boundary: OuterBoundary,
}

impl Default for SphereHeatParams {
    fn default() -> Self {
        Self {
            outer_radius_m: 0.05, // 5 cm outer radius
            inner_radius_m: 0.04, // 4 cm inner radius (1 cm thick shell)
            shell_temp_c: 100.0,  // 100°C shell
            liquid_temp_c: 200.0, // 200°C liquid (hotter than shell)
            room_temp_c: 20.0,    // 20°C room
            shell_material: HeatMaterial::steel(),
            liquid_material: HeatMaterial::water(),
            medium_material: HeatMaterial::air(),
            exposure: 0.0,          // Sealed lid by default
            convection_coeff: 10.0, // Natural convection in air ~10 W/(m²·K)
            r_max_multiplier: 3.0,
            outer_boundary: OuterBoundary::Ambient,
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

/// Number of backward-Euler steps before switching to Crank–Nicolson.
///
/// Crank–Nicolson's amplification factor tends to −1 for stiff modes, so the
/// initial temperature jump at a high-diffusivity shell would flip sign every
/// step instead of decaying. Backward Euler damps those modes strongly; a few
/// such steps at the start (Rannacher start-up) remove them, and the run then
/// keeps Crank–Nicolson's second-order accuracy.
const RANNACHER_STEPS: usize = 4;

/// Radial bounds of node `i`'s control volume: the half-points
/// [r_{i−½}, r_{i+½}], with r_{−½} = 0 and the last node's volume ending at
/// its own radius.
fn control_volume_bounds(radii: &[f64], i: usize) -> (f64, f64) {
    let lo = if i == 0 {
        0.0
    } else {
        0.5 * (radii[i - 1] + radii[i])
    };
    let hi = if i + 1 < radii.len() {
        0.5 * (radii[i] + radii[i + 1])
    } else {
        radii[i]
    };
    (lo, hi)
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
/// Uses a finite-volume discretisation with Crank–Nicolson time stepping
/// after a short backward-Euler start-up (`RANNACHER_STEPS`).
/// Three regions: liquid (0 to r_inner), shell (r_inner to r_outer), air (r_outer to r_max)
/// Returns snapshots of temperature distribution at specified times.
pub fn simulate_sphere_heat_numerical(
    params: &SphereHeatParams,
    duration_s: f64,
    num_snapshots: usize,
) -> Vec<HeatSnapshot> {
    // Spatial discretization with adaptive grid
    // Uses piecewise-uniform grids: coarse in bulk, fine near interfaces
    // Scale resolution inversely with duration to keep total work bounded
    // Domain size from multiplier; the boundary condition at r_max is `outer_boundary`.
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
    let inner_interface_idx = radii
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| {
            (*a - r_inner)
                .abs()
                .partial_cmp(&(*b - r_inner).abs())
                .unwrap()
        })
        .unwrap()
        .0;
    let outer_interface_idx = radii
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| {
            (*a - r_outer)
                .abs()
                .partial_cmp(&(*b - r_outer).abs())
                .unwrap()
        })
        .unwrap()
        .0;

    // Finite-volume discretisation of ρc ∂T/∂t = (1/r²) ∂/∂r(r² k ∂T/∂r).
    // Grid segment j = [r_j, r_{j+1}] lies in one region; the interface
    // nodes sit on its boundaries. Node i owns the spherical shell between
    // the half-points r_{i−½} and r_{i+½} (r_{−½} = 0, and the last node's
    // shell ends at r_max), so an interface node's control volume is part
    // liquid and part shell (or part shell and part medium).
    let segment_value = |j: usize, liquid: f64, shell: f64, medium: f64| {
        if j < inner_interface_idx {
            liquid
        } else if j < outer_interface_idx {
            shell
        } else {
            medium
        }
    };
    let segment_rho_c = |j: usize| {
        segment_value(
            j,
            params.liquid_material.density * params.liquid_material.specific_heat,
            params.shell_material.density * params.shell_material.specific_heat,
            params.medium_material.density * params.medium_material.specific_heat,
        )
    };
    let segment_initial_temp = |j: usize| {
        segment_value(
            j,
            params.liquid_temp_c,
            params.shell_temp_c,
            params.room_temp_c,
        )
    };

    // Half-point radii
    let r_half: Vec<f64> = (0..n_points - 1)
        .map(|i| 0.5 * (radii[i] + radii[i + 1]))
        .collect();

    // Heat capacity of each node's control volume, per steradian [J/(K·sr)],
    // split into the parts left and right of the node. Using ρc_i·r_i²·dr
    // instead gives an interface node one side's ρc over the whole volume,
    // which creates or destroys heat when the two sides exchange it.
    let capacity_parts: Vec<(f64, f64)> = (0..n_points)
        .map(|i| {
            let (lo, hi) = control_volume_bounds(&radii, i);
            let r3 = radii[i].powi(3);
            let left = if i > 0 {
                segment_rho_c(i - 1) * (r3 - lo.powi(3)) / 3.0
            } else {
                0.0
            };
            let right = if i < n_points - 1 {
                segment_rho_c(i) * (hi.powi(3) - r3) / 3.0
            } else {
                0.0
            };
            (left, right)
        })
        .collect();
    let capacity: Vec<f64> = capacity_parts.iter().map(|(l, r)| l + r).collect();

    // Initial temperature: the control-volume average of the piecewise
    // initial field, so the discrete heat content matches the physical one.
    let mut temp: Vec<f64> = (0..n_points)
        .map(|i| {
            let (left, right) = capacity_parts[i];
            // The clamped segment index is only reached where that part is 0 (i = 0, i = n−1).
            let t_left = segment_initial_temp(i.saturating_sub(1));
            let t_right = segment_initial_temp(i.min(n_points - 2));
            (left * t_left + right * t_right) / (left + right)
        })
        .collect();

    // Pre-compute spacing between adjacent points
    let dr_arr: Vec<f64> = (0..n_points - 1).map(|i| radii[i + 1] - radii[i]).collect();

    // Conductance of each face r_{j+½}, per steradian [W/(K·sr)]. Segment j
    // lies wholly in one material, so its conductivity is that material's;
    // flux continuity across an interface holds because the interface node's
    // two faces carry the two materials' fluxes into one shared control volume.
    let face_conductance: Vec<f64> = (0..n_points - 1)
        .map(|j| {
            let k = segment_value(
                j,
                params.liquid_material.conductivity,
                params.shell_material.conductivity,
                params.medium_material.conductivity,
            );
            r_half[j].powi(2) * k / dr_arr[j]
        })
        .collect();

    // Tridiagonal diffusion operator
    // L(T)_i = lower[i]*T_{i-1} + diag[i]*T_i + upper[i]*T_{i+1} + source[i]
    // Each face conductance r_{j+½}²·k_j/dr_j is divided by the capacity of
    // the node it feeds, so Σ capacity·T is conserved by the face fluxes.
    // At r = 0 the face has zero area and this reduces to 6k/(ρc·dr²).
    let mut lower = vec![0.0; n_points];
    let mut upper = vec![0.0; n_points];
    let mut diag_op = vec![0.0; n_points];
    let mut source = vec![0.0; n_points];

    for i in 0..n_points {
        if i > 0 {
            lower[i] = face_conductance[i - 1] / capacity[i];
        }
        if i < n_points - 1 {
            upper[i] = face_conductance[i] / capacity[i];
        }
        diag_op[i] = -(lower[i] + upper[i]);
    }

    // Surface losses at an interface node: a flux h·(T − T_ref) [W/m²]
    // through area r_i² (per steradian) into capacity C_i gives the rate
    // coefficient h·r_i²/C_i [1/s] on the diagonal and h·r_i²·T_ref/C_i in
    // the source.
    let area_per_capacity = |i: usize| radii[i].powi(2) / capacity[i];
    let area_per_capacity_outer = area_per_capacity(outer_interface_idx);
    let area_per_capacity_inner = area_per_capacity(inner_interface_idx);
    let open_fraction = sphere_open_fraction(params.exposure);

    // Convection at outer interface: implicit treatment (linear in T)
    if params.convection_coeff > 0.0 {
        let conv_coeff = params.convection_coeff * area_per_capacity_outer;
        diag_op[outer_interface_idx] -= conv_coeff;
        source[outer_interface_idx] += conv_coeff * params.room_temp_c;
    }

    // Store base operator values at outer interface (without radiation)
    let base_diag_outer = diag_op[outer_interface_idx];
    let base_source_outer = source[outer_interface_idx];

    // Base operator values at the inner interface (without open-top losses)
    let base_diag_inner = diag_op[inner_interface_idx];
    let base_source_inner = source[inner_interface_idx];

    // Time stepping: the θ-scheme is unconditionally stable, so dt is chosen for accuracy:
    // duration/500, capped at 0.5 s and at the snapshot interval. Runs longer than
    // about 4 minutes therefore take 2 steps per simulated second (3,600 for 30 min).
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

    let mut step = 0;
    while snapshots.len() < num_snapshots {
        // θ-scheme weight: backward Euler (θ = 1) for the start-up steps,
        // Crank–Nicolson (θ = ½) after.
        let theta = if step < RANNACHER_STEPS { 1.0 } else { 0.5 };
        let implicit_dt = theta * dt;
        let explicit_dt = (1.0 - theta) * dt;
        step += 1;

        // Linearized radiation at mug outer surface: h_rad = 4·ε·σ·T_m³
        // Re-linearized each step for accuracy, treated implicitly like convection.
        {
            let t_k = temp[outer_interface_idx] + CELSIUS_TO_KELVIN;
            let t_room_k = params.room_temp_c + CELSIUS_TO_KELVIN;
            let t_m = 0.5 * (t_k + t_room_k);
            let h_rad = 4.0 * params.shell_material.emissivity * STEFAN_BOLTZMANN * t_m.powi(3);
            let rad_coeff = h_rad * area_per_capacity_outer;
            diag_op[outer_interface_idx] = base_diag_outer - rad_coeff;
            source[outer_interface_idx] = base_source_outer + rad_coeff * params.room_temp_c;
        }

        // Coffee surface heat loss through the open top, over the share of the
        // liquid sphere's surface that the lid opening maps to: convection,
        // radiation and evaporation, linearised at the current surface
        // temperature as h·(T − t_ref) and treated implicitly. t_ref sits below
        // room while the coffee evaporates.
        if open_fraction > 0.0 {
            let loss = open_top_loss(
                temp[inner_interface_idx],
                params.room_temp_c,
                params.convection_coeff,
            );
            let coeff = open_fraction * loss.h * area_per_capacity_inner;
            diag_op[inner_interface_idx] = base_diag_inner - coeff;
            source[inner_interface_idx] = base_source_inner + coeff * loss.t_ref;
        }

        // θ-scheme: (I − θ·dt·L) T^{n+1} = (I + (1−θ)·dt·L) T^n + dt·source
        for i in 0..n_points {
            let mut r = (1.0 + explicit_dt * diag_op[i]) * temp[i] + dt * source[i];
            if i > 0 {
                r += explicit_dt * lower[i] * temp[i - 1];
            }
            if i < n_points - 1 {
                r += explicit_dt * upper[i] * temp[i + 1];
            }
            rhs[i] = r;
            sub_lhs[i] = -implicit_dt * lower[i];
            dia_lhs[i] = 1.0 - implicit_dt * diag_op[i];
            sup_lhs[i] = -implicit_dt * upper[i];
        }

        // Ambient: r_max is held at room temperature (identity row).
        // Insulated: the last row is the zero-flux finite-volume row built above.
        if params.outer_boundary == OuterBoundary::Ambient {
            let last = n_points - 1;
            rhs[last] = params.room_temp_c;
            sub_lhs[last] = 0.0;
            dia_lhs[last] = 1.0;
            sup_lhs[last] = 0.0;
        }

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
    snapshots
        .iter()
        .map(|s| (s.time_s, s.temperatures[0]))
        .collect()
}

/// Volume-weighted mean temperature over [r_lo, r_hi], using the solver's
/// control volumes clipped to that range: node i covers [r_{i−½}, r_{i+½}],
/// so an interface node contributes only its share inside the range. ρc is
/// uniform within a region, so these are the solver's own heat-capacity
/// weights. `None` when the range holds no volume.
fn volume_average(s: &HeatSnapshot, r_lo: f64, r_hi: f64) -> Option<f64> {
    let mut sum = 0.0;
    let mut volume = 0.0;
    for i in 0..s.radii.len() {
        let (cv_lo, cv_hi) = control_volume_bounds(&s.radii, i);
        let lo = cv_lo.max(r_lo);
        let hi = cv_hi.min(r_hi);
        if hi > lo {
            let v = hi.powi(3) - lo.powi(3);
            sum += s.temperatures[i] * v;
            volume += v;
        }
    }
    (volume > 0.0).then(|| sum / volume)
}

/// Get average liquid temperature over time
pub fn liquid_average_temperature(snapshots: &[HeatSnapshot]) -> Vec<(f64, f64)> {
    snapshots
        .iter()
        .map(|s| {
            let r_inner = s.radii[s.inner_interface_idx];
            let t = volume_average(s, 0.0, r_inner).unwrap_or(s.temperatures[0]);
            (s.time_s, t)
        })
        .collect()
}

/// Get volume-weighted liquid near-wall average temperature over time.
/// `wall_fraction` is the fraction of liquid radius from the wall inward
/// (e.g. 0.10 = outer 10% of the liquid radius).
pub fn liquid_near_wall_temperature(
    snapshots: &[HeatSnapshot],
    wall_fraction: f64,
) -> Vec<(f64, f64)> {
    snapshots
        .iter()
        .map(|s| {
            let interface = s.inner_interface_idx;
            let r_inner = s.radii[interface];
            let r_start = r_inner * (1.0 - wall_fraction);
            let t = volume_average(s, r_start, r_inner).unwrap_or(s.temperatures[interface]);
            (s.time_s, t)
        })
        .collect()
}

/// Get average shell temperature over time
pub fn shell_average_temperature(snapshots: &[HeatSnapshot]) -> Vec<(f64, f64)> {
    snapshots
        .iter()
        .map(|s| {
            let r_inner = s.radii[s.inner_interface_idx];
            let r_outer = s.radii[s.outer_interface_idx];
            let t = volume_average(s, r_inner, r_outer)
                .unwrap_or(s.temperatures[s.inner_interface_idx]);
            (s.time_s, t)
        })
        .collect()
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

        // Should approach room temperature (a larger air domain is slower to equilibrate)
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
            shell_temp_c: 22.0, // Cold shell, will show interface effects
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
            jump,
            t_liquid_side,
            t_shell_side
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
                    snap.radii[i],
                    snap.time_s,
                    t,
                    t_max_init
                );
                assert!(
                    t >= t_min_init - 1.0,
                    "Temperature at r={:.4} t={:.0}s is {:.1}°C, below room temp {:.1}°C",
                    snap.radii[i],
                    snap.time_s,
                    t,
                    t_min_init
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
            let gradient =
                (snap.temperatures[1] - snap.temperatures[0]) / (snap.radii[1] - snap.radii[0]);
            // Gradient at center should be small relative to overall temp range
            let temp_range = snap.temperatures[0] - params.room_temp_c;
            if temp_range.abs() > 1.0 {
                let normalized_gradient = gradient.abs() * snap.radii[1] / temp_range.abs();
                assert!(
                    normalized_gradient < 0.5,
                    "Center gradient too large at t={:.0}s: normalized={:.3}",
                    snap.time_s,
                    normalized_gradient
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
                h.1 >= c.1 - 0.5, // small tolerance for numerical noise
                "Preheated should be warmer: t={:.0}s hot={:.1} cold={:.1}",
                h.0,
                h.1,
                c.1
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
            exposure: 0.0,
            ..SphereHeatParams::default()
        };
        let with_evap = SphereHeatParams {
            exposure: 1.0,
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
            final_evap,
            final_no_evap
        );
    }

    #[test]
    fn open_top_cools_liquid_at_room_below_room() {
        // Evaporation has no cutoff: liquid starting at room temperature keeps
        // evaporating and drops below room, without NaN or blow-up.
        let params = SphereHeatParams {
            liquid_temp_c: 22.0,
            shell_temp_c: 22.0,
            room_temp_c: 22.0,
            shell_material: HeatMaterial::ceramic(),
            exposure: 1.0,
            ..SphereHeatParams::default()
        };
        let snaps = simulate_sphere_heat_numerical(&params, 600.0, 10);
        assert!(snaps
            .iter()
            .all(|s| s.temperatures.iter().all(|t| t.is_finite())));
        let final_avg = liquid_average_temperature(&snaps).last().unwrap().1;
        assert!(
            final_avg < 22.0 && final_avg > 10.0,
            "liquid average after 10 min = {final_avg:.3} °C, expected in (10, 22)"
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
            final_high,
            final_low
        );
    }

    /// Number of times the step-to-step change of `series` flips sign,
    /// ignoring steps with no change.
    fn sign_changes(series: impl Iterator<Item = f64>) -> usize {
        let values: Vec<f64> = series.collect();
        let signs: Vec<f64> = values
            .windows(2)
            .map(|w| w[1] - w[0])
            .filter(|d| *d != 0.0)
            .map(f64::signum)
            .collect();
        signs.windows(2).filter(|w| w[0] != w[1]).count()
    }

    /// A steel shell's interface nodes see α·dt/dr² ≈ 10⁴, where pure
    /// Crank–Nicolson barely damps the initial jump and the node alternates
    /// between two temperatures every step. Recording every step exposes it.
    #[test]
    fn interface_nodes_do_not_ring() {
        let inner_radius_m = (3.0 * 350e-6 / (4.0 * std::f64::consts::PI)).cbrt();
        let params = SphereHeatParams {
            inner_radius_m,
            outer_radius_m: inner_radius_m + 5e-3,
            liquid_temp_c: 95.0,
            shell_temp_c: 22.0,
            room_temp_c: 22.0,
            shell_material: HeatMaterial::from_mug(MugMaterial::StainlessSteel, 0.1),
            ..SphereHeatParams::default()
        };
        let duration_s = 1800.0;
        // dt is capped at 0.5 s, so one snapshot per step.
        let num_steps = (duration_s / 0.5) as usize;
        let snapshots = simulate_sphere_heat_numerical(&params, duration_s, num_steps + 1);
        assert_eq!(snapshots.len(), num_steps + 1);
        assert!(
            (snapshots[1].time_s - 0.5).abs() < 1e-9,
            "one snapshot per step"
        );

        let first = &snapshots[0];
        for (name, idx) in [
            ("inner", first.inner_interface_idx),
            ("outer", first.outer_interface_idx),
        ] {
            let series = snapshots[RANNACHER_STEPS..]
                .iter()
                .map(|s| s.temperatures[idx]);
            let changes = sign_changes(series);
            assert!(
                changes <= 2,
                "{name} interface node changes direction {changes} times"
            );
        }
    }

    /// Benchmark: all 4 panel simulations (hot/cold × enclosed/exposed)
    /// must complete in < 2 seconds to avoid WASM UI hangs. Wall-clock, so it
    /// runs with the other perf checks (`make bench`), not under coverage
    /// instrumentation alongside parallel tests.
    #[test]
    #[ignore = "timing-sensitive; run with make bench"]
    fn perf_panel_timing() {
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
            exposure: 0.0,
            convection_coeff: 10.0,
            r_max_multiplier: chart_max_r_cm / outer_r_cm,
            outer_boundary: OuterBoundary::Ambient,
        };

        let t0 = std::time::Instant::now();
        simulate_sphere_heat_numerical(&base, dur, snaps);
        simulate_sphere_heat_numerical(
            &SphereHeatParams {
                shell_temp_c: 22.0,
                ..base.clone()
            },
            dur,
            snaps,
        );
        simulate_sphere_heat_numerical(
            &SphereHeatParams {
                shell_material: HeatMaterial::air(),
                exposure: 1.0,
                convection_coeff: 15.0,
                ..base.clone()
            },
            dur,
            snaps,
        );
        simulate_sphere_heat_numerical(
            &SphereHeatParams {
                shell_material: HeatMaterial::air(),
                exposure: 1.0,
                convection_coeff: 15.0,
                ..base
            },
            dur,
            snaps,
        );
        let total = t0.elapsed();

        assert!(total.as_secs() < 2, "4 sims took {:?}, must be < 2s", total);
    }
}

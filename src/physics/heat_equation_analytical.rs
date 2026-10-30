//! Analytical solution for the three-region spherical heat equation
//!
//! Solves the heat equation for a sphere with three concentric regions:
//! - Region 1 (liquid): 0 < r < r_i
//! - Region 2 (shell): r_i < r < r_o
//! - Region 3 (air): r_o < r < r_max
//!
//! The eigenvalues are found numerically by solving the transcendental equation
//! arising from interface matching conditions.
//!
//! ## Finite-domain limitation
//!
//! This solver uses a finite domain [0, r_max] with the air held at ambient
//! temperature at r_max (Dirichlet). Convection and linearised radiation act at
//! r_o as the surface loss `h_conv`. Physically the air extends to infinity, so
//! the solution depends on r_max: a small r_max pins room-temperature air close
//! to the mug and speeds cooling. In practice, r_max >= 10x r_outer is
//! sufficient for convergence (the boundary barely affects the liquid/shell
//! temperatures).
//!
//! ## Future: infinite-domain continuous-spectrum approach
//!
//! For r → ∞, both sin(β₃r)/r and cos(β₃r)/r decay as 1/r, so the far-field
//! BC is automatically satisfied for *any* λ > 0. The interface conditions at
//! r_i and r_o are exactly determined (not over-determined), meaning the
//! eigenvalue spectrum is **continuous** rather than discrete. The solution
//! becomes a Fourier-integral expansion:
//!
//!   T(r,t) = T_amb + ∫₀^∞ A(λ) φ(r,λ) exp(-λ² t) ρ(λ) dλ
//!
//! where ρ(λ) is the spectral density and A(λ) is computed from the initial
//! condition. This eliminates r_max entirely but requires:
//! - Numerically evaluating the λ-integral (fine grid over λ)
//! - Computing the spectral density ρ(λ) from the interface conditions
//! - A different normalization convention (eigenfunctions are not L² on [0,∞))

use super::constants::{CELSIUS_TO_KELVIN, K_AIR, STEFAN_BOLTZMANN};
use super::heat_equation_numerical::HeatMaterial;
use super::surface::three_region_h_inner;
use std::f64::consts::PI;

/// Parameters for the three-region analytical solution
#[derive(Debug, Clone, PartialEq)]
pub struct ThreeRegionParams {
    /// Inner radius (liquid-shell interface) in meters
    pub r_inner: f64,
    /// Outer radius (shell-air interface) in meters
    pub r_outer: f64,
    /// Maximum radius (far boundary) in meters
    pub r_max: f64,
    /// Thermal diffusivity of liquid (region 1) m²/s
    pub alpha_1: f64,
    /// Thermal diffusivity of shell (region 2) m²/s
    pub alpha_2: f64,
    /// Thermal diffusivity of air (region 3) m²/s
    pub alpha_3: f64,
    /// Thermal conductivity of liquid W/(m·K)
    pub k_1: f64,
    /// Thermal conductivity of shell W/(m·K)
    pub k_2: f64,
    /// Thermal conductivity of air W/(m·K)
    pub k_3: f64,
    /// Initial temperature of liquid °C
    pub t_liquid_init: f64,
    /// Initial temperature of shell °C
    pub t_shell_init: f64,
    /// Ambient temperature (at r_max) °C
    pub t_ambient: f64,
    /// Number of eigenvalues/modes to use
    pub num_modes: usize,
    /// Surface loss coefficient at the shell's outer surface r_o in W/(m²·K)
    /// (convection plus linearised radiation in production). Enters the flux
    /// balance −k₂ ∂T/∂r = −k₃ ∂T/∂r + h_conv·(T − T_amb) at r = r_o.
    pub h_conv: f64,
    /// Inner surface heat loss coefficient in W/(m²·K)
    ///
    /// Models heat lost through the lid opening at r = r_inner as a sink.
    /// Modifies the flux condition from pure continuity to:
    ///   k₁ ∂T/∂r = k₂ ∂T/∂r − h_inner·(T − T_amb)
    /// (the heat arriving from the liquid, −k₁∂T/∂r, feeds the shell and the sink).
    ///
    /// Production sets it from `surface::three_region_h_inner`: the open-top
    /// loss (convection, radiation, evaporation) linearised for one segment
    /// and scaled by the open share of the sphere. Set to 0.0 for a sealed lid (pure flux
    /// continuity at r_inner).
    pub h_inner: f64,
}

impl Default for ThreeRegionParams {
    fn default() -> Self {
        Self {
            r_inner: 0.022,   // 2.2 cm (liquid radius)
            r_outer: 0.025,   // 2.5 cm (outer shell radius)
            r_max: 0.075,     // 7.5 cm (3x outer radius)
            alpha_1: 1.43e-7, // Water
            alpha_2: 4.73e-7, // Ceramic
            alpha_3: 2.19e-5, // Air
            k_1: 0.6,         // Water conductivity
            k_2: 1.0,         // Ceramic conductivity
            k_3: K_AIR,       // Air conductivity
            t_liquid_init: 95.0,
            t_shell_init: 80.0,
            t_ambient: 20.0,
            num_modes: 50, // More modes for better convergence
            h_conv: 0.0,   // No additional convection (pure conduction)
            h_inner: 0.0,  // No inner surface heat loss
        }
    }
}

/// A mug as the simulator describes it, in the terms both heat-equation
/// solvers share: geometry, the shell material, temperatures, the outer
/// convection coefficient and the lid openness.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ThreeRegionMug {
    /// Liquid radius in meters
    pub r_inner: f64,
    /// Shell outer radius in meters
    pub r_outer: f64,
    /// Far boundary radius in meters
    pub r_max: f64,
    /// Shell material, including its surface emissivity
    pub shell: HeatMaterial,
    /// Initial liquid temperature °C
    pub liquid_temp_c: f64,
    /// Initial shell temperature °C
    pub shell_temp_c: f64,
    /// Room temperature °C
    pub room_temp_c: f64,
    /// Natural-convection coefficient in W/(m²·K)
    pub h_conv: f64,
    /// Lid openness 0–1
    pub exposure: f64,
    /// Number of eigenmodes
    pub num_modes: usize,
}

impl ThreeRegionParams {
    /// Parameters for `mug` with the whole run linearised once from its
    /// initial temperatures: the open-top loss at the mean of the initial
    /// liquid and room temperatures, the outer radiation at the initial
    /// liquid temperature. The first segment of a run solved this way is the
    /// reference the kernel tests compare against.
    pub fn for_mug(mug: &ThreeRegionMug) -> Self {
        Self::linearised(
            mug,
            0.5 * (mug.liquid_temp_c + mug.room_temp_c),
            mug.liquid_temp_c,
        )
    }

    /// Parameters for `mug`: water inside, air outside, and the surface
    /// losses linearised at fixed temperatures, because the eigenfunction
    /// expansion takes constant coefficients. The lid opening loses
    /// `surface::three_region_h_inner` at liquid temperature `liquid_lin_c`
    /// (r_inner); the outer surface loses `outer_surface_h` at outer-surface
    /// temperature `outer_lin_c` (r_outer).
    pub fn linearised(mug: &ThreeRegionMug, liquid_lin_c: f64, outer_lin_c: f64) -> Self {
        let water = HeatMaterial::water();
        let air = HeatMaterial::air();
        Self {
            r_inner: mug.r_inner,
            r_outer: mug.r_outer,
            r_max: mug.r_max,
            alpha_1: water.diffusivity(),
            alpha_2: mug.shell.diffusivity(),
            alpha_3: air.diffusivity(),
            k_1: water.conductivity,
            k_2: mug.shell.conductivity,
            k_3: air.conductivity,
            t_liquid_init: mug.liquid_temp_c,
            t_shell_init: mug.shell_temp_c,
            t_ambient: mug.room_temp_c,
            num_modes: mug.num_modes,
            h_conv: outer_surface_h(
                mug.h_conv,
                mug.shell.emissivity,
                outer_lin_c,
                mug.room_temp_c,
            ),
            h_inner: three_region_h_inner(mug.exposure, liquid_lin_c, mug.room_temp_c, mug.h_conv),
        }
    }
}

/// Outer-surface loss coefficient for the analytical solver, in W/(m²·K):
/// convection plus the shell's radiation as a secant toward room at surface
/// temperature `surface_c`, h_conv + 4εσT_m³ with T_m the mean of surface
/// and room temperatures.
pub fn outer_surface_h(h_conv: f64, emissivity: f64, surface_c: f64, room_temp_c: f64) -> f64 {
    let t_mean_k = 0.5 * (surface_c + room_temp_c) + CELSIUS_TO_KELVIN;
    h_conv + 4.0 * emissivity * STEFAN_BOLTZMANN * t_mean_k.powi(3)
}

/// An eigenmode of the three-region problem
#[derive(Debug, Clone, PartialEq)]
pub struct Eigenmode {
    /// The eigenvalue λ
    pub(crate) lambda: f64,
    /// Coefficient A (region 1)
    pub(crate) coeff_a: f64,
    /// Coefficient B (region 2, sin term)
    pub(crate) coeff_b: f64,
    /// Coefficient C (region 2, cos term)
    pub(crate) coeff_c: f64,
    /// Pole-free coefficient for region 3 eigenfunction:
    ///   φ₃(r) = D·(gd·sin(β₃r) + gn·cos(β₃r))/r
    /// where gd = cos(β₃R), gn = -sin(β₃R).
    /// This avoids the γ = -tan(β₃R) singularity when cos(β₃R) ≈ 0.
    pub(crate) coeff_d: f64,
    /// Dirichlet BC numerator: gn = -sin(β₃·r_max)
    pub(crate) gamma_num: f64,
    /// Dirichlet BC denominator: gd = cos(β₃·r_max)
    pub(crate) gamma_den: f64,
    /// Expansion coefficient for this mode (from initial conditions)
    pub(crate) amplitude: f64,
}

/// Result of three-region analytical solution
#[derive(Debug, Clone, PartialEq)]
pub struct ThreeRegionSnapshot {
    pub time_s: f64,
    pub radii: Vec<f64>,
    pub temperatures: Vec<f64>,
    pub r_inner: f64,
    pub r_outer: f64,
}

/// Compute pole-free Dirichlet BC parameters at r_max.
///
/// Returns (gn, gd) where gn = -sin(β₃R) and gd = cos(β₃R).
/// The region 3 eigenfunction is D·(gd·sin(β₃r) + gn·cos(β₃r))/r
/// which satisfies φ(r_max) = 0 since gd·sin(β₃R) + gn·cos(β₃R) = 0.
///
/// This avoids the numerical instability of γ = -tan(β₃R) which
/// diverges when cos(β₃R) ≈ 0.
fn dirichlet_gamma_polefree(lambda: f64, params: &ThreeRegionParams) -> (f64, f64) {
    let beta_3 = lambda / params.alpha_3.sqrt();
    let r_max = params.r_max;
    let gd = (beta_3 * r_max).cos();
    let gn = -(beta_3 * r_max).sin();
    (gn, gd)
}

/// Compute the determinant of the interface conditions matrix
/// Returns 0 when λ is an eigenvalue
pub(crate) fn interface_determinant(lambda: f64, params: &ThreeRegionParams) -> f64 {
    if lambda <= 0.0 {
        return f64::MAX;
    }

    let r_i = params.r_inner;
    let r_o = params.r_outer;

    // Wave numbers in each region: β = λ/√α
    let beta_1 = lambda / params.alpha_1.sqrt();
    let beta_2 = lambda / params.alpha_2.sqrt();
    let beta_3 = lambda / params.alpha_3.sqrt();

    let k_1 = params.k_1;
    let k_2 = params.k_2;
    let k_3 = params.k_3;

    // Dirichlet BC at r_max: γ = -sin(β₃R)/cos(β₃R) = -tan(β₃R)
    // To avoid poles (when cos = 0), we multiply column 4 of the
    // interface matrix by cos(β₃R). This gives a pole-free determinant.
    // γ = gn/gd where gn = -sin(β₃R), gd = cos(β₃R)
    let r_max = params.r_max;
    let h = params.h_conv;
    let sin_br = (beta_3 * r_max).sin();
    let cos_br = (beta_3 * r_max).cos();
    let gn = -sin_br; // gamma numerator
    let gd = cos_br; // gamma denominator

    // Precompute trig values
    let s1i = (beta_1 * r_i).sin();
    let c1i = (beta_1 * r_i).cos();
    let s2i = (beta_2 * r_i).sin();
    let c2i = (beta_2 * r_i).cos();
    let s2o = (beta_2 * r_o).sin();
    let c2o = (beta_2 * r_o).cos();

    let s3o = (beta_3 * r_o).sin();
    let c3o = (beta_3 * r_o).cos();

    // Build the 4x4 matrix of interface conditions (pole-free form)
    // Column 4 is multiplied by gamma_denom to clear the γ = gn/gd fraction.
    // φ₃(r) = D̃ (gd sin(β₃r) + gn cos(β₃r)) / r  (where D̃ = D/gd)

    // Row 1: temperature continuity at r_i
    let m11 = s1i;
    let m12 = -s2i;
    let m13 = -c2i;
    let m14 = 0.0;

    // Row 2: flux balance at r_i (with optional inner surface heat loss)
    //   k₁ ∂φ₁/∂r = k₂ ∂φ₂/∂r − h_inner·φ   at r = r_i  (sink)
    // Multiplied by r_i²: the h_inner term adds +h_inner·r_i·sin(β₁r_i) to col 1.
    let h_ri = params.h_inner * r_i;
    let m21 = k_1 * (r_i * beta_1 * c1i - s1i) + h_ri * s1i;
    let m22 = -k_2 * (r_i * beta_2 * c2i - s2i);
    let m23 = -k_2 * (-r_i * beta_2 * s2i - c2i);
    let m24 = 0.0;

    // Row 3: temperature continuity at r_o (pole-free column 4)
    // Original: -(s3o + γ c3o) → multiply by gd: -(gd s3o + gn c3o)
    let m31 = 0.0;
    let m32 = s2o;
    let m33 = c2o;
    let m34 = -(gd * s3o + gn * c3o);

    // Row 4: flux + convection at r_o (pole-free column 4)
    // Energy balance: -k₂ φ' = -k₃ φ' + h φ  →  k₂ φ' + h φ - k₃ φ' = 0
    // After ×r_o²: h term becomes h r_o (B s2o + C c2o), added to shell cols.
    let h_ro = h * r_o;
    let m41 = 0.0;
    let m42 = k_2 * (r_o * beta_2 * c2o - s2o) + h_ro * s2o;
    let m43 = k_2 * (-r_o * beta_2 * s2o - c2o) + h_ro * c2o;
    let m44 = -k_3 * (r_o * beta_3 * (gd * c3o - gn * s3o) - (gd * s3o + gn * c3o));

    // Compute 4x4 determinant using cofactor expansion
    // det = m11 * M11 - m12 * M12 + m13 * M13 - m14 * M14
    // where Mij is the 3x3 minor

    let minor_11 = m22 * (m33 * m44 - m34 * m43) - m23 * (m32 * m44 - m34 * m42)
        + m24 * (m32 * m43 - m33 * m42);
    let minor_12 = m21 * (m33 * m44 - m34 * m43) - m23 * (m31 * m44 - m34 * m41)
        + m24 * (m31 * m43 - m33 * m41);
    let minor_13 = m21 * (m32 * m44 - m34 * m42) - m22 * (m31 * m44 - m34 * m41)
        + m24 * (m31 * m42 - m32 * m41);
    let minor_14 = m21 * (m32 * m43 - m33 * m42) - m22 * (m31 * m43 - m33 * m41)
        + m23 * (m31 * m42 - m32 * m41);

    m11 * minor_11 - m12 * minor_12 + m13 * minor_13 - m14 * minor_14
}

/// Find eigenvalues by searching for zeros of the determinant
fn find_eigenvalues(params: &ThreeRegionParams, num_modes: usize) -> Vec<f64> {
    let mut eigenvalues = Vec::with_capacity(num_modes);

    // The eigenvalue λ is defined such that temporal decay is exp(-λ²t).
    // The spatial wave number in region i is β_i = λ/√α_i.
    // For a homogeneous sphere of radius R with diffusivity α:
    //   β = nπ/R, so λ = nπ√α/R
    //
    // For our three-region problem, the fundamental eigenvalue is approximately:
    //   λ₁ ≈ π × √α_liquid / r_inner
    // This is typically around 0.05 for water.

    let lambda_fundamental = PI * params.alpha_1.sqrt() / params.r_inner;

    // Start search from a small positive λ to catch low-frequency modes.
    // The bulk liquid mode can have λ << π√α/R (the single-region fundamental),
    // because the three-region interface conditions allow slower modes.
    // We skip λ = 0 exactly (degenerate trivial solution).
    let lambda_min = 0.001;
    let lambda_max = lambda_fundamental * (num_modes as f64 + 5.0) * 2.0;

    // Search density: 100 points per mode is enough for reliable sign-change detection
    let num_search_points = num_modes * 100;
    let d_lambda = (lambda_max - lambda_min) / num_search_points as f64;

    let mut prev_lambda = lambda_min;
    let mut prev_det = interface_determinant(prev_lambda, params);
    // The sample before `prev`, to spot a dip of |det| toward zero.
    let mut before: Option<(f64, f64)> = None;

    let push = |root: f64, eigenvalues: &mut Vec<f64>| {
        // Brackets never overlap, so the same root can only come back when it
        // sits on a shared sample point; anything farther apart than the
        // bisection tolerance is a distinct mode, however close.
        let is_new = eigenvalues
            .last()
            .is_none_or(|&last: &f64| (root - last).abs() > ROOT_TOLERANCE);
        if is_new {
            eigenvalues.push(root);
        }
    };

    for i in 1..=num_search_points {
        let lambda = lambda_min + i as f64 * d_lambda;
        let det = interface_determinant(lambda, params);

        if prev_det.is_finite() && det.is_finite() {
            if prev_det * det < 0.0 {
                // Sign change: one root (or an odd number) in the step.
                if let Some(root) = bisect_root(prev_lambda, lambda, params) {
                    push(root, &mut eigenvalues);
                }
            } else if let Some((lambda_b, det_b)) = before {
                // Two roots closer than a step leave the sign unchanged; the
                // determinant then dips toward zero between them. Find the
                // dip's extremum and, if it crosses zero, bisect either side.
                if det_b.is_finite()
                    && det_b * prev_det > 0.0
                    && prev_det.abs() < det_b.abs()
                    && prev_det.abs() <= det.abs()
                {
                    let (lambda_dip, det_dip) = dip_extremum(lambda_b, lambda, prev_det, params);
                    if det_dip * prev_det < 0.0 {
                        for (a, b) in [(lambda_b, lambda_dip), (lambda_dip, lambda)] {
                            if let Some(root) = bisect_root(a, b, params) {
                                push(root, &mut eigenvalues);
                            }
                        }
                    }
                }
            }
            if eigenvalues.len() >= num_modes {
                eigenvalues.truncate(num_modes);
                break;
            }
        }

        before = Some((prev_lambda, prev_det));
        prev_det = det;
        prev_lambda = lambda;
    }

    eigenvalues
}

/// The λ in [a, b] where the determinant comes closest to (or furthest
/// past) zero from the side of `sign_ref`, by golden-section search on
/// sign(sign_ref)·det, with the determinant there.
fn dip_extremum(mut a: f64, mut b: f64, sign_ref: f64, params: &ThreeRegionParams) -> (f64, f64) {
    let s = sign_ref.signum();
    let f = |l: f64| s * interface_determinant(l, params);
    let ratio = (5.0_f64.sqrt() - 1.0) / 2.0;
    let mut c = b - ratio * (b - a);
    let mut d = a + ratio * (b - a);
    let (mut fc, mut fd) = (f(c), f(d));
    for _ in 0..60 {
        if fc < 0.0 || fd < 0.0 || (b - a) < 1e-12 {
            break;
        }
        if fc < fd {
            b = d;
            d = c;
            fd = fc;
            c = b - ratio * (b - a);
            fc = f(c);
        } else {
            a = c;
            c = d;
            fc = fd;
            d = a + ratio * (b - a);
            fd = f(d);
        }
    }
    let lambda = if fc < fd { c } else { d };
    (lambda, interface_determinant(lambda, params))
}

/// Half-width to which `bisect_root` brackets an eigenvalue.
const ROOT_TOLERANCE: f64 = 1e-10;

/// Bisection method to find a root of the determinant
fn bisect_root(mut a: f64, mut b: f64, params: &ThreeRegionParams) -> Option<f64> {
    let tol = ROOT_TOLERANCE;
    let max_iter = 50;

    let mut fa = interface_determinant(a, params);

    for _ in 0..max_iter {
        let mid = (a + b) / 2.0;
        if (b - a) / 2.0 < tol {
            return Some(mid);
        }

        let fmid = interface_determinant(mid, params);
        if fmid.abs() < tol {
            return Some(mid);
        }

        if fa * fmid < 0.0 {
            b = mid;
        } else {
            a = mid;
            fa = fmid;
        }
    }

    Some((a + b) / 2.0)
}

/// Compute eigenfunction coefficients for a given eigenvalue.
///
/// Returns (A, B, C, D, gn, gd) where D is the pole-free region 3
/// coefficient: φ₃(r) = D·(gd·sin(β₃r) + gn·cos(β₃r))/r.
fn compute_coefficients(lambda: f64, params: &ThreeRegionParams) -> (f64, f64, f64, f64, f64, f64) {
    let r_i = params.r_inner;
    let r_o = params.r_outer;

    let beta_1 = lambda / params.alpha_1.sqrt();
    let beta_2 = lambda / params.alpha_2.sqrt();
    let beta_3 = lambda / params.alpha_3.sqrt();

    let k_1 = params.k_1;
    let k_2 = params.k_2;

    let (gn, gd) = dirichlet_gamma_polefree(lambda, params);

    // Set A = 1 (normalization)
    let a = 1.0;

    // From condition 1: A sin(β₁r_i) = B sin(β₂r_i) + C cos(β₂r_i)
    // From condition 2: flux matching gives another equation
    // Solve 2x2 system for B, C

    let s1i = (beta_1 * r_i).sin();
    let c1i = (beta_1 * r_i).cos();
    let s2i = (beta_2 * r_i).sin();
    let c2i = (beta_2 * r_i).cos();
    let s2o = (beta_2 * r_o).sin();
    let c2o = (beta_2 * r_o).cos();
    let s3o = (beta_3 * r_o).sin();
    let c3o = (beta_3 * r_o).cos();

    // Condition 1: A s1i = B s2i + C c2i
    // Condition 2 (flux + inner surface sink):
    //   k₁ A (r_i β₁ c1i - s1i) + h_inner·r_i·A·s1i = k₂ (B(...) + C(...))

    let rhs1 = a * s1i;
    let h_ri = params.h_inner * r_i;
    let rhs2 = k_1 * a * (r_i * beta_1 * c1i - s1i) + h_ri * a * s1i;

    let a11 = s2i;
    let a12 = c2i;
    let a21 = k_2 * (r_i * beta_2 * c2i - s2i);
    let a22 = k_2 * (-r_i * beta_2 * s2i - c2i);

    let det = a11 * a22 - a12 * a21;
    let b = (rhs1 * a22 - rhs2 * a12) / det;
    let c = (a11 * rhs2 - a21 * rhs1) / det;

    // From condition 3 (pole-free form):
    //   B sin(β₂r_o) + C cos(β₂r_o) = D·(gd·sin(β₃r_o) + gn·cos(β₃r_o))
    // where gd·s3o + gn·c3o = sin(β₃(r_o - r_max)), which stays bounded.
    let d = (b * s2o + c * c2o) / (gd * s3o + gn * c3o);

    (a, b, c, d, gn, gd)
}

/// Evaluate eigenfunction at radius r for given mode.
///
/// Region 3 uses the pole-free form: D·(gd·sin(β₃r) + gn·cos(β₃r))/r
/// where gn, gd are stored instead of γ = -tan(β₃R).
#[allow(clippy::too_many_arguments)]
fn eigenfunction_at_r(
    r: f64,
    lambda: f64,
    a: f64,
    b: f64,
    c: f64,
    d: f64,
    gn: f64,
    gd: f64,
    params: &ThreeRegionParams,
) -> f64 {
    let r_i = params.r_inner;
    let r_o = params.r_outer;
    let r_max = params.r_max;

    let beta_1 = lambda / params.alpha_1.sqrt();
    let beta_2 = lambda / params.alpha_2.sqrt();
    let beta_3 = lambda / params.alpha_3.sqrt();

    if r < 1e-10 {
        // At r=0, use L'Hopital: sin(βr)/r → β
        a * beta_1
    } else if r <= r_i {
        // Region 1: liquid
        a * (beta_1 * r).sin() / r
    } else if r <= r_o {
        // Region 2: shell
        (b * (beta_2 * r).sin() + c * (beta_2 * r).cos()) / r
    } else if r <= r_max {
        // Region 3: air (pole-free Dirichlet form)
        d * (gd * (beta_3 * r).sin() + gn * (beta_3 * r).cos()) / r
    } else {
        // Beyond domain: φ(r_max) = 0 by Dirichlet BC, decay as 1/r
        0.0
    }
}

/// Analytical integral of ∫ r sin(βr) dr from a to b
/// Antiderivative: sin(βr)/β² - r cos(βr)/β
pub(crate) fn integral_r_sin(beta: f64, a: f64, b: f64) -> f64 {
    let at_b = (beta * b).sin() / (beta * beta) - b * (beta * b).cos() / beta;
    let at_a = (beta * a).sin() / (beta * beta) - a * (beta * a).cos() / beta;
    at_b - at_a
}

/// Analytical integral of ∫ r cos(βr) dr from a to b
/// Antiderivative: cos(βr)/β² + r sin(βr)/β
fn integral_r_cos(beta: f64, a: f64, b: f64) -> f64 {
    let at_b = (beta * b).cos() / (beta * beta) + b * (beta * b).sin() / beta;
    let at_a = (beta * a).cos() / (beta * beta) + a * (beta * a).sin() / beta;
    at_b - at_a
}

/// Analytical integral of ∫ sin²(βr) dr from a to b
/// Antiderivative: r/2 - sin(2βr)/(4β)
fn integral_sin2(beta: f64, a: f64, b: f64) -> f64 {
    let at_b = b / 2.0 - (2.0 * beta * b).sin() / (4.0 * beta);
    let at_a = a / 2.0 - (2.0 * beta * a).sin() / (4.0 * beta);
    at_b - at_a
}

/// Analytical integral of ∫ cos²(βr) dr from a to b
/// Antiderivative: r/2 + sin(2βr)/(4β)
fn integral_cos2(beta: f64, a: f64, b: f64) -> f64 {
    let at_b = b / 2.0 + (2.0 * beta * b).sin() / (4.0 * beta);
    let at_a = a / 2.0 + (2.0 * beta * a).sin() / (4.0 * beta);
    at_b - at_a
}

/// Analytical integral of ∫ sin(βr)cos(βr) dr from a to b
/// Antiderivative: sin²(βr)/(2β)
fn integral_sin_cos(beta: f64, a: f64, b: f64) -> f64 {
    let at_b = (beta * b).sin().powi(2) / (2.0 * beta);
    let at_a = (beta * a).sin().powi(2) / (2.0 * beta);
    at_b - at_a
}

/// Compute inner product <f, φ> with proper Sturm-Liouville weight for initial condition projection
/// Uses closed-form analytical integrals for exactness at any mode number.
///
/// The inner product is ∫₀^{r_max} (k/α) r² f(r) φ(r) dr
/// where φ(r) = eigenfunction (piecewise sin/cos divided by r)
/// and f(r) = T_init(r) - T_ambient (piecewise constant).
///
/// This reduces to integrals of r*sin(βr) and r*cos(βr) which have exact antiderivatives.
fn inner_product_initial(
    lambda: f64,
    a: f64,
    b: f64,
    c: f64,
    _d: f64,
    params: &ThreeRegionParams,
    _num_points: usize,
) -> f64 {
    let r_i = params.r_inner;
    let r_o = params.r_outer;

    let beta_1 = lambda / params.alpha_1.sqrt();
    let beta_2 = lambda / params.alpha_2.sqrt();

    let w1 = params.k_1 / params.alpha_1; // ρc for liquid
    let w2 = params.k_2 / params.alpha_2; // ρc for shell

    let dt_liquid = params.t_liquid_init - params.t_ambient;
    let dt_shell = params.t_shell_init - params.t_ambient;

    // Region 1 (liquid, 0 to r_i): ∫ w₁ r² (dt_liquid) (A sin(β₁r)/r) dr
    //   = w₁ dt_liquid A ∫₀^{r_i} r sin(β₁r) dr
    let region1 = w1 * dt_liquid * a * integral_r_sin(beta_1, 0.0, r_i);

    // Region 2 (shell, r_i to r_o): ∫ w₂ r² (dt_shell) (B sin(β₂r) + C cos(β₂r))/r dr
    //   = w₂ dt_shell [B ∫ r sin(β₂r) dr + C ∫ r cos(β₂r) dr]
    let region2 = w2
        * dt_shell
        * (b * integral_r_sin(beta_2, r_i, r_o) + c * integral_r_cos(beta_2, r_i, r_o));

    // Region 3 (air, r_o to r_max): f(r) = 0, so no contribution

    region1 + region2
}

/// Compute norm ||φ||² with proper Sturm-Liouville weight using analytical integrals.
///
/// ||φ||² = ∫₀^{r_max} (k/α) r² φ(r)² dr
///
/// Region 3 uses pole-free form: D·(gd·sin + gn·cos)/r, so the integral
/// becomes D²·(gd²·∫sin² + 2·gd·gn·∫sin·cos + gn²·∫cos²).
/// Since gd² + gn² = 1, this is always well-conditioned.
#[allow(clippy::too_many_arguments)]
fn eigenfunction_norm_squared(
    lambda: f64,
    a: f64,
    b: f64,
    c: f64,
    d: f64,
    gn: f64,
    gd: f64,
    params: &ThreeRegionParams,
    _num_points: usize,
) -> f64 {
    let r_i = params.r_inner;
    let r_o = params.r_outer;
    let r_max = params.r_max;

    let beta_1 = lambda / params.alpha_1.sqrt();
    let beta_2 = lambda / params.alpha_2.sqrt();
    let beta_3 = lambda / params.alpha_3.sqrt();

    let w1 = params.k_1 / params.alpha_1;
    let w2 = params.k_2 / params.alpha_2;
    let w3 = params.k_3 / params.alpha_3;

    // Region 1: ∫₀^{r_i} w₁ r² (A sin(β₁r)/r)² dr = w₁ A² ∫ sin²(β₁r) dr
    let norm1 = w1 * a * a * integral_sin2(beta_1, 0.0, r_i);

    // Region 2: ∫_{r_i}^{r_o} w₂ r² ((B sin(β₂r) + C cos(β₂r))/r)² dr
    //   = w₂ [B² ∫sin² + 2BC ∫sin·cos + C² ∫cos²]
    let norm2 = w2
        * (b * b * integral_sin2(beta_2, r_i, r_o)
            + 2.0 * b * c * integral_sin_cos(beta_2, r_i, r_o)
            + c * c * integral_cos2(beta_2, r_i, r_o));

    // Region 3 (pole-free): ∫_{r_o}^{r_max} w₃ r² (D(gd·sin + gn·cos)/r)² dr
    //   = w₃ D² [gd² ∫sin² + 2·gd·gn ∫sin·cos + gn² ∫cos²]
    let norm3 = w3
        * d
        * d
        * (gd * gd * integral_sin2(beta_3, r_o, r_max)
            + 2.0 * gd * gn * integral_sin_cos(beta_3, r_o, r_max)
            + gn * gn * integral_cos2(beta_3, r_o, r_max));

    norm1 + norm2 + norm3
}

/// Compute all eigenmodes for the three-region problem
pub fn compute_eigenmodes(params: &ThreeRegionParams) -> Vec<Eigenmode> {
    let eigenvalues = find_eigenvalues(params, params.num_modes);
    let num_integration_points = 500; // More points for accuracy

    eigenvalues
        .into_iter()
        .map(|lambda| {
            let (a_raw, b_raw, c_raw, d_raw, gn, gd) = compute_coefficients(lambda, params);

            // Compute the norm of the eigenfunction
            let norm_sq = eigenfunction_norm_squared(
                lambda,
                a_raw,
                b_raw,
                c_raw,
                d_raw,
                gn,
                gd,
                params,
                num_integration_points,
            );
            let norm = norm_sq.sqrt();

            // Normalize the coefficients so that ||φ|| = 1
            let (a, b, c, d) = if norm > 1e-20 {
                (a_raw / norm, b_raw / norm, c_raw / norm, d_raw / norm)
            } else {
                (a_raw, b_raw, c_raw, d_raw)
            };

            // Now amplitude = <f, φ> directly (since ||φ|| = 1)
            let amplitude =
                inner_product_initial(lambda, a, b, c, d, params, num_integration_points);

            Eigenmode {
                lambda,
                coeff_a: a,
                coeff_b: b,
                coeff_c: c,
                coeff_d: d,
                gamma_num: gn,
                gamma_den: gd,
                amplitude,
            }
        })
        .collect()
}

/// Compute temperature at (r, t) using the eigenfunction expansion.
/// At t=0, returns the exact initial conditions to avoid Gibbs phenomenon.
fn three_region_temperature(
    r: f64,
    t: f64,
    params: &ThreeRegionParams,
    modes: &[Eigenmode],
) -> f64 {
    // At t=0, return exact initial conditions (the series has Gibbs issues)
    if t <= 0.0 {
        if r <= params.r_inner {
            return params.t_liquid_init;
        } else if r <= params.r_outer {
            return params.t_shell_init;
        } else {
            return params.t_ambient;
        }
    }

    let mut sum = 0.0;

    for mode in modes {
        let phi = eigenfunction_at_r(
            r,
            mode.lambda,
            mode.coeff_a,
            mode.coeff_b,
            mode.coeff_c,
            mode.coeff_d,
            mode.gamma_num,
            mode.gamma_den,
            params,
        );
        let temporal = (-mode.lambda * mode.lambda * t).exp();
        sum += mode.amplitude * phi * temporal;
    }

    params.t_ambient + sum
}

/// Generate a snapshot of the three-region solution.
///
/// The radial grid is uniform with `num_points` entries spanning [0, display_r_max],
/// but the two nearest grid points are snapped to r_inner and r_outer so that
/// the interface radii land exactly on grid nodes.  This ensures that
/// `liquid_average_temperature` (which uses inner_interface_idx) never
/// accidentally includes a shell-temperature point in the liquid average.
fn three_region_snapshot(
    t: f64,
    params: &ThreeRegionParams,
    modes: &[Eigenmode],
    num_points: usize,
    display_r_max: f64,
) -> ThreeRegionSnapshot {
    let mut radii: Vec<f64> = (0..num_points)
        .map(|i| i as f64 / (num_points - 1) as f64 * display_r_max)
        .collect();

    // Snap nearest grid points to the exact interface radii.
    // This prevents liquid_average_temperature from including a shell-temp
    // point in the liquid region (which causes ~7°C error at t=0 for cold mugs).
    for &r_iface in &[params.r_inner, params.r_outer] {
        if let Some((idx, _)) = radii.iter().enumerate().min_by(|(_, a), (_, b)| {
            (**a - r_iface)
                .abs()
                .partial_cmp(&(**b - r_iface).abs())
                .unwrap()
        }) {
            radii[idx] = r_iface;
        }
    }

    let temperatures: Vec<f64> = radii
        .iter()
        .map(|&r| three_region_temperature(r, t, params, modes))
        .collect();

    ThreeRegionSnapshot {
        time_s: t,
        radii,
        temperatures,
        r_inner: params.r_inner,
        r_outer: params.r_outer,
    }
}

/// Snapshots of one constant-coefficient solution: `params` and its `modes`
/// over `duration_s`, evenly spaced. `display_r_max` extends the radial grid
/// beyond params.r_max for chart display. The single-linearisation reference
/// for the kernel tests; `simulate_three_region_relinearised` is the solver.
#[cfg(test)]
pub(crate) fn simulate_three_region_analytical(
    params: &ThreeRegionParams,
    modes: &[Eigenmode],
    duration_s: f64,
    num_snapshots: usize,
    num_radial_points: usize,
    display_r_max: f64,
) -> Vec<ThreeRegionSnapshot> {
    (0..num_snapshots)
        .map(|i| {
            let t = i as f64 / (num_snapshots - 1).max(1) as f64 * duration_s;
            three_region_snapshot(t, params, modes, num_radial_points, display_r_max)
        })
        .collect()
}

/// Number of equal time segments over which `simulate_three_region_relinearised`
/// re-estimates the surface losses.
///
/// Why 10: measured against the numerical solver (ceramic; defaults and the
/// 500 ml / 2 mm / 8 mm cases), 10 segments put every case's extra drinkable
/// minutes within 0.1 min of it; 5 miss the 8 mm case by 0.17 min, and 2–3
/// hit or miss the default case by luck. 15–40 segments move extra minutes
/// by at most 0.06 and leave the crossings' bias (0.2–0.4 min early, from
/// linearising at the liquid average rather than the open surface) in place,
/// while the cost grows linearly with the count.
pub const ANALYTICAL_SEGMENTS: usize = 10;

/// Simpson panels per region when projecting a segment's starting profile
/// onto that segment's eigenmodes.
///
/// Why 100: the run's metrics are insensitive (3 panels already keep extra
/// minutes, penalty and crossings within 0.01 of 250 panels, because a
/// segment's high modes decay within seconds), so the bound is the profile
/// right after a segment boundary: at the panel's maximum of 100 series
/// terms, 60 panels leave a 0.05 K RMS jump there and 100 keep it under
/// 0.01 K (`relinearised_is_continuous_across_segments`).
const PROJECTION_PANELS: usize = 100;

/// One time segment of the re-linearised solution: constant-coefficient
/// parameters and their modes, with amplitudes for the segment's start.
struct Segment {
    params: ThreeRegionParams,
    modes: Vec<Eigenmode>,
}

impl Segment {
    /// Temperature at radius `r` and time `t` since the segment's start.
    fn temperature(&self, r: f64, t: f64) -> f64 {
        three_region_temperature(r, t, &self.params, &self.modes)
    }

    /// Volume-weighted liquid average at time `t` since the segment's start.
    fn liquid_average(&self, t: f64) -> f64 {
        let p = &self.params;
        if t <= 0.0 {
            // Only a run's first segment is evaluated at its own start, where
            // the series has Gibbs error and the start is known exactly.
            return p.t_liquid_init;
        }
        let r_i = p.r_inner;
        let sum: f64 = self
            .modes
            .iter()
            .map(|m| {
                let beta_1 = m.lambda / p.alpha_1.sqrt();
                m.amplitude
                    * m.coeff_a
                    * integral_r_sin(beta_1, 0.0, r_i)
                    * (-m.lambda * m.lambda * t).exp()
            })
            .sum();
        p.t_ambient + 3.0 / (r_i * r_i * r_i) * sum
    }
}

/// r·φ(r) of `mode` in `region` (0 liquid, 1 shell, 2 air): the
/// eigenfunction without its 1/r, finite at r = 0.
fn mode_times_r(mode: &Eigenmode, region: usize, r: f64, params: &ThreeRegionParams) -> f64 {
    match region {
        0 => mode.coeff_a * (mode.lambda / params.alpha_1.sqrt() * r).sin(),
        1 => {
            let beta_2 = mode.lambda / params.alpha_2.sqrt();
            mode.coeff_b * (beta_2 * r).sin() + mode.coeff_c * (beta_2 * r).cos()
        }
        _ => {
            let beta_3 = mode.lambda / params.alpha_3.sqrt();
            mode.coeff_d
                * (mode.gamma_den * (beta_3 * r).sin() + mode.gamma_num * (beta_3 * r).cos())
        }
    }
}

/// A temperature profile sampled on Simpson nodes, one panel set per region
/// with the interfaces on panel edges, so each region's integrand is smooth.
/// Each sample is stored premultiplied by its quadrature weight, the
/// region's capacity k/α and r (the other r of the r² measure is in
/// `mode_times_r`).
struct SampledProfile {
    radii: [Vec<f64>; 3],
    weighted: [Vec<f64>; 3],
}

impl SampledProfile {
    /// Sample `profile(region, r)` (°C) for the geometry of `params`. The
    /// region index lets a discontinuous profile give each side of an
    /// interface its own value.
    fn sample(params: &ThreeRegionParams, profile: impl Fn(usize, f64) -> f64) -> Self {
        let bounds = [
            (0.0, params.r_inner, params.k_1 / params.alpha_1),
            (params.r_inner, params.r_outer, params.k_2 / params.alpha_2),
            (params.r_outer, params.r_max, params.k_3 / params.alpha_3),
        ];
        let n = 2 * PROJECTION_PANELS;
        let region = |idx: usize| {
            let (a, b, capacity) = bounds[idx];
            let h = (b - a) / n as f64;
            (0..=n)
                .map(|j| {
                    let r = a + j as f64 * h;
                    let simpson = if j == 0 || j == n {
                        1.0
                    } else if j % 2 == 1 {
                        4.0
                    } else {
                        2.0
                    };
                    let w = simpson * h / 3.0 * capacity * r;
                    (r, w * (profile(idx, r) - params.t_ambient))
                })
                .unzip::<f64, f64, Vec<f64>, Vec<f64>>()
        };
        let [(r0, w0), (r1, w1), (r2, w2)] = [region(0), region(1), region(2)];
        Self {
            radii: [r0, r1, r2],
            weighted: [w0, w1, w2],
        }
    }

    /// Expansion amplitudes of this profile on `modes`: the capacity-weighted
    /// inner product ⟨profile − T_amb, φ⟩ (the modes are normalised).
    fn project(&self, params: &ThreeRegionParams, modes: &[Eigenmode]) -> Vec<f64> {
        modes
            .iter()
            .map(|m| {
                (0..3)
                    .map(|region| {
                        self.radii[region]
                            .iter()
                            .zip(&self.weighted[region])
                            .map(|(&r, &s)| s * mode_times_r(m, region, r, params))
                            .sum::<f64>()
                    })
                    .sum()
            })
            .collect()
    }
}

/// Expansion amplitudes of `profile` (°C, by region and radius) on `modes`,
/// by Simpson's rule with `PROJECTION_PANELS` panels per region.
#[cfg(test)]
fn project_profile(
    params: &ThreeRegionParams,
    modes: &[Eigenmode],
    profile: impl Fn(usize, f64) -> f64,
) -> Vec<f64> {
    SampledProfile::sample(params, profile).project(params, modes)
}

/// The state a segment starts from: its liquid-average and outer-surface
/// temperatures, and its profile. A run's first segment starts from the
/// mug's uniform three-region state, whose amplitudes are closed-form
/// (`compute_eigenmodes`); later segments start from the previous
/// segment's end, projected onto their own modes.
struct SegmentStart {
    liquid_c: f64,
    outer_c: f64,
    profile: Option<SampledProfile>,
}

impl SegmentStart {
    fn uniform(mug: &ThreeRegionMug) -> Self {
        Self {
            liquid_c: mug.liquid_temp_c,
            outer_c: mug.shell_temp_c,
            profile: None,
        }
    }

    /// The end of `prev`, `len` seconds after its start.
    fn after(prev: &Segment, len: f64) -> Self {
        let p = &prev.params;
        Self {
            liquid_c: prev.liquid_average(len),
            outer_c: prev.temperature(p.r_outer, len),
            profile: Some(SampledProfile::sample(p, |_, r| prev.temperature(r, len))),
        }
    }

    /// The segment from this start with the losses linearised at liquid
    /// temperature `liquid_lin_c` and outer-surface temperature `outer_lin_c`.
    fn solve(&self, mug: &ThreeRegionMug, liquid_lin_c: f64, outer_lin_c: f64) -> Segment {
        let params = ThreeRegionParams::linearised(mug, liquid_lin_c, outer_lin_c);
        let mut modes = compute_eigenmodes(&params);
        if let Some(profile) = &self.profile {
            let amplitudes = profile.project(&params, &modes);
            for (m, a) in modes.iter_mut().zip(amplitudes) {
                m.amplitude = a;
            }
        }
        Segment { params, modes }
    }

    /// The segment of length `len` from this start, by the midpoint rule:
    /// linearise at the start temperatures to predict the segment's end, then
    /// re-solve linearised at the mean of start and predicted end.
    fn solve_midpoint(&self, mug: &ThreeRegionMug, len: f64) -> Segment {
        let predictor = self.solve(mug, self.liquid_c, self.outer_c);
        let liquid_end = predictor.liquid_average(len);
        let outer_end = predictor.temperature(mug.r_outer, len);
        self.solve(
            mug,
            0.5 * (self.liquid_c + liquid_end),
            0.5 * (self.outer_c + outer_end),
        )
    }
}

/// The `segments` equal segments of a `duration_s` run of `mug`, in order.
fn solve_segments(mug: &ThreeRegionMug, segments: usize, duration_s: f64) -> Vec<Segment> {
    let len = duration_s / segments as f64;
    let mut solved: Vec<Segment> = Vec::with_capacity(segments);
    for _ in 0..segments {
        let start = match solved.last() {
            None => SegmentStart::uniform(mug),
            Some(prev) => SegmentStart::after(prev, len),
        };
        solved.push(start.solve_midpoint(mug, len));
    }
    solved
}

/// The segment holding run time `t` (of `count` segments of length `len`)
/// and the time since that segment's start. A boundary time belongs to the
/// segment it ends, so no evaluation lands on a projected series' own start.
fn segment_at(t: f64, len: f64, count: usize) -> (usize, f64) {
    let k = ((t / len).ceil() as usize).saturating_sub(1).min(count - 1);
    (k, t - k as f64 * len)
}

/// Snapshots of a run made of `solved` equal segments spanning `duration_s`.
fn segment_snapshots(
    solved: &[Segment],
    duration_s: f64,
    num_snapshots: usize,
    n_r: usize,
    display_r_max: f64,
) -> Vec<ThreeRegionSnapshot> {
    let len = duration_s / solved.len() as f64;
    (0..num_snapshots)
        .map(|i| {
            let t = i as f64 / (num_snapshots - 1).max(1) as f64 * duration_s;
            let (k, local_t) = segment_at(t, len, solved.len());
            let seg = &solved[k];
            let mut snap =
                three_region_snapshot(local_t, &seg.params, &seg.modes, n_r, display_r_max);
            snap.time_s = t;
            snap
        })
        .collect()
}

/// Simulate `mug` with the heat-equation analytical solver. The run is cut
/// into `segments` equal segments (at least one); each is solved exactly by
/// eigenfunction expansion with its open-top loss and outer radiation
/// linearised by the midpoint rule for that segment, starting from the
/// previous segment's end profile. Snapshots are evenly spaced over
/// `duration_s` on `n_r` radii; `display_r_max` extends the radial grid for
/// chart display.
pub fn simulate_three_region_relinearised(
    mug: &ThreeRegionMug,
    segments: usize,
    duration_s: f64,
    num_snapshots: usize,
    n_r: usize,
    display_r_max: f64,
) -> Vec<ThreeRegionSnapshot> {
    let solved = solve_segments(mug, segments.max(1), duration_s);
    segment_snapshots(&solved, duration_s, num_snapshots, n_r, display_r_max)
}

/// Closed-form liquid average of `simulate_three_region_relinearised` at
/// `num_points` evenly spaced times, free of the radial grid.
#[cfg(test)]
pub(crate) fn relinearised_liquid_average(
    mug: &ThreeRegionMug,
    segments: usize,
    duration_s: f64,
    num_points: usize,
) -> Vec<(f64, f64)> {
    let solved = solve_segments(mug, segments.max(1), duration_s);
    let len = duration_s / solved.len() as f64;
    (0..num_points)
        .map(|i| {
            let t = i as f64 / (num_points - 1).max(1) as f64 * duration_s;
            let (k, local_t) = segment_at(t, len, solved.len());
            (t, solved[k].liquid_average(local_t))
        })
        .collect()
}

/// Closed-form volume-weighted liquid average temperature over time.
///
/// Evaluates  T_avg(t) = T_amb + (3/r³) Σ_n a_n A_n exp(-λ²_n t) ∫₀^r r'sin(β₁r')dr'
/// using `integral_r_sin`, avoiding all radial discretization.
#[cfg(test)]
pub(crate) fn three_region_liquid_average_closed_form(
    params: &ThreeRegionParams,
    modes: &[Eigenmode],
    duration_s: f64,
    num_points: usize,
) -> Vec<(f64, f64)> {
    let r_i = params.r_inner;
    let r_i3 = r_i * r_i * r_i;
    let coeff = 3.0 / r_i3;

    // Pre-compute the radial integral for each mode (time-independent)
    let mode_integrals: Vec<f64> = modes
        .iter()
        .map(|mode| {
            let beta1 = mode.lambda / params.alpha_1.sqrt();
            mode.coeff_a * integral_r_sin(beta1, 0.0, r_i)
        })
        .collect();

    (0..num_points)
        .map(|i| {
            let t = i as f64 / (num_points - 1).max(1) as f64 * duration_s;
            if t <= 0.0 {
                return (t, params.t_liquid_init);
            }
            let sum: f64 = modes
                .iter()
                .zip(mode_integrals.iter())
                .map(|(mode, &integral)| {
                    mode.amplitude * integral * (-mode.lambda * mode.lambda * t).exp()
                })
                .sum();
            (t, params.t_ambient + coeff * sum)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::super::constants::{CELSIUS_TO_KELVIN, STEFAN_BOLTZMANN};
    use super::*;

    #[test]
    fn test_finds_eigenvalues() {
        let params = ThreeRegionParams::default();
        let eigenvalues = find_eigenvalues(&params, 10);
        assert!(
            eigenvalues.len() >= 5,
            "Should find at least 5 eigenvalues, got {}",
            eigenvalues.len()
        );
        // Eigenvalues should be positive and increasing
        for i in 1..eigenvalues.len() {
            assert!(
                eigenvalues[i] > eigenvalues[i - 1],
                "Eigenvalues should be increasing"
            );
        }
    }

    /// Evaluate a stored mode's eigenfunction at `r`.
    fn phi(mode: &Eigenmode, r: f64, params: &ThreeRegionParams) -> f64 {
        eigenfunction_at_r(
            r,
            mode.lambda,
            mode.coeff_a,
            mode.coeff_b,
            mode.coeff_c,
            mode.coeff_d,
            mode.gamma_num,
            mode.gamma_den,
            params,
        )
    }

    /// Initial deviation from ambient in region 0 (liquid), 1 (shell) or 2 (air).
    fn initial_deviation(region: usize, params: &ThreeRegionParams) -> f64 {
        match region {
            0 => params.t_liquid_init - params.t_ambient,
            1 => params.t_shell_init - params.t_ambient,
            _ => 0.0,
        }
    }

    /// Weighted integral ∫₀^{r_max} (k/α)·r²·g(region, r) dr.
    ///
    /// Composite Simpson within each region, so the weight's and the initial
    /// condition's jumps at r_inner and r_outer fall on panel edges and the
    /// integrand is smooth inside every panel.
    fn weighted_integral(params: &ThreeRegionParams, g: impl Fn(usize, f64) -> f64) -> f64 {
        const INTERVALS: usize = 2000;
        let regions = [
            (0.0, params.r_inner, params.k_1 / params.alpha_1),
            (params.r_inner, params.r_outer, params.k_2 / params.alpha_2),
            (params.r_outer, params.r_max, params.k_3 / params.alpha_3),
        ];
        regions
            .iter()
            .enumerate()
            .map(|(region, &(a, b, weight))| {
                let h = (b - a) / INTERVALS as f64;
                let sum: f64 = (0..=INTERVALS)
                    .map(|i| {
                        let r = a + i as f64 * h;
                        let simpson = if i == 0 || i == INTERVALS {
                            1.0
                        } else if i % 2 == 1 {
                            4.0
                        } else {
                            2.0
                        };
                        simpson * r * r * g(region, r)
                    })
                    .sum();
                weight * sum * h / 3.0
            })
            .sum()
    }

    /// ‖f‖²_w of the initial deviation from ambient.
    fn initial_norm_squared(params: &ThreeRegionParams) -> f64 {
        weighted_integral(params, |region, _| {
            initial_deviation(region, params).powi(2)
        })
    }

    #[test]
    fn test_eigenfunction_orthonormality() {
        let params = ThreeRegionParams::default();
        let modes = compute_eigenmodes(&params);
        let n = 10.min(modes.len());

        for i in 0..n {
            for j in i..n {
                let inner = weighted_integral(&params, |_, r| {
                    phi(&modes[i], r, &params) * phi(&modes[j], r, &params)
                });
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (inner - expected).abs() < 1e-6,
                    "<φ_{}, φ_{}>_w = {inner:.3e}, expected {expected}",
                    i + 1,
                    j + 1
                );
            }
        }

        // The stored amplitude is the projection <f, φ>_w of the initial condition.
        for mode in modes.iter().take(n) {
            let projection = weighted_integral(&params, |region, r| {
                initial_deviation(region, &params) * phi(mode, r, &params)
            });
            assert!(
                (projection - mode.amplitude).abs() < 1e-6 * mode.amplitude.abs().max(1.0),
                "λ={:.6}: stored amplitude {:.6} != projection {projection:.6}",
                mode.lambda,
                mode.amplitude
            );
        }
    }

    #[test]
    fn test_parseval_identity() {
        // With orthonormal modes Σ a_n² approaches ‖f‖²_w: the expansion
        // captures the whole initial condition.
        let params = ThreeRegionParams {
            num_modes: 100,
            ..Default::default()
        };
        let modes = compute_eigenmodes(&params);
        assert_eq!(modes.len(), 100);

        let f_norm_sq = initial_norm_squared(&params);
        let sum_amp_sq: f64 = modes.iter().map(|m| m.amplitude * m.amplitude).sum();
        let ratio = sum_amp_sq / f_norm_sq;
        assert!(
            (ratio - 1.0).abs() < 0.02,
            "Σa_n² = {sum_amp_sq:.1} is {:.2}% of ‖f‖² = {f_norm_sq:.1}",
            ratio * 100.0
        );
    }

    #[test]
    fn test_high_frequency_eigenvalues() {
        // The initial condition is a step function, so its coefficients decay
        // like 1/n: |a_n|·n stays within a bounded band at high n, and the
        // partial sums of a_n² never exceed ‖f‖² (Bessel's inequality).
        let params = ThreeRegionParams {
            num_modes: 500,
            ..Default::default()
        };
        let modes = compute_eigenmodes(&params);
        assert_eq!(modes.len(), 500);

        let f_norm_sq = initial_norm_squared(&params);
        let total: f64 = modes.iter().map(|m| m.amplitude * m.amplitude).sum();
        assert!(
            total <= f_norm_sq * (1.0 + 1e-6),
            "Bessel's inequality violated: Σa_n² = {total:.3} > ‖f‖² = {f_norm_sq:.3}"
        );

        // Mode families from the three regions interleave, so |a_n| beats
        // and passes near zero every ~16 modes; the 1/n law holds for its
        // envelope. Take the envelope as the maximum of |a_n| over 32-mode
        // blocks and fit its log-log slope against the block's central n:
        // 1/n decay gives −1, 1/n² gives −2, no decay gives 0.
        let points: Vec<(f64, f64)> = (100..500)
            .collect::<Vec<usize>>()
            .chunks(32)
            .map(|block| {
                let n_mid = 0.5 * (block[0] + block[block.len() - 1]) as f64;
                let envelope = block
                    .iter()
                    .map(|&n| modes[n - 1].amplitude.abs())
                    .fold(0.0, f64::max);
                (n_mid.ln(), envelope.ln())
            })
            .collect();
        let count = points.len() as f64;
        let mean_x = points.iter().map(|p| p.0).sum::<f64>() / count;
        let mean_y = points.iter().map(|p| p.1).sum::<f64>() / count;
        let slope = points
            .iter()
            .map(|(x, y)| (x - mean_x) * (y - mean_y))
            .sum::<f64>()
            / points
                .iter()
                .map(|(x, _)| (x - mean_x).powi(2))
                .sum::<f64>();
        assert!(
            (-1.3..=-0.7).contains(&slope),
            "envelope of |a_n| over n=100..500 decays as n^{slope:.2}; expected n^-1"
        );
    }

    #[test]
    fn test_eigenvalue_scale() {
        // With h_conv = 0 the slowest mode is the liquid and shell losing heat
        // by conduction through the air shell to r_max. Its decay time is the
        // lumped estimate C/G: C the liquid + shell heat capacity, G the
        // conductance 4πk₃/(1/r_outer − 1/r_max) of a spherical air shell.
        let params = ThreeRegionParams::default();
        let eigenvalues = find_eigenvalues(&params, 15);
        assert_eq!(eigenvalues.len(), 15);

        let volume = |a: f64, b: f64| 4.0 / 3.0 * PI * (b.powi(3) - a.powi(3));
        let capacity = params.k_1 / params.alpha_1 * volume(0.0, params.r_inner)
            + params.k_2 / params.alpha_2 * volume(params.r_inner, params.r_outer);
        let conductance = 4.0 * PI * params.k_3 / (1.0 / params.r_outer - 1.0 / params.r_max);
        let tau_lumped = capacity / conductance;
        let tau_1 = 1.0 / eigenvalues[0].powi(2);
        assert!(
            (tau_1 / tau_lumped - 1.0).abs() < 0.05,
            "τ₁ = {tau_1:.0} s, lumped C/G = {tau_lumped:.0} s"
        );

        // Every eigenvalue is a root: the determinant changes sign across it.
        for &lambda in &eigenvalues {
            let eps = 1e-7 * lambda;
            let below = interface_determinant(lambda - eps, &params);
            let above = interface_determinant(lambda + eps, &params);
            assert!(
                below * above < 0.0,
                "det does not change sign across λ={lambda:.6}: {below:.3e}, {above:.3e}"
            );
        }
    }

    #[test]
    fn test_dense_eigenvalue_scan() {
        // A scan far denser than find_eigenvalues' own search brackets every
        // eigenvalue it returns, and finds no root it skipped.
        let params = ThreeRegionParams::default();
        let eigenvalues = find_eigenvalues(&params, params.num_modes);
        let last = *eigenvalues.last().unwrap();

        let lambda_min = 0.001;
        let num_points = 50_000;
        let d_lambda = (last * 1.001 - lambda_min) / num_points as f64;

        let mut sign_changes = Vec::new();
        let mut prev_det = interface_determinant(lambda_min, &params);
        for i in 1..=num_points {
            let lambda = lambda_min + i as f64 * d_lambda;
            let det = interface_determinant(lambda, &params);
            if prev_det.is_finite() && det.is_finite() && prev_det * det < 0.0 {
                sign_changes.push(lambda - d_lambda / 2.0);
            }
            prev_det = det;
        }

        assert_eq!(
            sign_changes.len(),
            eigenvalues.len(),
            "dense scan found {} roots up to λ={last:.4}, find_eigenvalues {}",
            sign_changes.len(),
            eigenvalues.len()
        );
        for (&scanned, &found) in sign_changes.iter().zip(&eigenvalues) {
            assert!(
                (scanned - found).abs() <= d_lambda,
                "dense-scan root λ≈{scanned:.5} does not match find_eigenvalues' λ={found:.5}"
            );
        }
    }

    #[test]
    fn test_eigenfunction_interface_conditions() {
        // compute_coefficients builds every λ's eigenfunction to satisfy
        // temperature and flux continuity at r_inner and temperature
        // continuity at r_outer, so those three are construction checks.
        // The energy balance at r_outer, k₂T₂′ + h_conv·T = k₃T₃′, holds only
        // at an eigenvalue: it holds at each one found, and fails midway
        // between two.
        let params = ThreeRegionParams::default();
        let eigenvalues = find_eigenvalues(&params, 5);

        let r_i = params.r_inner;
        let r_o = params.r_outer;
        // (T₁, T₂) and (k₁T₁′, k₂T₂′) at r_inner; (T₂, T₃) and
        // (k₂T₂′ + h_conv·T₂, k₃T₃′) at r_outer, for λ's coefficients.
        let interface_pairs = |lambda: f64| -> [(f64, f64); 4] {
            let (a, b, c, d, gn, gd) = compute_coefficients(lambda, &params);
            let beta_1 = lambda / params.alpha_1.sqrt();
            let beta_2 = lambda / params.alpha_2.sqrt();
            let beta_3 = lambda / params.alpha_3.sqrt();

            let u1_ri = a * (beta_1 * r_i).sin();
            let u2_ri = b * (beta_2 * r_i).sin() + c * (beta_2 * r_i).cos();
            let u2_ro = b * (beta_2 * r_o).sin() + c * (beta_2 * r_o).cos();
            let u3_ro = d * (gd * (beta_3 * r_o).sin() + gn * (beta_3 * r_o).cos());

            // dT/dr = (r·du/dr − u)/r²
            let du1_ri = a * beta_1 * (beta_1 * r_i).cos();
            let du2_ri = b * beta_2 * (beta_2 * r_i).cos() - c * beta_2 * (beta_2 * r_i).sin();
            let du2_ro = b * beta_2 * (beta_2 * r_o).cos() - c * beta_2 * (beta_2 * r_o).sin();
            let du3_ro = d * beta_3 * (gd * (beta_3 * r_o).cos() - gn * (beta_3 * r_o).sin());
            let grad = |r: f64, u: f64, du: f64| (r * du - u) / (r * r);

            [
                (u1_ri / r_i, u2_ri / r_i),
                (
                    params.k_1 * grad(r_i, u1_ri, du1_ri),
                    params.k_2 * grad(r_i, u2_ri, du2_ri),
                ),
                (u2_ro / r_o, u3_ro / r_o),
                (
                    params.k_2 * grad(r_o, u2_ro, du2_ro) + params.h_conv * u2_ro / r_o,
                    params.k_3 * grad(r_o, u3_ro, du3_ro),
                ),
            ]
        };
        let rel_diff = |(x, y): (f64, f64)| (x - y).abs() / x.abs().max(y.abs()).max(1e-12);
        let names = ["T at r_inner", "flux at r_inner", "T at r_outer"];

        for &lambda in &eigenvalues {
            let pairs = interface_pairs(lambda);
            for (&pair, name) in pairs.iter().zip(names) {
                assert!(
                    rel_diff(pair) < 1e-9,
                    "λ={lambda:.6}: {name} mismatched: {pair:?}"
                );
            }
            // find_eigenvalues brackets each root to ~1e-10 in λ, which leaves
            // a relative imbalance up to ~2e-6; off an eigenvalue it is O(1).
            assert!(
                rel_diff(pairs[3]) < 1e-4,
                "λ={lambda:.6}: energy balance at r_outer mismatched: {:?}",
                pairs[3]
            );
        }
        for pair in eigenvalues.windows(2) {
            let midway = 0.5 * (pair[0] + pair[1]);
            let balance = interface_pairs(midway)[3];
            assert!(
                rel_diff(balance) > 0.1,
                "λ={midway:.6} is not an eigenvalue, yet the r_outer balance holds: {balance:?}"
            );
        }
    }

    #[test]
    fn test_initial_condition_reconstruction() {
        // At t = 0 the series reproduces the piecewise-constant initial
        // condition in every region, away from the jumps where Gibbs
        // oscillation is expected.
        let params = ThreeRegionParams {
            num_modes: 100,
            ..Default::default()
        };
        let modes = compute_eigenmodes(&params);
        let series_at_t0 = |r: f64| -> f64 {
            params.t_ambient
                + modes
                    .iter()
                    .map(|m| m.amplitude * phi(m, r, &params))
                    .sum::<f64>()
        };

        // A truncated series rings for about one wavelength of its highest
        // mode on either side of a jump; that wavelength differs per region.
        let lambda_max = modes.last().unwrap().lambda;
        let margin = |alpha: f64| 2.0 * PI * alpha.sqrt() / lambda_max;
        let regions = [
            (
                0.0,
                params.r_inner - margin(params.alpha_1),
                params.t_liquid_init,
                "liquid",
            ),
            (
                params.r_inner + margin(params.alpha_2),
                params.r_outer - margin(params.alpha_2),
                params.t_shell_init,
                "shell",
            ),
            (
                params.r_outer + margin(params.alpha_3),
                params.r_max,
                params.t_ambient,
                "air",
            ),
        ];
        for (a, b, expected, name) in regions {
            const SAMPLES: usize = 50;
            let sum_sq: f64 = (0..SAMPLES)
                .map(|i| {
                    let r = a + (b - a) * i as f64 / (SAMPLES - 1) as f64;
                    (series_at_t0(r) - expected).powi(2)
                })
                .sum();
            let rms = (sum_sq / SAMPLES as f64).sqrt();
            assert!(
                rms < 2.0,
                "{name}: t=0 reconstruction RMS error {rms:.2} °C"
            );
        }
    }

    #[test]
    fn test_boundary_condition_satisfied() {
        let params = ThreeRegionParams::default();
        let modes = compute_eigenmodes(&params);
        // Each mode vanishes at r_max (Dirichlet), so T(r_max) is ambient.
        let temp = three_region_temperature(params.r_max, 100000.0, &params, &modes);
        assert!(
            (temp - params.t_ambient).abs() < 5.0,
            "Temperature at r_max at late time should be ~ambient, got {}",
            temp
        );
    }

    #[test]
    fn test_temperature_decreases() {
        let params = ThreeRegionParams::default();
        let modes = compute_eigenmodes(&params);
        let t1 = three_region_temperature(0.0, 100.0, &params, &modes);
        let t2 = three_region_temperature(0.0, 500.0, &params, &modes);
        let t3 = three_region_temperature(0.0, 1000.0, &params, &modes);
        // Should be cooling toward ambient
        assert!(
            t2 < t1 || (t1 - params.t_ambient).abs() < 5.0,
            "Should cool: t1={}, t2={}",
            t1,
            t2
        );
        assert!(
            t3 < t2 || (t2 - params.t_ambient).abs() < 5.0,
            "Should cool: t2={}, t3={}",
            t2,
            t3
        );
    }

    #[test]
    fn test_approaches_equilibrium() {
        let params = ThreeRegionParams::default();
        let modes = compute_eigenmodes(&params);
        // With τ₁ ~ 19000s, need t >> τ₁ for equilibrium. At t = 100000s (~28 hours):
        let temp = three_region_temperature(0.0, 100000.0, &params, &modes);
        assert!(
            (temp - params.t_ambient).abs() < 5.0,
            "Should approach ambient at long times, got {}",
            temp
        );
    }

    #[test]
    fn test_solution_at_various_times() {
        // The liquid sits above both shell and room, so its average never rises.
        let params = ThreeRegionParams::default();
        let modes = compute_eigenmodes(&params);
        let avg = three_region_liquid_average_closed_form(&params, &modes, 3600.0, 361);
        for pair in avg.windows(2) {
            let ((t0, temp0), (t1, temp1)) = (pair[0], pair[1]);
            assert!(
                temp1 <= temp0 + 1e-9,
                "liquid average rose from {temp0:.4} °C at {t0:.0} s to {temp1:.4} °C at {t1:.0} s"
            );
        }
        let (_, t_end) = *avg.last().unwrap();
        assert!(
            t_end < params.t_liquid_init && t_end > params.t_ambient,
            "liquid average {t_end:.2} °C at 3600 s should lie between ambient {} °C and \
             initial {} °C",
            params.t_ambient,
            params.t_liquid_init
        );
    }

    #[test]
    fn test_litmus_liquid_temp_early_seconds() {
        // Litmus test: liquid center temperature in the first few seconds.
        //
        // Physics: thermal diffusion length = √(α·t).
        // For water (α = 1.43e-7 m²/s):
        //   t=1s  → √(1.43e-7) ≈ 0.38mm  (liquid radius = 22mm)
        //   t=5s  → √(7.15e-7) ≈ 0.85mm
        //   t=10s → √(1.43e-6) ≈ 1.2mm
        //
        // The center is far from any interface, so it should remain at the
        // initial liquid temperature for at least the first ~10 seconds.
        let params = ThreeRegionParams {
            num_modes: 50,
            ..Default::default()
        };
        let modes = compute_eigenmodes(&params);

        let t_init = params.t_liquid_init; // 95°C

        println!("Litmus test: liquid center temperature in first seconds");
        println!("  Initial liquid temp: {:.1}°C", t_init);

        for &t in &[0.0, 1.0, 5.0, 10.0, 30.0, 60.0] {
            let temp = three_region_temperature(0.0, t, &params, &modes);
            let error = (temp - t_init).abs();
            println!(
                "  t={:5.1}s: T(center) = {:.2}°C, error = {:.2}°C",
                t, temp, error
            );
        }

        // At t=0: should be very close to initial temp
        let t0 = three_region_temperature(0.0, 0.0, &params, &modes);
        assert!(
            (t0 - t_init).abs() < 2.0,
            "At t=0, center should be ~{:.0}°C, got {:.1}°C",
            t_init,
            t0
        );

        // At t=1s: center should barely have changed
        let t1 = three_region_temperature(0.0, 1.0, &params, &modes);
        assert!(
            (t1 - t_init).abs() < 2.0,
            "At t=1s, center should still be ~{:.0}°C, got {:.1}°C",
            t_init,
            t1
        );

        // At t=5s: center should still be very close
        let t5 = three_region_temperature(0.0, 5.0, &params, &modes);
        assert!(
            (t5 - t_init).abs() < 3.0,
            "At t=5s, center should still be ~{:.0}°C, got {:.1}°C",
            t_init,
            t5
        );

        // At t=60s: center may have cooled slightly but should still be hot
        let t60 = three_region_temperature(0.0, 60.0, &params, &modes);
        assert!(
            t60 > t_init - 5.0,
            "At t=60s, center should still be above {:.0}°C, got {:.1}°C",
            t_init - 5.0,
            t60
        );
    }

    #[test]
    fn heat_eq_analytical_h_inner_cools_liquid() {
        // An inner-surface sink must cool the liquid, and no faster than the
        // sink can remove heat. With the sign wrong the interface term was a
        // source, its growing mode was missed, and the expansion collapsed to
        // ambient within minutes.
        const H_INNER: f64 = 20.0;
        let liquid_avg = |h_inner: f64| {
            let params = ThreeRegionParams {
                h_inner,
                ..ThreeRegionParams::default()
            };
            let modes = compute_eigenmodes(&params);
            three_region_liquid_average_closed_form(&params, &modes, 1800.0, 31)
        };
        let sealed = liquid_avg(0.0);
        let open = liquid_avg(H_INNER);

        let (t_end, sealed_end) = *sealed.last().unwrap();
        let open_end = open.last().unwrap().1;
        assert!(
            open_end <= sealed_end - 2.0,
            "h_inner={H_INNER}: liquid {open_end:.2} °C at {t_end} s should be ≥ 2 °C below \
             the h_inner=0 liquid {sealed_end:.2} °C"
        );

        // Energy budget: the sink removes at most h·A·(T_liquid_init − T_amb)
        // per second, and it lowers the wall temperature, which only reduces
        // the liquid's loss to the shell. So the extra drop versus the sealed
        // liquid is bounded by that budget over the liquid's heat capacity.
        let p = ThreeRegionParams::default();
        let area = 4.0 * PI * p.r_inner.powi(2);
        let capacity = p.k_1 / p.alpha_1 * 4.0 / 3.0 * PI * p.r_inner.powi(3);
        let max_power = H_INNER * area * (p.t_liquid_init - p.t_ambient);
        for (&(t, t_sealed), &(_, t_open)) in sealed.iter().zip(open.iter()) {
            let max_extra_drop = max_power * t / capacity;
            assert!(
                t_open >= t_sealed - max_extra_drop - 1e-6,
                "t={t:.0} s: liquid {t_open:.2} °C dropped more than the sink allows \
                 (sealed {t_sealed:.2} °C, max extra drop {max_extra_drop:.2} K)"
            );
        }
    }

    #[test]
    fn test_panel_default_analytical() {
        // Reproduce the EXACT parameters the panel uses with default settings.
        // Defaults: volume=180, wall=5mm, ceramic, liquid=90, shell=80, room=22,
        //           h_conv=20, exposure=0.2, r_max_mult=10, num_modes=50
        let water = super::super::heat_equation_numerical::HeatMaterial::water();
        let ceramic = super::super::heat_equation_numerical::HeatMaterial::ceramic();
        let air = super::super::heat_equation_numerical::HeatMaterial::air();

        let volume_ml = 180.0;
        let inner_r_cm = (3.0 * volume_ml / (4.0 * std::f64::consts::PI)).cbrt();
        let wall_mm = 5.0;
        let outer_r_cm = inner_r_cm + wall_mm / 10.0;
        let r_max_mult = 10.0;
        let h_conv = 20.0;
        let liquid_temp = 90.0;
        let shell_temp = 80.0;
        let room_temp = 22.0;

        // Compute h_conv_eff exactly as run_heat_sim does:
        // Only outer-surface radiation (inner-surface effects destroy eigenvalues)
        let sigma = STEFAN_BOLTZMANN;
        let t_avg_k: f64 = 0.5 * (liquid_temp + room_temp) + CELSIUS_TO_KELVIN;
        let r_inner = inner_r_cm / 100.0;
        let r_outer = outer_r_cm / 100.0;

        let h_rad_outer = 4.0 * ceramic.emissivity * sigma * t_avg_k.powi(3);
        let h_conv_eff = h_conv + h_rad_outer;

        println!("inner_r = {:.4}m, outer_r = {:.4}m", r_inner, r_outer);
        println!("h_rad_outer={:.2}", h_rad_outer);
        println!("h_conv_eff = {:.2} W/m²K", h_conv_eff);

        let params = ThreeRegionParams {
            r_inner,
            r_outer,
            r_max: r_outer * r_max_mult,
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
            h_conv: h_conv_eff,
            h_inner: 0.0,
        };

        let modes = compute_eigenmodes(&params);
        println!("\nFound {} modes", modes.len());

        // Print first 5 modes and last 5
        for (i, m) in modes.iter().enumerate() {
            if i < 5 || i >= modes.len().saturating_sub(5) {
                println!(
                    "  mode {:3}: λ={:.6}, amp={:10.3}, τ={:.2}s, gn={:.4}, gd={:.4}",
                    i + 1,
                    m.lambda,
                    m.amplitude,
                    1.0 / (m.lambda * m.lambda),
                    m.gamma_num,
                    m.gamma_den
                );
            }
        }

        // Check liquid average at several times
        let avg = three_region_liquid_average_closed_form(&params, &modes, 1800.0, 20);
        println!("\nLiquid average temperature:");
        for &(t, temp) in &avg {
            println!("  t={:5.1}min: T_avg = {:.1}°C", t / 60.0, temp);
        }

        // The liquid should start near 90°C and stay above 60°C for some minutes
        assert!(
            avg[0].1 > 85.0,
            "t=0: liquid avg should be ~90°C, got {:.1}",
            avg[0].1
        );
        let above_60_min = avg
            .iter()
            .rfind(|(_, t)| *t >= 60.0)
            .map(|(t, _)| t / 60.0)
            .unwrap_or(0.0);
        println!("\nTime above 60°C: {:.1} min", above_60_min);
        assert!(
            above_60_min > 1.0,
            "Should stay above 60°C for more than 1 minute"
        );
    }

    /// The panel's default mug (180 ml, 5 mm ceramic, coffee 90 °C, room
    /// 22 °C, lid 20% open) with its shell at `shell_c`.
    fn default_mug(shell_c: f64) -> ThreeRegionMug {
        let r_inner = (3.0 * 180.0 / (4.0 * PI)).cbrt() / 100.0;
        let r_outer = r_inner + 0.005;
        let ceramic = super::super::heat_equation_numerical::HeatMaterial::ceramic();
        ThreeRegionMug {
            r_inner,
            r_outer,
            r_max: r_outer * 10.0,
            shell: ceramic,
            liquid_temp_c: 90.0,
            shell_temp_c: shell_c,
            room_temp_c: 22.0,
            h_conv: 20.0,
            exposure: 0.2,
            num_modes: 50,
        }
    }

    fn rms(a: &[f64], b: &[f64]) -> f64 {
        let sum: f64 = a.iter().zip(b).map(|(x, y)| (x - y).powi(2)).sum();
        (sum / a.len() as f64).sqrt()
    }

    #[test]
    fn first_segment_kernel_reproduces_single_linearisation() {
        let mug = default_mug(80.0);
        let duration = 1800.0;
        let params = ThreeRegionParams::for_mug(&mug);
        let modes = compute_eigenmodes(&params);
        let today = simulate_three_region_analytical(&params, &modes, duration, 100, 80, 0.10);

        let kernel = SegmentStart::uniform(&mug).solve(
            &mug,
            0.5 * (mug.liquid_temp_c + mug.room_temp_c),
            mug.liquid_temp_c,
        );
        let segmented = segment_snapshots(&[kernel], duration, 100, 80, 0.10);

        assert_eq!(today.len(), segmented.len());
        for (a, b) in today.iter().zip(&segmented) {
            assert_eq!(a.time_s, b.time_s);
            for (ta, tb) in a.temperatures.iter().zip(&b.temperatures) {
                assert!((ta - tb).abs() < 1e-9, "t={}: {ta} vs {tb}", a.time_s);
            }
        }
    }

    #[test]
    fn projection_reconstructs_uniform_start() {
        // The uniform three-region start has closed-form amplitudes; Simpson
        // projection must give the same series. Compared at 30 s, where the
        // series no longer carries the start's interface discontinuities.
        let mug = default_mug(80.0);
        let params = ThreeRegionParams::for_mug(&mug);
        let exact = compute_eigenmodes(&params);
        let amplitudes = project_profile(&params, &exact, |region, _| match region {
            0 => mug.liquid_temp_c,
            1 => mug.shell_temp_c,
            _ => mug.room_temp_c,
        });
        let mut projected = exact.clone();
        for (m, a) in projected.iter_mut().zip(amplitudes) {
            m.amplitude = a;
        }
        let radii: Vec<f64> = (0..200).map(|i| i as f64 / 199.0 * params.r_max).collect();
        let eval = |modes: &[Eigenmode]| -> Vec<f64> {
            radii
                .iter()
                .map(|&r| three_region_temperature(r, 30.0, &params, modes))
                .collect()
        };
        let err = rms(&eval(&exact), &eval(&projected));
        assert!(err < 0.01, "RMS {err:.4} K");
    }

    #[test]
    fn relinearised_is_continuous_across_segments() {
        // The panel's maximum series terms: the highest modes are the
        // hardest to project.
        let mug = ThreeRegionMug {
            num_modes: 100,
            ..default_mug(22.0)
        };
        let segments = solve_segments(&mug, ANALYTICAL_SEGMENTS, 1800.0);
        let len = 1800.0 / ANALYTICAL_SEGMENTS as f64;
        let radii: Vec<f64> = (0..200).map(|i| i as f64 / 199.0 * mug.r_max).collect();
        for (k, pair) in segments.windows(2).enumerate() {
            let (prev, next) = (&pair[0], &pair[1]);
            let jump = (next.liquid_average(1e-3) - prev.liquid_average(len)).abs();
            assert!(
                jump < 0.05,
                "boundary {k}: liquid average jumps {jump:.4} K"
            );
            let before: Vec<f64> = radii.iter().map(|&r| prev.temperature(r, len)).collect();
            let after: Vec<f64> = radii.iter().map(|&r| next.temperature(r, 1e-3)).collect();
            let err = rms(&before, &after);
            assert!(err < 0.01, "boundary {k}: profile RMS jump {err:.4} K");
        }
    }

    #[test]
    fn relinearised_handles_near_room_start() {
        let mug = ThreeRegionMug {
            liquid_temp_c: 30.0,
            ..default_mug(22.0)
        };
        let snaps =
            simulate_three_region_relinearised(&mug, ANALYTICAL_SEGMENTS, 1800.0, 61, 80, 0.10);
        let avg: Vec<f64> = snaps
            .iter()
            .map(|s| {
                let liquid: Vec<f64> = s
                    .radii
                    .iter()
                    .zip(&s.temperatures)
                    .filter(|(r, _)| **r <= mug.r_inner)
                    .map(|(_, t)| *t)
                    .collect();
                liquid.iter().sum::<f64>() / liquid.len() as f64
            })
            .collect();
        for s in &snaps {
            assert!(
                s.temperatures.iter().all(|t| t.is_finite()),
                "t={}",
                s.time_s
            );
        }
        for (i, w) in avg.windows(2).enumerate() {
            assert!(
                w[1] <= w[0] + 1e-6,
                "snapshot {i}: liquid {:.4} → {:.4} °C",
                w[0],
                w[1]
            );
        }
        let end = *avg.last().unwrap();
        assert!(
            end > mug.room_temp_c - 1e-6 && end < 30.0,
            "end {end:.3} °C"
        );
    }

    #[test]
    fn relinearised_respects_snapshot_times() {
        let mug = default_mug(80.0);
        for duration in [600.0, 7200.0] {
            let n = 37;
            let snaps = simulate_three_region_relinearised(
                &mug,
                ANALYTICAL_SEGMENTS,
                duration,
                n,
                40,
                0.10,
            );
            assert_eq!(snaps.len(), n);
            for (i, s) in snaps.iter().enumerate() {
                let expected = i as f64 / (n - 1) as f64 * duration;
                assert!(
                    (s.time_s - expected).abs() < 1e-9,
                    "{duration} s run, snapshot {i}"
                );
            }
            assert_eq!(snaps[0].temperatures[0], mug.liquid_temp_c);
        }
    }

    #[test]
    fn sealed_lid_has_no_inner_loss() {
        let mug = ThreeRegionMug {
            exposure: 0.0,
            ..default_mug(80.0)
        };
        for (k, seg) in solve_segments(&mug, ANALYTICAL_SEGMENTS, 1800.0)
            .iter()
            .enumerate()
        {
            assert_eq!(seg.params.h_inner, 0.0, "segment {k}");
        }
    }

    #[test]
    fn finds_close_root_pair() {
        // A 350 ml mug with the lid 56% open, its open-top loss linearised at
        // 90 °C: the two slowest modes sit 0.0014 apart, inside one scan
        // step, so the determinant does not change sign between samples.
        // Missing them dropped the coffee to room within one segment.
        let ceramic = super::super::heat_equation_numerical::HeatMaterial::ceramic();
        let r_inner = (3.0 * 350.0e-6 / (4.0 * PI)).cbrt();
        let r_outer = r_inner + 0.005;
        let mug = ThreeRegionMug {
            r_inner,
            r_outer,
            r_max: r_outer * 10.0,
            shell: ceramic,
            liquid_temp_c: 90.0,
            shell_temp_c: 22.0,
            room_temp_c: 22.0,
            h_conv: 20.0,
            exposure: 0.56,
            num_modes: 50,
        };
        let params = ThreeRegionParams::linearised(&mug, 90.0, 22.0);
        let lambdas = find_eigenvalues(&params, 4);
        for (lo, hi) in [(0.0330, 0.0336), (0.0343, 0.0349)] {
            assert!(
                lambdas.iter().any(|&l| (lo..hi).contains(&l)),
                "no eigenvalue in [{lo}, {hi}]: {lambdas:?}"
            );
        }
        assert!(
            lambdas.windows(2).all(|w| w[0] < w[1]),
            "eigenvalues out of order: {lambdas:?}"
        );
    }

    #[test]
    fn exposure_sweep_keeps_every_slow_mode() {
        // 500 ml, 10 mm wall, h_conv 50: as the lid opens, the gap between
        // the two slowest modes of some segment passes through 1% of the
        // fundamental; dropping one of them cooled the coffee 20–30 K in a
        // single segment. Coffee at 10 and 30 min must fall as the lid opens.
        let ceramic = super::super::heat_equation_numerical::HeatMaterial::ceramic();
        let r_inner = (3.0 * 500.0e-6 / (4.0 * PI)).cbrt();
        let r_outer = r_inner + 0.010;
        let at = |exposure: f64| {
            let mug = ThreeRegionMug {
                r_inner,
                r_outer,
                r_max: r_outer * 10.0,
                shell: ceramic,
                liquid_temp_c: 90.0,
                shell_temp_c: 22.0,
                room_temp_c: 22.0,
                h_conv: 50.0,
                exposure,
                num_modes: 50,
            };
            let avg = relinearised_liquid_average(&mug, ANALYTICAL_SEGMENTS, 1800.0, 4);
            (avg[1].1, avg[3].1)
        };
        let sweep: Vec<(f64, f64)> = (0..=25).map(|k| at(k as f64 / 100.0)).collect();
        for (k, w) in sweep.windows(2).enumerate() {
            assert!(
                w[1].0 <= w[0].0 + 1e-6 && w[1].1 <= w[0].1 + 1e-6,
                "lid {}% → {}%: 10 min {:.2} → {:.2} °C, 30 min {:.2} → {:.2} °C",
                k,
                k + 1,
                w[0].0,
                w[1].0,
                w[0].1,
                w[1].1
            );
        }
    }
}

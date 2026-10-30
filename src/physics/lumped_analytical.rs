use super::constants::{CELSIUS_TO_KELVIN, STEFAN_BOLTZMANN};
use super::surface::open_top_loss;
use super::{CoolingDataPoint, MugParameters};

/// Re-linearize every this many minutes. Shorter segments improve accuracy
/// at the cost of more eigenvalue decompositions. 0.2 min (= record interval)
/// gives < 0.15°C agreement with the numerical solver.
const SEGMENT_MINUTES: f64 = 0.2;

/// Linearised environment losses at the given temperatures.
///
/// Returns the coffee's open-top loss as `(h_top_w_per_k, t_ref_c)` — a loss
/// `h_top·(T_coffee − t_ref)` in W, with `t_ref` below room while the coffee
/// evaporates — and the mug's outer-wall coefficient toward room in W/K.
fn linearized_env_coefficients(
    params: &MugParameters,
    t_coffee_c: f64,
    t_mug_c: f64,
) -> (f64, f64, f64) {
    let t_room_c = params.room_temp_c;
    let a_outer = params.outer_surface_area_m2();

    // Open top: convection, radiation and the evaporation tangent at the
    // current coffee temperature, over the open part of the top
    let top = open_top_loss(t_coffee_c, t_room_c, params.h_conv);
    let h_top = params.open_top_area_m2() * top.h;

    // Mug outer walls: convection plus linearised radiation h_rad = 4·ε·σ·T_m³
    let t_m_mug = 0.5 * (t_mug_c + t_room_c) + CELSIUS_TO_KELVIN;
    let h_mug_env =
        (params.h_conv + 4.0 * params.emissivity * STEFAN_BOLTZMANN * t_m_mug.powi(3)) * a_outer;

    (h_top, top.t_ref, h_mug_env)
}

/// Analytical (closed-form) lumped capacitance cooling simulation.
///
/// Uses the same 2-body model as `simulate_lumped_numerical` but solves the
/// coupled ODEs via eigenvalue decomposition of the 2×2 system matrix:
///
/// ```text
/// M = [ -(a+b),  a  ]
///     [   c,   -(c+d)]
/// ```
///
/// where a,b,c,d are rate constants derived from the same heat transfer
/// coefficients used in the numerical solver. The open-top loss pulls the
/// coffee toward a reference temperature below room (evaporation continues at
/// room temperature), so the system is affine, `x′ = M·x + s`, and each segment
/// solves as `x(t) = x_ss + e^{Mt}(x₀ − x_ss)` with `x_ss = −M⁻¹s`. The
/// eigenvalues of M give the two exponential decay modes.
///
/// To handle the nonlinear evaporation and radiation terms, the solver uses
/// **piecewise linearization**: the losses are re-linearized at the segment's
/// start temperatures every `SEGMENT_MINUTES` (evaporation by its tangent),
/// and the eigenvalue problem is re-solved for each segment. This keeps the analytical structure while
/// closely tracking the nonlinear numerical solution.
pub fn simulate_lumped_analytical(
    params: &MugParameters,
    duration_minutes: usize,
) -> Vec<CoolingDataPoint> {
    let record_interval = 0.2; // Same as numerical: every 0.2 minutes
    let mut results = Vec::with_capacity((duration_minutes as f64 / record_interval) as usize + 1);

    // Thermal masses (J/K)
    let c_coffee = params.coffee_thermal_mass();
    let c_mug = params.mug_thermal_mass();

    // Material constants (fixed for the whole simulation)
    let conductance = params.coffee_mug_conductance();
    let t_room = params.room_temp_c;

    // Current temperatures — updated at each segment boundary
    let mut cur_coffee = params.coffee_temp_c;
    let mut cur_mug = if params.preheated {
        params.preheat_temp_c
    } else {
        t_room
    };

    // Record initial state
    results.push(CoolingDataPoint {
        time_minutes: 0.0,
        coffee_temp_c: cur_coffee,
        mug_temp_c: cur_mug,
    });

    let mut time_min = 0.0;
    let end_min = duration_minutes as f64;

    while time_min < end_min - 1e-9 {
        // Determine this segment's duration
        let seg_end = (time_min + SEGMENT_MINUTES).min(end_min);

        // Linearize at current temperatures
        let (h_top, t_ref_top, h_mug_env) =
            linearized_env_coefficients(params, cur_coffee, cur_mug);

        // Rate constants (1/s)
        let a = conductance / c_coffee;
        let b = h_top / c_coffee;
        let c = conductance / c_mug;
        let d = h_mug_env / c_mug;

        // Relative to room, x = (u, v) obeys x′ = M·x + (b·u_ref, 0): the open
        // top pulls the coffee toward t_ref, not room. Its steady state
        // x_ss = −M⁻¹·(b·u_ref, 0); the deviation x − x_ss decays by M's modes.
        let source = b * (t_ref_top - t_room);
        let (u_ss, v_ss) = if source == 0.0 {
            (0.0, 0.0)
        } else {
            let det = (a + b) * (c + d) - a * c;
            (source * (c + d) / det, source * c / det)
        };

        // Eigenvalues of M = [-(a+b), a; c, -(c+d)]
        let trace = a + b + c + d;
        let diff = (a + b) - (c + d);
        let discriminant = (diff * diff + 4.0 * a * c).sqrt();
        let lambda1 = 0.5 * (-trace + discriminant); // slower decay
        let lambda2 = 0.5 * (-trace - discriminant); // faster decay

        // Eigenvector ratios
        let r1 = (lambda1 + a + b) / a;
        let r2 = (lambda2 + a + b) / a;

        // Initial deviation from the segment's steady state
        let u0 = cur_coffee - t_room - u_ss;
        let v0 = cur_mug - t_room - v_ss;

        // Solve for coefficients
        let coeff2 = (v0 - u0 * r1) / (r2 - r1);
        let coeff1 = u0 - coeff2;

        // Emit data points within this segment
        let mut t_local = time_min + record_interval;
        while t_local <= seg_end + 1e-9 && t_local <= end_min + 1e-9 {
            let dt_sec = (t_local - time_min) * 60.0;
            let e1 = (lambda1 * dt_sec).exp();
            let e2 = (lambda2 * dt_sec).exp();

            let coffee_temp = coeff1 * e1 + coeff2 * e2 + t_room + u_ss;
            let mug_temp = coeff1 * r1 * e1 + coeff2 * r2 * e2 + t_room + v_ss;

            results.push(CoolingDataPoint {
                time_minutes: (t_local * 10.0).round() / 10.0, // round to 0.1 to avoid float drift
                coffee_temp_c: coffee_temp,
                mug_temp_c: mug_temp,
            });

            t_local += record_interval;
        }

        // Advance to segment boundary — evaluate at exact seg_end
        let dt_seg_sec = (seg_end - time_min) * 60.0;
        let e1 = (lambda1 * dt_seg_sec).exp();
        let e2 = (lambda2 * dt_seg_sec).exp();
        cur_coffee = coeff1 * e1 + coeff2 * e2 + t_room + u_ss;
        cur_mug = coeff1 * r1 * e1 + coeff2 * r2 * e2 + t_room + v_ss;

        time_min = seg_end;
    }

    results
}

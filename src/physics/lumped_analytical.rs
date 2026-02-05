use super::{MugParameters, CoolingDataPoint};

const STEFAN_BOLTZMANN: f64 = 5.670_374e-8; // W/(m²·K⁴)
const WATER_EMISSIVITY: f64 = 0.96;
const LATENT_HEAT_VAPORIZATION: f64 = 2.26e6; // J/kg
const R_VAPOR: f64 = 461.5; // J/(kg·K)
const RHO_AIR: f64 = 1.18; // kg/m³
const CP_AIR: f64 = 1005.0; // J/(kg·K)
const RELATIVE_HUMIDITY: f64 = 0.50;
const EVAP_ENHANCEMENT: f64 = 2.5;

/// Re-linearize every this many minutes. Shorter segments improve accuracy
/// at the cost of more eigenvalue decompositions. 0.2 min (= record interval)
/// gives < 0.15°C agreement with the numerical solver.
const SEGMENT_MINUTES: f64 = 0.2;

fn saturation_vapor_pressure(t_celsius: f64) -> f64 {
    610.78 * (17.27 * t_celsius / (t_celsius + 237.3)).exp()
}

/// Compute linearized environment heat loss coefficients at the given temperatures.
/// Returns (h_coffee_env, h_mug_env) in W/K.
fn linearized_env_coefficients(
    t_coffee_c: f64,
    t_mug_c: f64,
    t_room_c: f64,
    h_conv: f64,
    mug_emissivity: f64,
    a_top: f64,
    a_outer: f64,
    exposure: f64,
) -> (f64, f64) {
    let t_room_k = t_room_c + 273.15;
    let t_coffee_k = t_coffee_c + 273.15;
    let t_mug_k = t_mug_c + 273.15;

    // Base convection — top surface scaled by lid opening
    let mut h_coffee_env = exposure * h_conv * a_top;
    let mut h_mug_env = h_conv * a_outer;

    // Linearized radiation: h_rad = 4 * epsilon * sigma * T_m^3
    let t_m_coffee = 0.5 * (t_coffee_k + t_room_k);
    h_coffee_env += exposure * 4.0 * WATER_EMISSIVITY * STEFAN_BOLTZMANN * t_m_coffee.powi(3) * a_top;

    let t_m_mug = 0.5 * (t_mug_k + t_room_k);
    h_mug_env += 4.0 * mug_emissivity * STEFAN_BOLTZMANN * t_m_mug.powi(3) * a_outer;

    // Linearized evaporation — scaled by lid opening
    let h_mass = h_conv / (RHO_AIR * CP_AIR);
    let t_ref = t_coffee_c; // linearize at current coffee temperature
    let t_ref_k = t_ref + 273.15;
    let rho_v_ref = saturation_vapor_pressure(t_ref) / (R_VAPOR * t_ref_k);
    let rho_v_ambient = RELATIVE_HUMIDITY * saturation_vapor_pressure(t_room_c) / (R_VAPOR * t_room_k);
    if rho_v_ref > rho_v_ambient {
        let q_evap_ref = EVAP_ENHANCEMENT * h_mass * a_top * LATENT_HEAT_VAPORIZATION
            * (rho_v_ref - rho_v_ambient);
        let delta_t = t_ref - t_room_c;
        if delta_t > 1.0 {
            h_coffee_env += exposure * q_evap_ref / delta_t;
        }
    }

    (h_coffee_env, h_mug_env)
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
/// coefficients used in the numerical solver. The eigenvalues give the two
/// exponential decay modes, and the solution is a superposition of these modes.
///
/// To handle the nonlinear evaporation and radiation terms, the solver uses
/// **piecewise linearization**: the coefficients are re-linearized at the
/// current temperature every `SEGMENT_MINUTES`, and the eigenvalue problem
/// is re-solved for each segment. This keeps the analytical structure while
/// closely tracking the nonlinear numerical solution.
pub fn simulate_lumped_analytical(params: &MugParameters, duration_minutes: usize) -> Vec<CoolingDataPoint> {
    let record_interval = 0.2; // Same as numerical: every 0.2 minutes
    let mut results = Vec::with_capacity((duration_minutes as f64 / record_interval) as usize + 1);

    // Thermal masses (J/K)
    let c_coffee = params.coffee_thermal_mass();
    let c_mug = params.mug_thermal_mass();

    // Geometry and material constants (fixed for the whole simulation)
    let k_eff = 2.0 * params.material.thermal_conductivity();
    let a_top = params.top_surface_area_m2();
    let a_outer = params.outer_surface_area_m2();
    let mug_emissivity = params.emissivity;
    let t_room = params.room_temp_c;

    // Current temperatures — updated at each segment boundary
    let mut cur_coffee = params.coffee_temp_c;
    let mut cur_mug = if params.preheated { params.preheat_temp_c } else { t_room };

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
        let (h_coffee_env, h_mug_env) = linearized_env_coefficients(
            cur_coffee, cur_mug, t_room, params.h_conv,
            mug_emissivity, a_top, a_outer, params.exposure_pct,
        );

        // Rate constants (1/s)
        let a = k_eff / c_coffee;
        let b = h_coffee_env / c_coffee;
        let c = k_eff / c_mug;
        let d = h_mug_env / c_mug;

        // Eigenvalues of M = [-(a+b), a; c, -(c+d)]
        let trace = a + b + c + d;
        let diff = (a + b) - (c + d);
        let discriminant = (diff * diff + 4.0 * a * c).sqrt();
        let lambda1 = 0.5 * (-trace + discriminant); // slower decay
        let lambda2 = 0.5 * (-trace - discriminant); // faster decay

        // Eigenvector ratios
        let r1 = (lambda1 + a + b) / a;
        let r2 = (lambda2 + a + b) / a;

        // Initial conditions for this segment (relative to room temp)
        let u0 = cur_coffee - t_room;
        let v0 = cur_mug - t_room;

        // Solve for coefficients
        let coeff2 = (v0 - u0 * r1) / (r2 - r1);
        let coeff1 = u0 - coeff2;

        // Emit data points within this segment
        let mut t_local = time_min + record_interval;
        while t_local <= seg_end + 1e-9 && t_local <= end_min + 1e-9 {
            let dt_sec = (t_local - time_min) * 60.0;
            let e1 = (lambda1 * dt_sec).exp();
            let e2 = (lambda2 * dt_sec).exp();

            let coffee_temp = coeff1 * e1 + coeff2 * e2 + t_room;
            let mug_temp = coeff1 * r1 * e1 + coeff2 * r2 * e2 + t_room;

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
        cur_coffee = coeff1 * e1 + coeff2 * e2 + t_room;
        cur_mug = coeff1 * r1 * e1 + coeff2 * r2 * e2 + t_room;

        time_min = seg_end;
    }

    results
}

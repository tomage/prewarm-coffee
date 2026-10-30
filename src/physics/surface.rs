//! Heat loss from the open coffee surface to the room, shared by every solver.
//!
//! The open top loses heat by convection, radiation and evaporation. All
//! three are per unit area (W/m²); callers scale by the exposed area and the
//! lid opening. Evaporation is driven by the vapour-density difference between
//! the saturated surface and the room air, so it has no cutoff: coffee keeps
//! evaporating at and below room temperature and can settle a few degrees
//! below room. `EVAP_ENHANCEMENT` scales mass transfer but not convection, so
//! the open surface's zero-flux temperature (about 14 °C at a 22 °C room,
//! 50% RH) sits below the wet-bulb temperature (15.4 °C); that follows from
//! the assumed enhancement.
//!
//! Exposure is lid openness, 0 (sealed lid) to 1 (open mug with no lid); the
//! mug is always present. At 1 the open area is the top of the lumped model's
//! cylinder, and `sphere_open_fraction` maps it onto the heat-equation sphere.

use super::constants::{
    saturation_vapor_pressure, CELSIUS_TO_KELVIN, CP_AIR, EVAP_ENHANCEMENT,
    LATENT_HEAT_VAPORIZATION, RELATIVE_HUMIDITY, RHO_AIR, R_VAPOR, STEFAN_BOLTZMANN,
    WATER_EMISSIVITY,
};

/// Half-width of the central difference for the evaporation tangent, in K.
const EVAP_DERIVATIVE_STEP: f64 = 0.01;

/// Smallest temperature excess over room, in K, at which
/// `three_region_h_inner` takes the secant toward room.
///
/// Assumption: 5 K is a modelling choice, not a physical quantity.
/// Evaporation is positive at room temperature, so the secant `q/ΔT` grows
/// without bound as ΔT → 0. The analytical solver can only express a loss
/// toward room, so near or below room it linearises at room + this excess
/// instead, keeping h_inner finite, positive and continuous in the liquid
/// temperature. 5 K is small next to typical coffee-to-room differences
/// (50–70 K), so the floor only matters for coffee that starts near room.
const MIN_SECANT_EXCESS_K: f64 = 5.0;

/// Open area of a fully open mug as a share of the heat-equation sphere's
/// inner surface.
///
/// The lumped model's mug is a cylinder with h = 2r whose open top is πr_cyl².
/// The heat-equation model is a sphere of the same volume with inner surface
/// 4πr_sphere². For volume V:
///   r_cyl = (V/(2π))^(1/3),  r_sphere = (3V/(4π))^(1/3)
///   ratio = πr_cyl² / (4πr_sphere²) = (2/3)^(2/3) / 4 ≈ 0.191
pub const LID_TO_SPHERE_AREA_RATIO: f64 = 0.191;

/// Share of the heat-equation sphere's liquid surface that is open to the
/// room at lid openness `exposure` (0–1).
pub fn sphere_open_fraction(exposure: f64) -> f64 {
    exposure * LID_TO_SPHERE_AREA_RATIO
}

/// Inner-surface loss coefficient `h_inner` for the three-region analytical
/// solver, in W/(m²·K) over the whole liquid sphere surface.
///
/// The eigenfunction expansion takes one constant coefficient toward room,
/// so the open-top loss is linearised as a secant toward room at liquid
/// temperature `t_lin_c`, and scaled by the open share of the sphere. When
/// `t_lin_c` is less than `MIN_SECANT_EXCESS_K` above room the secant is
/// taken at room plus that excess. The solver cannot settle below room.
pub fn three_region_h_inner(exposure: f64, t_lin_c: f64, room_temp_c: f64, h_conv: f64) -> f64 {
    let t_lin = t_lin_c.max(room_temp_c + MIN_SECANT_EXCESS_K);
    sphere_open_fraction(exposure) * open_top_secant_h(t_lin, room_temp_c, h_conv)
}

/// Water-vapour density of saturated air at `t_c`, in kg/m³.
fn saturated_vapor_density(t_c: f64) -> f64 {
    saturation_vapor_pressure(t_c) / (R_VAPOR * (t_c + CELSIUS_TO_KELVIN))
}

/// Evaporative heat flux from a water surface at `t_surf_c` into room air, in W/m².
///
/// Mass transfer follows the Lewis analogy with the buoyancy enhancement
/// `EVAP_ENHANCEMENT`; the driving force is the saturated vapour density at
/// the surface minus the room's vapour density at `RELATIVE_HUMIDITY`.
/// Positive means the surface loses heat. It is positive at room temperature.
pub fn evaporation_flux(t_surf_c: f64, t_room_c: f64, h_conv: f64) -> f64 {
    let h_mass = h_conv / (RHO_AIR * CP_AIR);
    let rho_v_room = RELATIVE_HUMIDITY * saturated_vapor_density(t_room_c);
    EVAP_ENHANCEMENT
        * h_mass
        * LATENT_HEAT_VAPORIZATION
        * (saturated_vapor_density(t_surf_c) - rho_v_room)
}

/// Exact (nonlinear) open-top heat flux in W/m²: convection, T⁴ radiation and evaporation.
pub fn open_top_flux(t_surf_c: f64, t_room_c: f64, h_conv: f64) -> f64 {
    let t_surf_k = t_surf_c + CELSIUS_TO_KELVIN;
    let t_room_k = t_room_c + CELSIUS_TO_KELVIN;
    h_conv * (t_surf_c - t_room_c)
        + WATER_EMISSIVITY * STEFAN_BOLTZMANN * (t_surf_k.powi(4) - t_room_k.powi(4))
        + evaporation_flux(t_surf_c, t_room_c, h_conv)
}

/// A surface loss linear in the surface temperature: `q = h·(T − t_ref)` in W/m².
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RobinLoss {
    /// Heat-transfer coefficient in W/(m²·K).
    pub h: f64,
    /// Temperature the loss drives the surface toward, in °C. Below room when
    /// evaporation is active.
    pub t_ref: f64,
}

/// Linearised open-top loss at surface temperature `t_surf_c`.
///
/// Convection `h_conv·(T − room)` and radiation `4εσT_m³·(T − room)` (T_m the
/// mean of surface and room, in K) are already linear through room
/// temperature. Evaporation uses its tangent at `t_surf_c`,
/// `q₀ + q′·(T − T₀)`, so the combined loss stays finite as T → room.
pub fn open_top_loss(t_surf_c: f64, t_room_c: f64, h_conv: f64) -> RobinLoss {
    let h_through_room = h_conv + radiation_h(t_surf_c, t_room_c);
    let q0 = evaporation_flux(t_surf_c, t_room_c, h_conv);
    let dq = (evaporation_flux(t_surf_c + EVAP_DERIVATIVE_STEP, t_room_c, h_conv)
        - evaporation_flux(t_surf_c - EVAP_DERIVATIVE_STEP, t_room_c, h_conv))
        / (2.0 * EVAP_DERIVATIVE_STEP);
    let h = h_through_room + dq;
    if h == 0.0 {
        return RobinLoss {
            h: 0.0,
            t_ref: t_room_c,
        };
    }
    // h·(T − t_ref) = h_through_room·(T − room) + q₀ + q′·(T − T₀)
    let t_ref = (h_through_room * t_room_c + dq * t_surf_c - q0) / h;
    RobinLoss { h, t_ref }
}

/// Secant coefficient of the open-top loss toward room: total loss at
/// `t_lin_c` divided by `t_lin_c − t_room_c`, in W/(m²·K).
///
/// For solvers that can only express a loss toward room temperature. When
/// the surface is at room temperature the secant is undefined and the
/// tangent coefficient of `open_top_loss` is returned. Just above room the
/// evaporation term makes the secant large, and below room it can be
/// negative, so callers linearise at a temperature above room.
pub fn open_top_secant_h(t_lin_c: f64, t_room_c: f64, h_conv: f64) -> f64 {
    let delta_t = t_lin_c - t_room_c;
    if delta_t.abs() < 1e-6 {
        return open_top_loss(t_lin_c, t_room_c, h_conv).h;
    }
    h_conv + radiation_h(t_lin_c, t_room_c) + evaporation_flux(t_lin_c, t_room_c, h_conv) / delta_t
}

/// Linearised water-surface radiation coefficient `4εσT_m³`, in W/(m²·K).
fn radiation_h(t_surf_c: f64, t_room_c: f64) -> f64 {
    let t_m = 0.5 * (t_surf_c + t_room_c) + CELSIUS_TO_KELVIN;
    4.0 * WATER_EMISSIVITY * STEFAN_BOLTZMANN * t_m.powi(3)
}

#[cfg(test)]
mod tests {
    use super::*;

    const ROOM: f64 = 22.0;
    const H_CONV: f64 = 10.0;

    /// Convection plus linearised radiation, the non-evaporative part of `open_top_loss`.
    fn conv_rad_flux(t_surf: f64) -> f64 {
        let t_m = 0.5 * (t_surf + ROOM) + CELSIUS_TO_KELVIN;
        let h_rad = 4.0 * WATER_EMISSIVITY * STEFAN_BOLTZMANN * t_m.powi(3);
        (H_CONV + h_rad) * (t_surf - ROOM)
    }

    #[test]
    fn evaporation_flux_positive_at_room_temp() {
        assert!(evaporation_flux(ROOM, ROOM, H_CONV) > 0.0);
    }

    #[test]
    fn open_top_loss_matches_flux_at_linearisation_point() {
        for t0 in [30.0, 60.0, 95.0] {
            let loss = open_top_loss(t0, ROOM, H_CONV);
            let linear = loss.h * (t0 - loss.t_ref);
            let expected = conv_rad_flux(t0) + evaporation_flux(t0, ROOM, H_CONV);
            assert!(
                ((linear - expected) / expected).abs() < 1e-6,
                "T0={t0}: h·(T0 − t_ref) = {linear}, flux = {expected}"
            );
        }
    }

    #[test]
    fn open_top_loss_finite_near_room() {
        for t0 in [21.9, 22.0, 22.1] {
            let loss = open_top_loss(t0, ROOM, H_CONV);
            assert!(loss.h.is_finite() && loss.h > 0.0, "T0={t0}: h={}", loss.h);
            assert!(loss.t_ref.is_finite(), "T0={t0}: t_ref={}", loss.t_ref);
            assert!(loss.t_ref < ROOM, "T0={t0}: t_ref={} ≥ room", loss.t_ref);
        }
    }

    #[test]
    fn open_top_loss_zero_coefficients() {
        let loss = open_top_loss(60.0, ROOM, 0.0);
        assert!(loss.h.is_finite() && loss.t_ref.is_finite());
        // No air flow means no convection and no evaporation; radiation remains.
        assert!((loss.t_ref - ROOM).abs() < 1e-9);
        assert_eq!(evaporation_flux(60.0, ROOM, 0.0), 0.0);
        let secant = open_top_secant_h(60.0, ROOM, 0.0);
        assert!(secant.is_finite() && secant > 0.0);
    }

    #[test]
    fn lid_to_sphere_ratio_matches_geometry() {
        let exact = (2.0_f64 / 3.0).powf(2.0 / 3.0) / 4.0;
        assert!((LID_TO_SPHERE_AREA_RATIO - exact).abs() < 1e-3);
        assert_eq!(sphere_open_fraction(0.0), 0.0);
        assert_eq!(sphere_open_fraction(1.0), LID_TO_SPHERE_AREA_RATIO);
    }

    #[test]
    fn three_region_h_inner_is_open_share_of_secant_at_t_lin() {
        let h = three_region_h_inner(0.5, 56.0, ROOM, H_CONV);
        let expected = 0.5 * LID_TO_SPHERE_AREA_RATIO * open_top_secant_h(56.0, ROOM, H_CONV);
        assert!(((h - expected) / expected).abs() < 1e-12);
        assert_eq!(three_region_h_inner(0.0, 56.0, ROOM, H_CONV), 0.0);
    }

    #[test]
    fn three_region_h_inner_finite_near_and_below_room() {
        let floor = three_region_h_inner(1.0, ROOM + MIN_SECANT_EXCESS_K, ROOM, H_CONV);
        for liquid in [ROOM - 5.0, ROOM, ROOM + 0.1, ROOM + MIN_SECANT_EXCESS_K] {
            let h = three_region_h_inner(1.0, liquid, ROOM, H_CONV);
            assert!(h.is_finite() && h > 0.0, "liquid={liquid}: h_inner={h}");
            assert_eq!(
                h, floor,
                "liquid={liquid}: below the floor h_inner is constant"
            );
        }
    }

    #[test]
    fn open_top_secant_h_is_total_flux_over_delta_t() {
        let t = 80.0;
        let expected = (conv_rad_flux(t) + evaporation_flux(t, ROOM, H_CONV)) / (t - ROOM);
        let secant = open_top_secant_h(t, ROOM, H_CONV);
        assert!(((secant - expected) / expected).abs() < 1e-12);
        // At room temperature the secant is undefined; it falls back to the tangent.
        let at_room = open_top_secant_h(ROOM, ROOM, H_CONV);
        assert_eq!(at_room, open_top_loss(ROOM, ROOM, H_CONV).h);
    }
}

//! Physical constants and named modelling assumptions shared by every solver.
//!
//! Each value is defined exactly once here. Constants marked "assumption" are
//! modelling choices rather than measured physical quantities; the doc comment
//! says why the value was chosen.

/// Stefan-Boltzmann constant in W/(m²·K⁴).
pub const STEFAN_BOLTZMANN: f64 = 5.670_374e-8;

/// Emissivity of a water surface (dimensionless).
pub const WATER_EMISSIVITY: f64 = 0.96;

/// Latent heat of vaporization of water in J/kg.
pub const LATENT_HEAT_VAPORIZATION: f64 = 2.26e6;

/// Specific gas constant of water vapor in J/(kg·K).
pub const R_VAPOR: f64 = 461.5;

/// Density of air at ~25 °C in kg/m³.
pub const RHO_AIR: f64 = 1.18;

/// Specific heat of air at ~25 °C in J/(kg·K).
pub const CP_AIR: f64 = 1005.0;

/// Thermal conductivity of air at ~25 °C in W/(m·K).
///
/// Used for boundary-layer profile rendering in the profile chart and for
/// convection calculations at the mug outer surface.
pub const K_AIR: f64 = 0.026;

/// Relative humidity of the room air (0–1).
///
/// Assumption: 50% is a typical indoor value; the app has no humidity input.
pub const RELATIVE_HUMIDITY: f64 = 0.50;

/// Buoyancy enhancement of evaporative mass transfer over an open cup.
///
/// Assumption: the Lewis analogy underestimates mass transfer because the
/// rising thermal plume above the cup carries vapor away faster than the
/// flat-plate analogy predicts; 2.5 is an order-of-magnitude estimate for an
/// open cup, not a fitted value. It scales mass transfer only, not the
/// convective heat transfer beside it, which puts the open surface's
/// zero-flux temperature a little below the wet-bulb temperature.
pub const EVAP_ENHANCEMENT: f64 = 2.5;

/// Density of water in kg/m³ (1 mL of water is 1 g).
pub const WATER_DENSITY: f64 = 1000.0;

/// Specific heat of water in J/(kg·K).
pub const WATER_SPECIFIC_HEAT: f64 = 4186.0;

/// Effective thermal conductivity of the coffee in W/(m·K).
///
/// Assumption: molecular conductivity of water is 0.6 W/(m·K), but hot coffee
/// has vigorous buoyancy-driven convection (Ra ~ 10⁹) that keeps the bulk
/// nearly well mixed. About 10× molecular approximates that internal mixing
/// without a full Navier-Stokes solve.
pub const WATER_EFFECTIVE_CONDUCTIVITY: f64 = 6.0;

/// Offset between degrees Celsius and kelvin.
pub const CELSIUS_TO_KELVIN: f64 = 273.15;

/// Saturation vapor pressure of water via the Tetens formula, in Pa.
pub fn saturation_vapor_pressure(t_celsius: f64) -> f64 {
    610.78 * (17.27 * t_celsius / (t_celsius + 237.3)).exp()
}

/// Water-side film coefficient at the inner mug wall in W/(m²·K).
///
/// Assumption: natural convection of hot water on a vertical wall,
/// Churchill–Chu / vertical-plate correlation (Incropera & DeWitt, ch. 9).
/// With ΔT ≈ 5 K and L ≈ 0.1 m this gives Nu ≈ 80, h ≈ 500 W/(m²·K), the same
/// order as the heat-equation model's effective water conductivity over the
/// coffee radius.
pub const H_WATER_SIDE: f64 = 500.0;

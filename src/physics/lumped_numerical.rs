use super::constants::{CELSIUS_TO_KELVIN, STEFAN_BOLTZMANN};
use super::surface::open_top_flux;
use super::MugParameters;

/// Maximum simulation duration in minutes for the lumped model.
/// Used to ensure metrics aren't clamped by chart duration.
pub const MAX_SIMULATION_MINUTES: usize = 120;

/// A single data point in the cooling simulation
#[derive(Debug, Clone, PartialEq)]
pub struct CoolingDataPoint {
    pub time_minutes: f64,
    pub coffee_temp_c: f64,
    pub mug_temp_c: f64,
}

impl CoolingDataPoint {
    /// Extract (time_minutes, coffee_temp_c) pairs from a slice.
    pub fn coffee_series(data: &[Self]) -> Vec<(f64, f64)> {
        data.iter()
            .map(|p| (p.time_minutes, p.coffee_temp_c))
            .collect()
    }

    /// Extract (time_minutes, mug_temp_c) pairs from a slice.
    pub fn mug_series(data: &[Self]) -> Vec<(f64, f64)> {
        data.iter()
            .map(|p| (p.time_minutes, p.mug_temp_c))
            .collect()
    }
}

/// Simulate the cooling of coffee in a mug over time.
///
/// # Physics Model
///
/// This simulation uses Newton's Law of Cooling with a lumped capacitance model,
/// which assumes uniform temperature within each body (coffee and mug). This is
/// a reasonable approximation for typical mug geometries and time scales.
///
/// ## Heat Transfer Equations
///
/// The model tracks three heat flows:
///
/// 1. **Coffee → Mug (conduction)**: `Q = G × (T_coffee - T_mug)`
///    where G is the water-side film and half the wall in series
///    (`MugParameters::coffee_mug_conductance`)
///
/// 2. **Coffee → Environment (open top)**: `Q = A_open × q(T_coffee)`
///    Convection, T⁴ radiation and evaporation from the exposed coffee surface,
///    with `q` from `surface::open_top_flux` and `A_open` the lid opening's
///    share of the top (`MugParameters::open_top_area_m2`). Evaporation continues at and below
///    room temperature, so the coffee can settle below room.
///
/// 3. **Mug → Environment (convection)**: `Q = h × A_walls × (T_mug - T_room)`
///    Heat lost from the mug's outer walls
///
/// ## Numerical Method
///
/// Uses forward Euler integration with a 0.6-second time step:
/// ```text
/// dT/dt = -Q / (m × c)
/// T_new = T_old + dT × dt
/// ```
/// where `m × c` is the thermal mass (mass × specific heat capacity).
///
/// ## Assumptions and Limitations
///
/// - Uniform temperature distribution within coffee and mug (lumped capacitance)
/// - The lid opening (`exposure`) sets the open top area; the mug is always present
/// - Constant room temperature
/// - 50% relative humidity assumed for evaporative cooling
/// - Simplified geometry (cylindrical approximation)
///
/// Returns temperature data points, 5 per minute for smooth curves.
pub fn simulate_lumped_numerical(
    params: &MugParameters,
    duration_minutes: usize,
) -> Vec<CoolingDataPoint> {
    let record_interval = 0.2; // Record every 0.2 minutes (12 seconds)
    let mut results = Vec::with_capacity((duration_minutes as f64 / record_interval) as usize + 1);

    // Thermal masses (J/K)
    let coffee_thermal_mass = params.coffee_thermal_mass();
    let mug_thermal_mass = params.mug_thermal_mass();

    // Initial temperatures
    let mut t_coffee = params.coffee_temp_c;
    let t_room = params.room_temp_c;
    let mut t_mug = if params.preheated {
        params.preheat_temp_c
    } else {
        t_room
    };

    // Heat transfer coefficients (W/K, i.e., J/(s·K))
    //
    // Natural convection in air: ~5-15 W/(m²·K)
    // For a mug with ~0.03 m² surface area, that's ~0.15-0.45 W/K total
    //
    // Coffee-to-mug conductance from the wetted wall area, water-side film
    // and half the wall thickness (see MugParameters::coffee_mug_conductance)
    let coffee_mug_coefficient = params.coffee_mug_conductance();

    // Environment heat loss coefficients using actual surface areas; the open
    // top area carries every coffee-surface loss (evaporation, convection, radiation)
    let a_open = params.open_top_area_m2();
    let a_outer = params.outer_surface_area_m2();
    let mug_env_coefficient = params.h_conv * a_outer;

    // Mug emissivity for radiation
    let mug_emissivity = params.emissivity;

    // Time step for numerical integration (in minutes)
    let dt = 0.01; // 0.6 seconds
    let dt_seconds = dt * 60.0;

    // Record initial state
    results.push(CoolingDataPoint {
        time_minutes: 0.0,
        coffee_temp_c: t_coffee,
        mug_temp_c: t_mug,
    });

    let mut time = 0.0;
    let mut next_record = record_interval;

    while time < duration_minutes as f64 {
        // Heat flow from coffee to mug (W = J/s)
        // Newton's Law of Cooling: Q = k × ΔT
        let q_coffee_to_mug = coffee_mug_coefficient * (t_coffee - t_mug);

        // Heat flow from mug to environment (convection from outer walls)
        let q_mug_to_env = mug_env_coefficient * (t_mug - t_room);

        // Heat flow from the coffee through the lid opening: convection,
        // radiation and evaporation, evaluated exactly each step
        let q_coffee_to_env = a_open * open_top_flux(t_coffee, t_room, params.h_conv);

        // Radiative heat loss from the mug's outer walls
        let t_mug_k = t_mug + CELSIUS_TO_KELVIN;
        let t_room_k = t_room + CELSIUS_TO_KELVIN;
        let q_rad_mug =
            mug_emissivity * STEFAN_BOLTZMANN * a_outer * (t_mug_k.powi(4) - t_room_k.powi(4));

        // Temperature changes using forward Euler integration
        let dt_coffee = -(q_coffee_to_mug + q_coffee_to_env) * dt_seconds / coffee_thermal_mass;
        let dt_mug = (q_coffee_to_mug - q_mug_to_env - q_rad_mug) * dt_seconds / mug_thermal_mass;

        t_coffee += dt_coffee;
        t_mug += dt_mug;
        time += dt;

        // Record at each interval
        if time >= next_record {
            results.push(CoolingDataPoint {
                time_minutes: next_record,
                coffee_temp_c: t_coffee,
                mug_temp_c: t_mug,
            });
            next_record += record_interval;
        }
    }

    results
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::MugMaterial;

    #[test]
    fn test_preheated_warmer_than_cold() {
        let base_params = MugParameters {
            material: MugMaterial::Ceramic,
            volume_ml: 350.0,
            wall_thickness_mm: 5.0,
            coffee_temp_c: 95.0,
            room_temp_c: 22.0,
            preheated: false,
            preheat_temp_c: 90.0,
            h_conv: 10.0,
            emissivity: MugMaterial::Ceramic.emissivity(),
            exposure: 0.2,
        };

        let cold_data = simulate_lumped_numerical(&base_params, 30);

        let hot_params = MugParameters {
            preheated: true,
            ..base_params
        };
        let hot_data = simulate_lumped_numerical(&hot_params, 30);

        // Preheated should be warmer at every point
        for (cold, hot) in cold_data.iter().zip(hot_data.iter()) {
            assert!(
                hot.coffee_temp_c >= cold.coffee_temp_c,
                "At t={}, hot ({:.1}) should be >= cold ({:.1})",
                cold.time_minutes,
                hot.coffee_temp_c,
                cold.coffee_temp_c
            );
        }
    }

    #[test]
    fn test_cooling_towards_room_temp() {
        let params = MugParameters {
            material: MugMaterial::Ceramic,
            volume_ml: 350.0,
            wall_thickness_mm: 5.0,
            coffee_temp_c: 95.0,
            room_temp_c: 22.0,
            preheated: false,
            preheat_temp_c: 90.0,
            h_conv: 10.0,
            emissivity: MugMaterial::Ceramic.emissivity(),
            exposure: 0.2,
        };

        let data = simulate_lumped_numerical(&params, 60);

        // Temperature should always decrease
        for i in 1..data.len() {
            assert!(
                data[i].coffee_temp_c <= data[i - 1].coffee_temp_c,
                "Temperature should not increase over time"
            );
        }

        // Should approach room temp (with realistic convection, coffee is still warm after 60 min)
        let final_temp = data.last().unwrap().coffee_temp_c;
        assert!(
            final_temp < 60.0,
            "After 60 min, temp ({:.1}) should be cooling towards room temp",
            final_temp
        );
    }

    #[test]
    fn test_all_materials_cool_properly() {
        let materials = [
            MugMaterial::Ceramic,
            MugMaterial::Glass,
            MugMaterial::StainlessSteel,
            MugMaterial::Plastic,
        ];

        for material in materials {
            let params = MugParameters {
                material,
                volume_ml: 350.0,
                wall_thickness_mm: 5.0,
                coffee_temp_c: 95.0,
                room_temp_c: 22.0,
                preheated: false,
                preheat_temp_c: 90.0,
                h_conv: 10.0,
                emissivity: material.emissivity(),
                exposure: 0.2,
            };

            let data = simulate_lumped_numerical(&params, 30);

            // Temperature should decrease
            assert!(
                data.last().unwrap().coffee_temp_c < data.first().unwrap().coffee_temp_c,
                "{:?} mug should cool over time",
                material
            );

            // Should stay above room temp
            assert!(
                data.last().unwrap().coffee_temp_c > params.room_temp_c,
                "{:?} mug should stay above room temp after 30 min",
                material
            );
        }
    }

    #[test]
    fn test_stainless_steel_retains_heat_better() {
        // Stainless steel has higher thermal mass, should retain heat better
        let base = MugParameters {
            material: MugMaterial::Ceramic,
            volume_ml: 350.0,
            wall_thickness_mm: 5.0,
            coffee_temp_c: 95.0,
            room_temp_c: 22.0,
            preheated: false,
            preheat_temp_c: 90.0,
            h_conv: 10.0,
            emissivity: MugMaterial::Ceramic.emissivity(),
            exposure: 0.2,
        };

        let ceramic_data = simulate_lumped_numerical(&base, 60);
        let steel_data = simulate_lumped_numerical(
            &MugParameters {
                material: MugMaterial::StainlessSteel,
                ..base
            },
            60,
        );

        // After 60 minutes, stainless steel should be warmer (or similar)
        // due to higher thermal mass
        let ceramic_final = ceramic_data.last().unwrap().coffee_temp_c;
        let steel_final = steel_data.last().unwrap().coffee_temp_c;

        // Steel should retain at least as much heat
        assert!(
            steel_final >= ceramic_final - 5.0,
            "Steel ({:.1}°C) should retain heat at least as well as ceramic ({:.1}°C)",
            steel_final,
            ceramic_final
        );
    }

    #[test]
    fn test_extreme_cold_room() {
        let params = MugParameters {
            material: MugMaterial::Ceramic,
            volume_ml: 350.0,
            wall_thickness_mm: 5.0,
            coffee_temp_c: 95.0,
            room_temp_c: -10.0, // Very cold room
            preheated: false,
            preheat_temp_c: 90.0,
            h_conv: 10.0,
            emissivity: MugMaterial::Ceramic.emissivity(),
            exposure: 0.2,
        };

        let data = simulate_lumped_numerical(&params, 60);

        // Should cool faster in cold room
        let final_temp = data.last().unwrap().coffee_temp_c;
        assert!(
            final_temp > params.room_temp_c,
            "Coffee should stay above room temp"
        );
        assert!(
            final_temp < 50.0,
            "Coffee should cool significantly in cold room, got {:.1}°C",
            final_temp
        );
    }

    #[test]
    fn test_hot_room() {
        let params = MugParameters {
            material: MugMaterial::Ceramic,
            volume_ml: 350.0,
            wall_thickness_mm: 5.0,
            coffee_temp_c: 95.0,
            room_temp_c: 40.0, // Hot room
            preheated: false,
            preheat_temp_c: 90.0,
            h_conv: 10.0,
            emissivity: MugMaterial::Ceramic.emissivity(),
            exposure: 0.2,
        };

        let data = simulate_lumped_numerical(&params, 60);

        // Should cool slower in hot room. With evaporation, coffee can cool
        // slightly below room temp (wet-bulb effect), but not drastically.
        let final_temp = data.last().unwrap().coffee_temp_c;
        assert!(
            final_temp > params.room_temp_c - 10.0,
            "Coffee should stay near room temp, got {:.1}°C",
            final_temp
        );
        assert!(
            final_temp > 30.0,
            "Coffee should stay warmer in hot room, got {:.1}°C",
            final_temp
        );
    }

    #[test]
    fn test_small_vs_large_mug() {
        let small_params = MugParameters {
            material: MugMaterial::Ceramic,
            volume_ml: 100.0, // Espresso cup
            wall_thickness_mm: 3.0,
            coffee_temp_c: 95.0,
            room_temp_c: 22.0,
            preheated: false,
            preheat_temp_c: 90.0,
            h_conv: 10.0,
            emissivity: MugMaterial::Ceramic.emissivity(),
            exposure: 0.2,
        };

        let large_params = MugParameters {
            volume_ml: 500.0, // Large mug
            wall_thickness_mm: 6.0,
            ..small_params
        };

        let small_data = simulate_lumped_numerical(&small_params, 30);
        let large_data = simulate_lumped_numerical(&large_params, 30);

        // Small cup should cool faster (less thermal mass in coffee)
        let small_final = small_data.last().unwrap().coffee_temp_c;
        let large_final = large_data.last().unwrap().coffee_temp_c;

        assert!(
            large_final > small_final,
            "Large mug ({:.1}°C) should stay warmer than small cup ({:.1}°C)",
            large_final,
            small_final
        );
    }

    #[test]
    fn test_crossing_time_accuracy() {
        let params = MugParameters {
            material: MugMaterial::Ceramic,
            volume_ml: 350.0,
            wall_thickness_mm: 5.0,
            coffee_temp_c: 95.0,
            room_temp_c: 22.0,
            preheated: false,
            preheat_temp_c: 90.0,
            h_conv: 10.0,
            emissivity: MugMaterial::Ceramic.emissivity(),
            exposure: 0.2,
        };

        let data = simulate_lumped_numerical(&params, 120);
        let threshold = 60.0;

        // Find crossing time using linear interpolation (same algorithm as UI)
        let mut crossing_time: Option<f64> = None;
        for i in 1..data.len() {
            let prev = &data[i - 1];
            let curr = &data[i];
            if prev.coffee_temp_c >= threshold && curr.coffee_temp_c < threshold {
                let t =
                    (threshold - prev.coffee_temp_c) / (curr.coffee_temp_c - prev.coffee_temp_c);
                crossing_time = Some(prev.time_minutes + t);
                break;
            }
        }

        // Should find a crossing
        assert!(crossing_time.is_some(), "Should find threshold crossing");

        let cross_t = crossing_time.unwrap();

        // Crossing should be reasonable (between 10-60 minutes for these params)
        assert!(
            cross_t > 10.0 && cross_t < 60.0,
            "Crossing time ({:.1} min) should be reasonable",
            cross_t
        );

        // Verify by checking temperature around crossing time
        let before_idx = data
            .iter()
            .position(|p| p.time_minutes > cross_t - 1.0)
            .unwrap();
        let after_idx = data
            .iter()
            .position(|p| p.time_minutes > cross_t + 1.0)
            .unwrap();

        assert!(
            data[before_idx].coffee_temp_c >= threshold - 2.0,
            "Temperature before crossing should be near threshold"
        );
        assert!(
            data[after_idx].coffee_temp_c <= threshold + 2.0,
            "Temperature after crossing should be near threshold"
        );
    }

    #[test]
    fn test_mug_temp_equilibrium() {
        // Mug temperature should approach coffee temperature initially,
        // then both should cool together
        let params = MugParameters {
            material: MugMaterial::Ceramic,
            volume_ml: 350.0,
            wall_thickness_mm: 5.0,
            coffee_temp_c: 95.0,
            room_temp_c: 22.0,
            preheated: false,
            preheat_temp_c: 90.0,
            h_conv: 10.0,
            emissivity: MugMaterial::Ceramic.emissivity(),
            exposure: 0.2,
        };

        let data = simulate_lumped_numerical(&params, 30);

        // Initially, mug is at room temp
        assert!(
            (data[0].mug_temp_c - params.room_temp_c).abs() < 1.0,
            "Cold mug should start at room temp"
        );

        // Mug should heat up initially
        assert!(
            data[10].mug_temp_c > data[0].mug_temp_c,
            "Mug should heat up from coffee"
        );

        // Eventually, mug and coffee temps should be close
        let last = data.last().unwrap();
        assert!(
            (last.coffee_temp_c - last.mug_temp_c).abs() < 10.0,
            "Coffee ({:.1}°C) and mug ({:.1}°C) should approach equilibrium",
            last.coffee_temp_c,
            last.mug_temp_c
        );
    }

    #[test]
    fn test_preheated_mug_starts_hot() {
        let params = MugParameters {
            material: MugMaterial::Ceramic,
            volume_ml: 350.0,
            wall_thickness_mm: 5.0,
            coffee_temp_c: 95.0,
            room_temp_c: 22.0,
            preheated: true,
            preheat_temp_c: 90.0,
            h_conv: 10.0,
            emissivity: MugMaterial::Ceramic.emissivity(),
            exposure: 0.2,
        };

        let data = simulate_lumped_numerical(&params, 30);

        // Preheated mug should start at ~90°C
        assert!(
            data[0].mug_temp_c > 85.0,
            "Preheated mug should start hot, got {:.1}°C",
            data[0].mug_temp_c
        );
    }

    #[test]
    fn test_data_point_timing() {
        let params = MugParameters {
            material: MugMaterial::Ceramic,
            volume_ml: 350.0,
            wall_thickness_mm: 5.0,
            coffee_temp_c: 95.0,
            room_temp_c: 22.0,
            preheated: false,
            preheat_temp_c: 90.0,
            h_conv: 10.0,
            emissivity: MugMaterial::Ceramic.emissivity(),
            exposure: 0.2,
        };

        let data = simulate_lumped_numerical(&params, 10);

        // Should have ~5 points per minute (every 0.2 min)
        // For 10 minutes: ~50 points + initial = ~51
        assert!(
            data.len() >= 50 && data.len() <= 55,
            "Should have ~5 data points per minute, got {}",
            data.len()
        );

        // First point should be at t=0
        assert!(
            data[0].time_minutes.abs() < 0.001,
            "First point should be at t=0"
        );

        // Points should be evenly spaced (~0.2 min apart)
        for i in 1..data.len() {
            let dt = data[i].time_minutes - data[i - 1].time_minutes;
            assert!(
                (dt - 0.2).abs() < 0.01,
                "Points should be 0.2 min apart, got {:.3}",
                dt
            );
        }
    }

    #[test]
    fn test_evaporation_increases_cooling() {
        // With evaporation (hot coffee in dry air), the model should cool faster
        // than a hypothetical convection-only model. We verify by checking that
        // 95°C coffee reaches 60°C in under 30 minutes (convection-only took ~35+).
        let params = MugParameters {
            material: MugMaterial::Ceramic,
            volume_ml: 350.0,
            wall_thickness_mm: 5.0,
            coffee_temp_c: 95.0,
            room_temp_c: 22.0,
            preheated: false,
            preheat_temp_c: 90.0,
            h_conv: 10.0,
            emissivity: MugMaterial::Ceramic.emissivity(),
            exposure: 0.2,
        };

        let data = simulate_lumped_numerical(&params, 60);

        // Find when coffee crosses 60°C
        let mut crossing_time: Option<f64> = None;
        for i in 1..data.len() {
            if data[i - 1].coffee_temp_c >= 60.0 && data[i].coffee_temp_c < 60.0 {
                let frac = (60.0 - data[i - 1].coffee_temp_c)
                    / (data[i].coffee_temp_c - data[i - 1].coffee_temp_c);
                crossing_time = Some(data[i - 1].time_minutes + frac * 0.2);
                break;
            }
        }

        let cross = crossing_time.expect("Should cross 60°C");
        assert!(
            cross < 30.0,
            "With evaporation, 60°C crossing should be under 30 min, got {:.1}",
            cross
        );
    }

    #[test]
    fn test_radiation_order_of_magnitude() {
        // Radiative heat loss from a mug-sized surface at 80°C should be ~5–15W
        let t_hot_k: f64 = 80.0 + CELSIUS_TO_KELVIN;
        let t_room_k: f64 = 22.0 + CELSIUS_TO_KELVIN;
        let emissivity = 0.90; // ceramic
        let a_outer = 0.03; // ~0.03 m² for a typical mug

        let q_rad = emissivity * STEFAN_BOLTZMANN * a_outer * (t_hot_k.powi(4) - t_room_k.powi(4));

        assert!(
            q_rad > 2.0 && q_rad < 20.0,
            "Radiative loss should be ~5-15W, got {:.2}W",
            q_rad
        );
    }

    #[test]
    fn lumped_settles_below_room() {
        // Coffee keeps evaporating at room temperature, so with the mug warmed
        // by the room it settles a few degrees below room, toward wet-bulb.
        let params = MugParameters {
            material: MugMaterial::Ceramic,
            volume_ml: 350.0,
            wall_thickness_mm: 5.0,
            coffee_temp_c: 30.0,
            room_temp_c: 22.0,
            preheated: false,
            preheat_temp_c: 22.0,
            h_conv: 10.0,
            emissivity: MugMaterial::Ceramic.emissivity(),
            exposure: 1.0,
        };
        let numerical = simulate_lumped_numerical(&params, 600);
        let analytical = crate::physics::simulate_lumped_analytical(&params, 600);
        for (name, data) in [("numerical", &numerical), ("analytical", &analytical)] {
            let final_coffee = data.last().unwrap().coffee_temp_c;
            assert!(
                final_coffee < 22.0 && final_coffee > 10.0,
                "{name}: coffee after 600 min = {final_coffee:.3} °C, expected in (10, 22)"
            );
        }
    }

    #[test]
    fn lumped_numerical_steel_thin_wall_monotone() {
        // Thin steel has a tiny mug heat capacity and a large conductance:
        // the explicit step must not make the mug temperature oscillate.
        let params = MugParameters {
            material: MugMaterial::StainlessSteel,
            volume_ml: 350.0,
            wall_thickness_mm: 1.0,
            coffee_temp_c: 95.0,
            room_temp_c: 22.0,
            preheated: false,
            preheat_temp_c: 22.0,
            h_conv: 10.0,
            emissivity: MugMaterial::StainlessSteel.emissivity(),
            exposure: 0.2,
        };
        let data = simulate_lumped_numerical(&params, 5);
        let mug: Vec<f64> = data.iter().map(|p| p.mug_temp_c).collect();
        let diffs: Vec<f64> = mug.windows(2).map(|w| w[1] - w[0]).collect();
        assert!(diffs[0] > 0.0, "mug should warm first");
        let sign_changes = diffs.windows(2).filter(|w| w[0] * w[1] < 0.0).count();
        assert!(
            sign_changes <= 1,
            "mug temp oscillates: {sign_changes} sign changes"
        );
    }
}

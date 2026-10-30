//! The numbers the home page states as its conclusion ("The Science" and
//! "What it adds up to" in `components/app.rs`), checked against both
//! numerical models. When physics changes move a result outside its band,
//! update the page text and the band together.

use super::*;

const COFFEE_C: f64 = 90.0;
const PREHEAT_C: f64 = 80.0;
const ROOM_C: f64 = 22.0;
const DRINKABLE_C: f64 = 60.0;
const H_CONV: f64 = 20.0;
const LID_OPENING: f64 = 0.2;
const RUN_MINUTES: usize = 30;

struct Outcome {
    /// How much warmer the preheated coffee is on average while the cold one is drinkable (°C).
    penalty_c: f64,
    cold_crossing_min: f64,
    extra_min: f64,
}

fn mug(volume_ml: f64, wall_mm: f64, preheated: bool) -> MugParameters {
    let material = MugMaterial::Ceramic;
    MugParameters {
        material,
        volume_ml,
        wall_thickness_mm: wall_mm,
        coffee_temp_c: COFFEE_C,
        room_temp_c: ROOM_C,
        preheated,
        preheat_temp_c: PREHEAT_C,
        h_conv: H_CONV,
        emissivity: material.emissivity(),
        exposure: LID_OPENING,
    }
}

fn outcome(cold: &[(f64, f64)], hot: &[(f64, f64)]) -> Outcome {
    let m = compute_comparison(cold, hot, DRINKABLE_C);
    let cold_t = m
        .crossing_time_cold
        .expect("cold coffee crosses within the run");
    let hot_t = m
        .crossing_time_hot
        .expect("preheated coffee crosses within the run");
    Outcome {
        penalty_c: m.temp_diff,
        cold_crossing_min: cold_t,
        extra_min: hot_t - cold_t,
    }
}

fn lumped(volume_ml: f64, wall_mm: f64) -> Outcome {
    let run = |pre| {
        CoolingDataPoint::coffee_series(&simulate_lumped_numerical(
            &mug(volume_ml, wall_mm, pre),
            RUN_MINUTES,
        ))
    };
    outcome(&run(false), &run(true))
}

fn heat_equation(volume_ml: f64, wall_mm: f64) -> Outcome {
    // Same construction as the heat-equation panel's numerical path.
    let inner_cm = (3.0 * volume_ml / (4.0 * std::f64::consts::PI)).cbrt();
    let outer_cm = inner_cm + wall_mm / 10.0;
    let material = MugMaterial::Ceramic;
    let run = |shell_c| {
        let params = SphereHeatParams {
            outer_radius_m: outer_cm / 100.0,
            inner_radius_m: inner_cm / 100.0,
            shell_temp_c: shell_c,
            liquid_temp_c: COFFEE_C,
            room_temp_c: ROOM_C,
            shell_material: HeatMaterial::from_mug(material, material.emissivity()),
            liquid_material: HeatMaterial::water(),
            medium_material: HeatMaterial::air(),
            exposure: LID_OPENING,
            convection_coeff: H_CONV,
            r_max_multiplier: 10.0,
            outer_boundary: OuterBoundary::Ambient,
        };
        let snaps = simulate_sphere_heat_numerical(&params, RUN_MINUTES as f64 * 60.0, 200);
        seconds_to_minutes(&liquid_average_temperature(&snaps))
    };
    outcome(&run(ROOM_C), &run(PREHEAT_C))
}

fn heat_equation_analytical(volume_ml: f64, wall_mm: f64) -> Outcome {
    // Same construction as the heat-equation panel's analytical path.
    let r_inner = (3.0 * volume_ml / (4.0 * std::f64::consts::PI)).cbrt() / 100.0;
    let r_outer = r_inner + wall_mm / 1000.0;
    let material = MugMaterial::Ceramic;
    let run = |shell_c| {
        let mug = ThreeRegionMug {
            r_inner,
            r_outer,
            r_max: r_outer * 10.0,
            shell: HeatMaterial::from_mug(material, material.emissivity()),
            liquid_temp_c: COFFEE_C,
            shell_temp_c: shell_c,
            room_temp_c: ROOM_C,
            h_conv: H_CONV,
            exposure: LID_OPENING,
            num_modes: 50,
        };
        let snaps = simulate_three_region_relinearised(
            &mug,
            ANALYTICAL_SEGMENTS,
            RUN_MINUTES as f64 * 60.0,
            100,
            80,
            0.10,
        );
        seconds_to_minutes(&liquid_average_temperature(&convert_snapshots(snaps)))
    };
    outcome(&run(ROOM_C), &run(PREHEAT_C))
}

fn both(volume_ml: f64, wall_mm: f64) -> [(&'static str, Outcome); 2] {
    [
        ("lumped", lumped(volume_ml, wall_mm)),
        ("heat equation", heat_equation(volume_ml, wall_mm)),
    ]
}

/// Theory page: the analytical heat-equation solver agrees with the
/// numerical one to within 0.1 min of extra drinkable time and 0.2 °C of
/// penalty; its crossings land early (0.2–0.4 min, larger mugs more),
/// because it linearises the open-top loss at the liquid average.
#[test]
fn analytical_heat_equation_agrees_with_numerical() {
    for (volume, wall) in [(180.0, 5.0), (500.0, 5.0), (180.0, 2.0), (180.0, 8.0)] {
        let a = heat_equation_analytical(volume, wall);
        let n = heat_equation(volume, wall);
        let case = format!("{volume} ml {wall} mm");
        let extra = a.extra_min - n.extra_min;
        assert!(extra.abs() < 0.1, "{case}: extra differs by {extra:.3} min");
        let penalty = a.penalty_c - n.penalty_c;
        assert!(
            penalty.abs() < 0.2,
            "{case}: penalty differs by {penalty:.3} °C"
        );
        let cold = a.cold_crossing_min - n.cold_crossing_min;
        let hot = cold + extra;
        for (which, d) in [("cold", cold), ("preheated", hot)] {
            assert!(
                (-0.4..0.0).contains(&d),
                "{case}: {which} crossing {d:.3} min from numerical"
            );
        }
    }
}

/// "The Science": a cold 180 ml ceramic mug takes about 13 °C out of the
/// coffee as the two even out; a preheated one about 2 °C.
#[test]
fn equalisation_drop_matches_intro() {
    let cold = mug(180.0, 5.0, false);
    let (c_coffee, c_mug) = (cold.coffee_thermal_mass(), cold.mug_thermal_mass());
    let drop = |mug_c: f64| (COFFEE_C - mug_c) * c_mug / (c_coffee + c_mug);
    assert!(
        (12.5..13.5).contains(&drop(ROOM_C)),
        "cold drop {:.2}",
        drop(ROOM_C)
    );
    assert!(
        (1.5..2.5).contains(&drop(PREHEAT_C)),
        "preheated drop {:.2}",
        drop(PREHEAT_C)
    );
}

/// "about 3½ more minutes … barely depends on mug size".
#[test]
fn extra_time_is_about_three_and_a_half_minutes_across_sizes() {
    for volume in [180.0, 250.0, 350.0, 500.0] {
        for (model, o) in both(volume, 5.0) {
            assert!(
                (3.0..4.0).contains(&o.extra_min),
                "{model} {volume} ml: extra {:.2} min",
                o.extra_min
            );
        }
    }
}

/// "stays 5–8 °C warmer on average while you drink it".
#[test]
fn penalty_is_five_to_eight_degrees_for_typical_mugs() {
    for volume in [180.0, 250.0, 350.0, 500.0] {
        for (model, o) in both(volume, 5.0) {
            assert!(
                (4.5..8.0).contains(&o.penalty_c),
                "{model} {volume} ml: penalty {:.2} °C",
                o.penalty_c
            );
        }
    }
}

/// "For a 180 ml cup that's about a third more time".
#[test]
fn extra_time_is_about_a_third_for_a_small_cup() {
    for (model, o) in both(180.0, 5.0) {
        let share = o.extra_min / o.cold_crossing_min;
        assert!(
            (0.3..0.45).contains(&share),
            "{model}: {:.0}% more",
            share * 100.0
        );
    }
}

/// "A thick-walled mug (8 mm) gains about 6 minutes; a thin one (2 mm) only about one."
#[test]
fn wall_thickness_sets_how_much_preheating_helps() {
    for (model, o) in both(180.0, 8.0) {
        assert!(
            (5.5..7.0).contains(&o.extra_min),
            "{model} 8 mm: extra {:.2} min",
            o.extra_min
        );
    }
    for (model, o) in both(180.0, 2.0) {
        assert!(
            (0.8..1.7).contains(&o.extra_min),
            "{model} 2 mm: extra {:.2} min",
            o.extra_min
        );
    }
}

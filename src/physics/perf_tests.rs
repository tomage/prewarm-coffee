//! Performance regression tests for physics solvers.
//!
//! Measures execution time of the four main solver paths and compares
//! against stored baselines.  Run with `UPDATE_BASELINES=1` to regenerate
//! the baselines file (`tests/perf-baselines-unit.json`).

use super::*;
use std::collections::{BTreeMap, HashMap};
use std::time::Instant;

const BASELINES_PATH: &str = "tests/perf-baselines-unit.json";
const TOLERANCE: f64 = 0.10; // 10%
const ABS_TOLERANCE_MS: f64 = 20.0; // minimum absolute tolerance for fast benchmarks
const SAMPLE_COUNT_BASELINE: usize = 10;
const SAMPLE_COUNT_CHECK: usize = 5;
const WARMUP_COUNT: usize = 2;

fn median(values: &[f64]) -> f64 {
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mid = sorted.len() / 2;
    if sorted.len() % 2 == 0 {
        (sorted[mid - 1] + sorted[mid]) / 2.0
    } else {
        sorted[mid]
    }
}

fn load_baselines() -> HashMap<String, f64> {
    let Ok(raw) = std::fs::read_to_string(BASELINES_PATH) else {
        return HashMap::new();
    };
    let Ok(map) = serde_json::from_str::<HashMap<String, serde_json::Value>>(&raw) else {
        return HashMap::new();
    };
    map.into_iter()
        .filter_map(|(k, v)| v.get("median_ms")?.as_f64().map(|ms| (k, ms)))
        .collect()
}

fn save_baselines(data: &HashMap<String, f64>) {
    let map: BTreeMap<&str, serde_json::Value> = data
        .iter()
        .map(|(k, &v)| {
            let mut inner = serde_json::Map::new();
            inner.insert("median_ms".to_string(), serde_json::json!((v * 100.0).round() / 100.0));
            (k.as_str(), serde_json::Value::Object(inner))
        })
        .collect();
    let json = serde_json::to_string_pretty(&map).unwrap();
    std::fs::write(BASELINES_PATH, format!("{json}\n")).unwrap();
}

/// Run a benchmark: warmup, then measure runs, return median ms.
fn bench<F: FnMut()>(mut f: F) -> f64 {
    let updating = std::env::var("UPDATE_BASELINES").is_ok();
    let samples = if updating { SAMPLE_COUNT_BASELINE } else { SAMPLE_COUNT_CHECK };
    for _ in 0..WARMUP_COUNT {
        f();
    }
    let times: Vec<f64> = (0..samples)
        .map(|_| {
            let start = Instant::now();
            f();
            start.elapsed().as_secs_f64() * 1000.0
        })
        .collect();
    median(&times)
}

fn assert_perf(name: &str, measured_ms: f64) {
    let updating = std::env::var("UPDATE_BASELINES").is_ok();

    if updating {
        let mut baselines = load_baselines();
        baselines.insert(name.to_string(), measured_ms);
        save_baselines(&baselines);
        eprintln!("  {name} — measured: {:.1}ms (baseline updated)", measured_ms);
        return;
    }

    let baselines = load_baselines();
    let Some(&baseline) = baselines.get(name) else {
        eprintln!("  ⚠ No baseline for {name} — skipping assertion (run UPDATE_BASELINES=1 cargo test to create)");
        return;
    };

    let threshold = f64::max(baseline * (1.0 + TOLERANCE), baseline + ABS_TOLERANCE_MS);
    let delta = (measured_ms - baseline) / baseline * 100.0;
    let status = if measured_ms <= threshold { "✓" } else { "✗" };
    eprintln!(
        "  {name} — baseline: {:.1}ms, measured: {:.1}ms ({:+.1}%) {status}",
        baseline, measured_ms, delta
    );

    assert!(
        measured_ms <= threshold,
        "{name} regressed: {:.1}ms vs baseline {:.1}ms ({:+.1}%, threshold {:.1}ms)",
        measured_ms, baseline, delta, threshold
    );
}

// ─── Lumped capacitance ─────────────────────────────────────────────

fn lumped_params() -> MugParameters {
    MugParameters {
        material: MugMaterial::Ceramic,
        volume_ml: 350.0,
        wall_thickness_mm: 5.0,
        coffee_temp_c: 95.0,
        room_temp_c: 22.0,
        preheated: false,
        preheat_temp_c: 90.0,
        h_conv: 10.0,
        emissivity: MugMaterial::Ceramic.emissivity(),
        exposure_pct: 0.2,
    }
}

#[test]
fn perf_lumped_numerical() {
    let params = lumped_params();
    let med = bench(|| {
        let _ = simulate_lumped_numerical(&params, 60);
    });
    assert_perf("lumped_numerical", med);
}

#[test]
fn perf_lumped_analytical() {
    let params = lumped_params();
    let med = bench(|| {
        let _ = simulate_lumped_analytical(&params, 60);
    });
    assert_perf("lumped_analytical", med);
}

// ─── Heat equation (3-region sphere) ────────────────────────────────

fn heat_eq_params() -> SphereHeatParams {
    let inner_r_cm = (3.0 * 180.0_f64 / (4.0 * std::f64::consts::PI)).cbrt();
    let outer_r_cm = inner_r_cm + 5.0 / 10.0;
    let chart_max_r_cm = 10.0_f64;
    SphereHeatParams {
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
    }
}

#[test]
fn perf_heat_equation_numerical() {
    let params = heat_eq_params();
    let med = bench(|| {
        let _ = simulate_sphere_heat_numerical(&params, 3600.0, 100);
    });
    assert_perf("heat_equation_numerical", med);
}

#[test]
fn perf_heat_equation_analytical() {
    let three_params = ThreeRegionParams::default();
    let med = bench(|| {
        let modes = compute_eigenmodes(&three_params);
        let _ = simulate_three_region_analytical(
            &three_params,
            &modes,
            3600.0,
            100,
            200,
            three_params.r_max,
        );
    });
    assert_perf("heat_equation_analytical", med);
}

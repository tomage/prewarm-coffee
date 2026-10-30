//! Agreement report: all four solvers over a grid of mugs and conditions.
//!
//! Prints, for each solver pair, the bias and spread of extra drinkable
//! minutes, penalty and crossing times, plus the worst scenarios; with
//! `GRID_CSV` set it also writes every scenario's results there. A report,
//! not a gate: it asserts nothing, so it is `#[ignore]`d and runs with
//! `make grid`.

use super::*;
use std::fmt::Write as _;

const COFFEE_C: f64 = 90.0;
const PREHEAT_C: f64 = 80.0;
const ROOM_C: f64 = 22.0;
const DRINKABLE_C: f64 = 60.0;
const RUN_MINUTES: usize = 60;

#[derive(Clone, Copy, Debug)]
struct Scenario {
    material: MugMaterial,
    volume_ml: f64,
    wall_mm: f64,
    exposure: f64,
    h_conv: f64,
}

#[derive(Clone, Copy, Debug, Default)]
struct Out {
    extra: Option<f64>,
    penalty: f64,
    cold: Option<f64>,
    hot: Option<f64>,
}

fn out(cold: &[(f64, f64)], hot: &[(f64, f64)]) -> Out {
    let m = compute_comparison(cold, hot, DRINKABLE_C);
    Out {
        extra: m
            .crossing_time_cold
            .zip(m.crossing_time_hot)
            .map(|(c, h)| h - c),
        penalty: m.temp_diff,
        cold: m.crossing_time_cold,
        hot: m.crossing_time_hot,
    }
}

fn mug(s: &Scenario, preheated: bool) -> MugParameters {
    MugParameters {
        material: s.material,
        volume_ml: s.volume_ml,
        wall_thickness_mm: s.wall_mm,
        coffee_temp_c: COFFEE_C,
        room_temp_c: ROOM_C,
        preheated,
        preheat_temp_c: PREHEAT_C,
        h_conv: s.h_conv,
        emissivity: s.material.emissivity(),
        exposure: s.exposure,
    }
}

fn lumped_numerical(s: &Scenario) -> Out {
    let run =
        |p| CoolingDataPoint::coffee_series(&simulate_lumped_numerical(&mug(s, p), RUN_MINUTES));
    out(&run(false), &run(true))
}

fn lumped_analytical(s: &Scenario) -> Out {
    let run =
        |p| CoolingDataPoint::coffee_series(&simulate_lumped_analytical(&mug(s, p), RUN_MINUTES));
    out(&run(false), &run(true))
}

fn heat_numerical(s: &Scenario) -> Out {
    let inner_cm = (3.0 * s.volume_ml / (4.0 * std::f64::consts::PI)).cbrt();
    let outer_cm = inner_cm + s.wall_mm / 10.0;
    let run = |shell_c| {
        let params = SphereHeatParams {
            outer_radius_m: outer_cm / 100.0,
            inner_radius_m: inner_cm / 100.0,
            shell_temp_c: shell_c,
            liquid_temp_c: COFFEE_C,
            room_temp_c: ROOM_C,
            shell_material: HeatMaterial::from_mug(s.material, s.material.emissivity()),
            liquid_material: HeatMaterial::water(),
            medium_material: HeatMaterial::air(),
            exposure: s.exposure,
            convection_coeff: s.h_conv,
            r_max_multiplier: 10.0,
            outer_boundary: OuterBoundary::Ambient,
        };
        let snaps = simulate_sphere_heat_numerical(&params, RUN_MINUTES as f64 * 60.0, 400);
        seconds_to_minutes(&liquid_average_temperature(&snaps))
    };
    out(&run(ROOM_C), &run(PREHEAT_C))
}

fn heat_analytical(s: &Scenario) -> Out {
    let r_inner = (3.0 * s.volume_ml / (4.0 * std::f64::consts::PI)).cbrt() / 100.0;
    let r_outer = r_inner + s.wall_mm / 1000.0;
    let run = |shell_c| {
        let m = ThreeRegionMug {
            r_inner,
            r_outer,
            r_max: r_outer * 10.0,
            shell: HeatMaterial::from_mug(s.material, s.material.emissivity()),
            liquid_temp_c: COFFEE_C,
            shell_temp_c: shell_c,
            room_temp_c: ROOM_C,
            h_conv: s.h_conv,
            exposure: s.exposure,
            num_modes: 50,
        };
        let snaps = simulate_three_region_relinearised(
            &m,
            ANALYTICAL_SEGMENTS,
            RUN_MINUTES as f64 * 60.0,
            200,
            80,
            0.10,
        );
        seconds_to_minutes(&liquid_average_temperature(&convert_snapshots(snaps)))
    };
    out(&run(ROOM_C), &run(PREHEAT_C))
}

fn scenarios() -> Vec<Scenario> {
    let mut v = Vec::new();
    for material in [
        MugMaterial::Ceramic,
        MugMaterial::Glass,
        MugMaterial::StainlessSteel,
        MugMaterial::Plastic,
    ] {
        for volume_ml in [150.0, 250.0, 350.0, 500.0] {
            for wall_mm in [2.0, 5.0, 8.0] {
                for exposure in [0.0, 0.2, 0.5, 1.0] {
                    for h_conv in [10.0, 20.0, 40.0] {
                        v.push(Scenario {
                            material,
                            volume_ml,
                            wall_mm,
                            exposure,
                            h_conv,
                        });
                    }
                }
            }
        }
    }
    v
}

fn f(x: Option<f64>) -> String {
    x.map_or("".into(), |x| format!("{x:.3}"))
}

struct Stat {
    name: &'static str,
    diffs: Vec<(f64, usize)>,
    missing: usize,
}

impl Stat {
    fn new(name: &'static str) -> Self {
        Self {
            name,
            diffs: Vec::new(),
            missing: 0,
        }
    }
    fn add(&mut self, a: Option<f64>, b: Option<f64>, i: usize) {
        match (a, b) {
            (Some(a), Some(b)) => self.diffs.push((a - b, i)),
            (None, None) => {}
            _ => self.missing += 1,
        }
    }
    fn line(&self) -> String {
        let n = self.diffs.len().max(1) as f64;
        let mean = self.diffs.iter().map(|d| d.0).sum::<f64>() / n;
        let mae = self.diffs.iter().map(|d| d.0.abs()).sum::<f64>() / n;
        let mut abs: Vec<f64> = self.diffs.iter().map(|d| d.0.abs()).collect();
        abs.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let p95 = abs
            .get(((abs.len() as f64) * 0.95) as usize)
            .copied()
            .unwrap_or(0.0);
        let max = abs.last().copied().unwrap_or(0.0);
        format!(
            "{:<34} n={:4} bias={:+.3} mae={:.3} p95={:.3} max={:.3} one-sided-missing={}",
            self.name,
            self.diffs.len(),
            mean,
            mae,
            p95,
            max,
            self.missing
        )
    }
}

#[test]
#[ignore = "report, not a gate; run with make grid"]
fn solver_grid_report() {
    let sc = scenarios();
    let threads = std::thread::available_parallelism().map_or(4, |n| n.get());
    let start = std::time::Instant::now();
    let mut results: Vec<Option<[Out; 4]>> = vec![None; sc.len()];
    std::thread::scope(|scope| {
        for (chunk_idx, chunk) in results.chunks_mut(sc.len().div_ceil(threads)).enumerate() {
            let sc = &sc;
            let base = chunk_idx * sc.len().div_ceil(threads);
            scope.spawn(move || {
                for (j, slot) in chunk.iter_mut().enumerate() {
                    let s = &sc[base + j];
                    *slot = Some([
                        lumped_analytical(s),
                        lumped_numerical(s),
                        heat_analytical(s),
                        heat_numerical(s),
                    ]);
                }
            });
        }
    });
    let results: Vec<[Out; 4]> = results.into_iter().map(Option::unwrap).collect();
    println!(
        "{} scenarios in {:.1}s on {threads} threads",
        sc.len(),
        start.elapsed().as_secs_f64()
    );

    let mut csv = String::from(
        "material,volume_ml,wall_mm,exposure,h_conv,\
la_extra,la_pen,la_cold,la_hot,ln_extra,ln_pen,ln_cold,ln_hot,\
ha_extra,ha_pen,ha_cold,ha_hot,hn_extra,hn_pen,hn_cold,hn_hot\n",
    );
    for (s, r) in sc.iter().zip(&results) {
        let _ = write!(
            csv,
            "{:?},{},{},{},{}",
            s.material, s.volume_ml, s.wall_mm, s.exposure, s.h_conv
        );
        for o in r {
            let _ = write!(
                csv,
                ",{},{:.3},{},{}",
                f(o.extra),
                o.penalty,
                f(o.cold),
                f(o.hot)
            );
        }
        csv.push('\n');
    }
    if let Ok(path) = std::env::var("GRID_CSV") {
        std::fs::write(&path, csv).unwrap();
        println!("csv -> {path}");
    }

    // (a, b, label): indices into [LA, LN, HA, HN]
    for (a, b, label) in [
        (0, 1, "lumped A−N"),
        (2, 3, "heat-eq A−N"),
        (1, 3, "lumped N − heat-eq N"),
    ] {
        let mut extra = Stat::new("extra min");
        let mut pen = Stat::new("penalty °C");
        let mut cold = Stat::new("cold crossing min");
        let mut hot = Stat::new("hot crossing min");
        for (i, r) in results.iter().enumerate() {
            extra.add(r[a].extra, r[b].extra, i);
            pen.add(Some(r[a].penalty), Some(r[b].penalty), i);
            cold.add(r[a].cold, r[b].cold, i);
            hot.add(r[a].hot, r[b].hot, i);
        }
        println!("\n== {label} ==");
        for st in [&extra, &pen, &cold, &hot] {
            println!("  {}", st.line());
        }
        let mut worst = extra.diffs.clone();
        worst.sort_by(|x, y| y.0.abs().partial_cmp(&x.0.abs()).unwrap());
        println!("  worst extra-min scenarios:");
        for (d, i) in worst.iter().take(8) {
            let s = &sc[*i];
            println!(
                "    {:+.3}  {:?} {}ml {}mm exp={} h={}  ({} vs {})",
                d,
                s.material,
                s.volume_ml,
                s.wall_mm,
                s.exposure,
                s.h_conv,
                f(results[*i][a].extra),
                f(results[*i][b].extra)
            );
        }
    }
}

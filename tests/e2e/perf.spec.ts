import { test, expect, type Page } from "@playwright/test";
import * as fs from "fs";
import * as path from "path";

// Run tests serially — perf measurements are timing-sensitive
test.describe.configure({ mode: "serial" });

const BASELINES_PATH = path.join(__dirname, "perf-baselines.json");
const TOLERANCE = 0.10; // 10%
const ABS_TOLERANCE_MS = 20; // minimum absolute tolerance to avoid flaky failures on fast pages
const SAMPLE_COUNT_BASELINE = 10;
const SAMPLE_COUNT_CHECK = 5;

// Routes with simulations set hasSimulation=true to wait for WASM computation
// to complete before stopping the timer. The .compute-time span shows "0.0 ms"
// initially and updates to a real value once the simulation finishes.
const routes: { path: string; name: string; hasSimulation?: boolean }[] = [
  { path: "/", name: "Home", hasSimulation: true },
  { path: "/about", name: "About" },
  { path: "/theory", name: "Theory" },
  { path: "/theory/lumped", name: "Theory: Lumped", hasSimulation: true },
  { path: "/theory/lumped?solver=numerical", name: "Theory: Lumped (Numerical)", hasSimulation: true },
  { path: "/theory/heat-equation", name: "Theory: Heat Equation", hasSimulation: true },
  { path: "/theory/heat-equation?solver=analytical", name: "Theory: Heat Equation (Analytical)", hasSimulation: true },
  { path: "/theory/heat-equation?solver=numerical", name: "Theory: Heat Equation (Numerical)", hasSimulation: true },
];

type Baseline = { median_ms: number; long_task_ms?: number };

function median(values: number[]): number {
  const sorted = [...values].sort((a, b) => a - b);
  const mid = Math.floor(sorted.length / 2);
  return sorted.length % 2 === 0
    ? (sorted[mid - 1] + sorted[mid]) / 2
    : sorted[mid];
}

function loadBaselines(): Record<string, Baseline> {
  try {
    const raw = fs.readFileSync(BASELINES_PATH, "utf-8");
    return JSON.parse(raw);
  } catch {
    return {};
  }
}

function saveBaselines(data: Record<string, Baseline>) {
  const sorted = Object.fromEntries(
    Object.entries(data).sort(([a], [b]) => a.localeCompare(b))
  );
  fs.writeFileSync(BASELINES_PATH, JSON.stringify(sorted, null, 2) + "\n");
}

const updating = !!process.env.UPDATE_BASELINES;
const baselines = loadBaselines();
const results: Record<string, number> = {};

for (const route of routes) {
  test(`perf: ${route.name} (${route.path})`, async ({ page }) => {
    // Wait for WASM simulation to complete (data-sim-ready attribute set by Leptos)
    const waitForSimulation = async () => {
      if (!route.hasSimulation) return;
      await page.locator(".solver-toggle[data-sim-ready]").waitFor();
    };

    // Warmup: uncounted load to compile WASM and cache route-specific resources
    await page.goto(route.path, { waitUntil: "load" });
    await waitForSimulation();

    const times: number[] = [];
    const samples = updating ? SAMPLE_COUNT_BASELINE : SAMPLE_COUNT_CHECK;
    for (let i = 0; i < samples; i++) {
      const start = Date.now();
      await page.goto(route.path, { waitUntil: "load" });
      await waitForSimulation();
      times.push(Date.now() - start);
    }

    const med = median(times);
    const rounded = Math.round(med);
    results[route.path] = rounded;

    if (updating) {
      // Read-merge-write to handle incremental updates
      const current = loadBaselines();
      current[route.path] = { median_ms: rounded };
      saveBaselines(current);
      console.log(`  ${route.path} — measured: ${rounded}ms (baseline updated)`);
      return;
    }

    const baseline = baselines[route.path];
    if (!baseline) {
      console.warn(`  ⚠ No baseline for ${route.path} — skipping assertion (run UPDATE_BASELINES=1 to create)`);
      return;
    }

    const threshold = Math.max(baseline.median_ms * (1 + TOLERANCE), baseline.median_ms + ABS_TOLERANCE_MS);
    const delta = ((med - baseline.median_ms) / baseline.median_ms) * 100;
    const status = med <= threshold ? "✓" : "✗";
    console.log(
      `  ${route.path} — baseline: ${baseline.median_ms}ms, measured: ${rounded}ms (${delta >= 0 ? "+" : ""}${delta.toFixed(1)}%) ${status}`
    );

    expect(
      med,
      `${route.name} (${route.path}) regressed: ${rounded}ms vs baseline ${baseline.median_ms}ms (+${delta.toFixed(1)}%, threshold ${Math.round(threshold)}ms)`
    ).toBeLessThanOrEqual(threshold);
  });
}

// ---------------------------------------------------------------------------
// Drag interaction — what a user feels while dragging a slider: from the first
// ArrowRight press until the chart's main line shows the final value's solve.
// Solve time dominates, so this is a before/after comparison on one machine,
// not an absolute framework number.
// ---------------------------------------------------------------------------

const SETTINGS_KEY = "prewarm_settings";
const DRAG_FROM_C = 70;
const DRAG_TO_C = 85;
// #coffee-temp has step="1" (mug_controls.rs), so 15 presses go 70 -> 85.
const DRAG_PRESSES = DRAG_TO_C - DRAG_FROM_C;
const SETTLE_TIMEOUT_MS = 10_000;

const interactionTabs: { tab: string; line: string }[] = [
  { tab: "lumped", line: ".desktop-sim polyline.temp-line.hot:not(.mug-temp)" },
  { tab: "heat_equation", line: ".desktop-sim polyline.temp-line.hot:not(.mug-temp)" },
  { tab: "compare", line: ".desktop-sim polyline.heateq-numerical" },
];

const interactionResults: Record<string, { median_ms: number; long_task_ms: number }> = {};

async function loadWithCoffeeTemp(page: Page, tab: string, tempC: number) {
  const seed = { version: 1, coffee_temp_c: tempC, use_fahrenheit: false };
  await page.goto(`/?tab=${tab}`, { waitUntil: "networkidle" });
  await page.evaluate(
    ([key, value]) => localStorage.setItem(key as string, JSON.stringify(value)),
    [SETTINGS_KEY, seed] as const
  );
  await page.goto(`/?tab=${tab}`, { waitUntil: "networkidle" });
  await expect(page.locator('.desktop-sim label[for="coffee-temp"]')).toContainText(String(tempC));
}

/** Resolves with the polyline's points once they are non-empty and unchanged
 * for STABLE_FRAMES consecutive frames: the label shows a new value as soon as
 * the setting changes, before the (possibly deferred) solve has redrawn the
 * chart, so a read right after the label check can capture a stale line. */
const STABLE_FRAMES = 5;
async function stableLine(page: Page, line: string): Promise<string> {
  const handle = await page.waitForFunction(
    ([sel, frames]) => {
      const w = window as unknown as { __stable?: { sel: string; last: string | null; n: number } };
      const pts = document.querySelector(sel as string)?.getAttribute("points") ?? "";
      if (!w.__stable || w.__stable.sel !== sel) w.__stable = { sel: sel as string, last: null, n: 0 };
      const st = w.__stable;
      st.n = pts !== "" && pts === st.last ? st.n + 1 : 0;
      st.last = pts;
      return st.n >= (frames as number) ? pts : false;
    },
    [line, STABLE_FRAMES] as const,
    { polling: "raf", timeout: SETTLE_TIMEOUT_MS }
  );
  return (await handle.jsonValue()) as string;
}

for (const { tab, line } of interactionTabs) {
  test(`perf: interaction: ${tab}`, async ({ page }) => {
    await page.setViewportSize({ width: 1280, height: 800 });

    // The resting line for the final value, from a fresh load.
    await loadWithCoffeeTemp(page, tab, DRAG_TO_C);
    const expected = await stableLine(page, line);

    const sample = async () => {
      await loadWithCoffeeTemp(page, tab, DRAG_FROM_C);
      await stableLine(page, line); // the initial solve is not part of the drag
      await page.evaluate(() => {
        const w = window as unknown as { __longMs: number };
        w.__longMs = 0;
        new PerformanceObserver((list) => {
          for (const e of list.getEntries()) w.__longMs += e.duration;
        }).observe({ type: "longtask", buffered: false });
      });
      await page.locator(".desktop-sim #coffee-temp").focus();

      const start = Date.now();
      for (let i = 0; i < DRAG_PRESSES; i++) {
        await page.keyboard.press("ArrowRight");
      }
      try {
        await page.waitForFunction(
          ([sel, want]) => document.querySelector(sel)?.getAttribute("points") === want,
          [line, expected] as [string, string],
          { polling: "raf", timeout: SETTLE_TIMEOUT_MS }
        );
      } catch {
        throw new Error(
          `interaction ${tab}: chart never settled on the ${DRAG_TO_C} °C line within ${SETTLE_TIMEOUT_MS}ms`
        );
      }
      const elapsed = Date.now() - start;
      // Let trailing long tasks of the settle land in the observer's queue.
      const longMs = await page.evaluate(
        () => new Promise<number>((r) => requestAnimationFrame(() =>
          r((window as unknown as { __longMs: number }).__longMs)))
      );
      return { elapsed, longMs };
    };

    await sample(); // warmup, uncounted

    const times: number[] = [];
    const longs: number[] = [];
    const samples = updating ? SAMPLE_COUNT_BASELINE : SAMPLE_COUNT_CHECK;
    for (let i = 0; i < samples; i++) {
      const { elapsed, longMs } = await sample();
      times.push(elapsed);
      longs.push(longMs);
    }

    const med = median(times);
    const rounded = Math.round(med);
    const longRounded = Math.round(median(longs));
    const key = `interaction:${tab}`;
    interactionResults[key] = { median_ms: rounded, long_task_ms: longRounded };

    if (updating) {
      const current = loadBaselines();
      current[key] = { median_ms: rounded, long_task_ms: longRounded };
      saveBaselines(current);
      console.log(`  ${key} — measured: ${rounded}ms, long tasks ${longRounded}ms (baseline updated)`);
      return;
    }

    const baseline = baselines[key];
    if (!baseline) {
      console.warn(`  ⚠ No baseline for ${key} — skipping assertion (run UPDATE_BASELINES=1 to create)`);
      return;
    }
    const threshold = Math.max(baseline.median_ms * (1 + TOLERANCE), baseline.median_ms + ABS_TOLERANCE_MS);
    const delta = ((med - baseline.median_ms) / baseline.median_ms) * 100;
    console.log(
      `  ${key} — baseline: ${baseline.median_ms}ms, measured: ${rounded}ms (${delta >= 0 ? "+" : ""}${delta.toFixed(1)}%); long tasks ${longRounded}ms (baseline ${baseline.long_task_ms ?? "—"}ms)`
    );
    expect(
      med,
      `${key} regressed: ${rounded}ms vs baseline ${baseline.median_ms}ms (+${delta.toFixed(1)}%, threshold ${Math.round(threshold)}ms)`
    ).toBeLessThanOrEqual(threshold);
  });
}

// Print summary table after all tests
test.afterAll(async () => {
  console.log("\n  ┌─────────────────────────────────────────────┬────────────┬────────────┬─────────┐");
  console.log("  │ Route                                       │ Baseline   │ Measured   │ Delta   │");
  console.log("  ├─────────────────────────────────────────────┼────────────┼────────────┼─────────┤");
  for (const route of routes) {
    const measured = results[route.path];
    const baseline = baselines[route.path];
    const routeCol = route.path.padEnd(43);
    const baseCol = baseline ? `${baseline.median_ms}ms`.padStart(10) : "     —    ";
    const measCol = measured != null ? `${measured}ms`.padStart(10) : "     —    ";
    let deltaCol = "    —   ";
    if (baseline && measured != null) {
      const d = ((measured - baseline.median_ms) / baseline.median_ms) * 100;
      deltaCol = `${d >= 0 ? "+" : ""}${d.toFixed(1)}%`.padStart(7);
    }
    console.log(`  │ ${routeCol} │ ${baseCol} │ ${measCol} │ ${deltaCol} │`);
  }
  console.log("  └─────────────────────────────────────────────┴────────────┴────────────┴─────────┘");

  const shown = interactionTabs.filter(({ tab }) => interactionResults[`interaction:${tab}`]);
  if (shown.length === 0) return;
  console.log("\n  Drag interaction (first press -> chart settled; long tasks are informational)");
  console.log("  ┌───────────────────────────┬────────────┬────────────┬─────────┬────────────┬────────────┐");
  console.log("  │ Tab                       │ Baseline   │ Measured   │ Delta   │ Long (base)│ Long (meas)│");
  console.log("  ├───────────────────────────┼────────────┼────────────┼─────────┼────────────┼────────────┤");
  for (const { tab } of shown) {
    const key = `interaction:${tab}`;
    const m = interactionResults[key];
    const b = baselines[key];
    const ms = (v?: number) => (v != null ? `${v}ms` : "—").padStart(10);
    let deltaCol = "    —   ";
    if (b) {
      const d = ((m.median_ms - b.median_ms) / b.median_ms) * 100;
      deltaCol = `${d >= 0 ? "+" : ""}${d.toFixed(1)}%`.padStart(7);
    }
    console.log(`  │ ${tab.padEnd(25)} │ ${ms(b?.median_ms)} │ ${ms(m.median_ms)} │ ${deltaCol} │ ${ms(b?.long_task_ms)} │ ${ms(m.long_task_ms)} │`);
  }
  console.log("  └───────────────────────────┴────────────┴────────────┴─────────┴────────────┴────────────┘");
});

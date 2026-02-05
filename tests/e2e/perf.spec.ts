import { test, expect } from "@playwright/test";
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
  { path: "/theory/heat-equation?solver=numerical", name: "Theory: Heat Equation (Numerical)", hasSimulation: true },
];

function median(values: number[]): number {
  const sorted = [...values].sort((a, b) => a - b);
  const mid = Math.floor(sorted.length / 2);
  return sorted.length % 2 === 0
    ? (sorted[mid - 1] + sorted[mid]) / 2
    : sorted[mid];
}

function loadBaselines(): Record<string, { median_ms: number }> {
  try {
    const raw = fs.readFileSync(BASELINES_PATH, "utf-8");
    return JSON.parse(raw);
  } catch {
    return {};
  }
}

function saveBaselines(data: Record<string, { median_ms: number }>) {
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
});

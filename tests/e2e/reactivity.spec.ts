import { test, expect, type Page } from "@playwright/test";

// ---------------------------------------------------------------------------
// Reactivity — expensive solves are throttled while a slider is dragged, but
// the chart stays live during the drag and the resting state is exact.
// ---------------------------------------------------------------------------

const SETTINGS_KEY = "prewarm_settings";

// Start the coffee at 70 °C so every one of the drag's ArrowRight presses
// changes the value (the default 90 °C would clamp at 100 °C after ten); pin
// Celsius so the labels read in the units the assertions use.
const SEED = { version: 1, coffee_temp_c: 70, use_fahrenheit: false };

const TABS = [
  {
    tab: "heat_equation",
    // Comparison metrics: the cold/preheated averages and the difference.
    metrics: ".desktop-sim .metrics",
    line: ".desktop-sim polyline.temp-line.hot:not(.mug-temp)",
  },
  {
    tab: "compare",
    // The temperature column only; the compute-time column varies run to run.
    metrics: ".desktop-sim .compare-stats tbody td:nth-child(2)",
    line: ".desktop-sim polyline.heateq-numerical",
  },
];

async function openTab(page: Page, tab: string) {
  await page.goto(`/?tab=${tab}`, { waitUntil: "networkidle" });
  await expect(page.locator(".desktop-sim #coffee-temp")).toBeVisible();
}

async function seededTab(page: Page, tab: string) {
  await openTab(page, tab);
  await page.evaluate(
    ([key, seed]) => localStorage.setItem(key as string, JSON.stringify(seed)),
    [SETTINGS_KEY, SEED] as const
  );
  await openTab(page, tab);
  await expect(page.locator('.desktop-sim label[for="coffee-temp"]')).toContainText("70");
}

/** Everything the resting state shows: the metrics text and the chart line. */
async function restingState(page: Page, metrics: string, line: string) {
  const text = (await page.locator(metrics).allTextContents()).join(" | ");
  const points = await page.locator(line).getAttribute("points");
  return { text, points };
}

for (const { tab, metrics, line } of TABS) {
  test(`resting_state_matches_fresh_load (${tab})`, async ({ page }) => {
    await seededTab(page, tab);

    await page.locator(".desktop-sim #coffee-temp").focus();
    for (let i = 0; i < 25; i++) {
      await page.keyboard.press("ArrowRight");
    }
    await expect(page.locator('.desktop-sim label[for="coffee-temp"]')).toContainText("95");
    await page.waitForTimeout(500);
    const afterDrag = await restingState(page, metrics, line);

    // Settings persist, so a reload solves the same inputs from scratch.
    await openTab(page, tab);
    await expect(page.locator('.desktop-sim label[for="coffee-temp"]')).toContainText("95");
    const fresh = await restingState(page, metrics, line);

    expect(afterDrag).toEqual(fresh);
  });

  test(`chart_updates_during_drag (${tab})`, async ({ page }) => {
    await seededTab(page, tab);

    const slider = page.locator(".desktop-sim #coffee-temp");
    const chartLine = page.locator(line);
    await slider.focus();

    // ~20 inputs/s for 1 s, sampling the chart after each input. The polyline's
    // points come from the throttled solve and an axis derived from the same
    // throttled inputs, so they change only when the throttle publishes: a
    // throttle that stops publishing during a drag leaves them frozen.
    const seen = new Set<string | null>();
    for (let i = 0; i < 20; i++) {
      await page.keyboard.press("ArrowRight");
      seen.add(await chartLine.getAttribute("points"));
      await page.waitForTimeout(50);
    }

    expect(seen.size, "chart does not freeze during a drag").toBeGreaterThanOrEqual(3);
  });
}

// A tab switch while the throttled solve is still catching up unmounts the
// dragged panel mid-flight: its pending timers and effects must die with it,
// and the newly mounted panel must solve the final inputs from scratch.
test("tab_switch_mid_drag_settles", async ({ page, context }) => {
  const errors: string[] = [];
  page.on("console", (msg) => {
    if (msg.type() === "error") errors.push(msg.text());
  });
  page.on("pageerror", (err) => errors.push(err.message));

  await seededTab(page, "heat_equation");
  const { metrics, line } = TABS.find((t) => t.tab === "compare")!;

  await page.locator(".desktop-sim #coffee-temp").focus();
  for (let i = 0; i < 5; i++) {
    await page.keyboard.press("ArrowRight");
  }
  await page.locator(".desktop-sim .solver-toggle button", { hasText: "Compare" }).click();
  await expect(page.locator('.desktop-sim label[for="coffee-temp"]')).toContainText("75");

  // The fresh load reads the persisted settings, so wait for the save first.
  await expect
    .poll(() =>
      page.evaluate((key) => {
        const stored = JSON.parse(localStorage.getItem(key) ?? "{}");
        return [stored.coffee_temp_c, stored.active_tab];
      }, SETTINGS_KEY)
    )
    .toEqual([75, "compare"]);

  const freshPage = await context.newPage();
  await openTab(freshPage, "compare");
  await expect(freshPage.locator('.desktop-sim label[for="coffee-temp"]')).toContainText("75");
  const fresh = await restingState(freshPage, metrics, line);
  await freshPage.close();
  expect(fresh.points, "the fresh load plotted a line").toBeTruthy();

  await expect.poll(() => restingState(page, metrics, line)).toEqual(fresh);
  expect(errors).toEqual([]);
});

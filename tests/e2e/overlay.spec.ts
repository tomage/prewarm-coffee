import { test, expect, type Page } from "@playwright/test";

// ---------------------------------------------------------------------------
// Overlay — exactly one simulator is mounted (desktop inline or the mobile
// overlay), the overlay keeps its controls across tab switches, and it is a
// modal dialog: labelled, focus moved in, trapped and restored, closed by
// Escape and by the browser Back button.
// ---------------------------------------------------------------------------

const DESKTOP = { width: 1280, height: 800 };
const MOBILE = { width: 390, height: 844 };

// Same selector and visibility rule as the overlay's focus trap.
const FOCUSABLE =
  'button:not([disabled]), a[href], input:not([disabled]), select:not([disabled]), textarea:not([disabled]), [tabindex]:not([tabindex="-1"])';

async function focusIsInOverlay(page: Page): Promise<boolean> {
  return page.evaluate(() => {
    const overlay = document.querySelector(".simulator-overlay");
    return !!overlay && overlay.contains(document.activeElement);
  });
}

async function openOverlay(page: Page) {
  await page.getByRole("button", { name: "Open Simulator" }).click();
  await expect(page.locator(".simulator-overlay")).toBeVisible();
}

test.describe("desktop", () => {
  test.use({ viewport: DESKTOP });

  test("one_simulator_mounted_desktop", async ({ page }) => {
    await page.goto("/", { waitUntil: "networkidle" });
    await expect(page.locator(".desktop-sim")).toBeVisible();
    await expect(page.locator(".mobile-sim-entry")).toHaveCount(0);
    await expect(page.locator(".simulator-overlay")).toHaveCount(0);
    await expect(page.locator("#cooling-chart-title")).toHaveCount(1);
  });
});

test.describe("mobile", () => {
  test.use({ viewport: MOBILE });

  test("one_simulator_mounted_mobile", async ({ page }) => {
    await page.goto("/", { waitUntil: "networkidle" });
    await openOverlay(page);
    await expect(page.locator("#cooling-chart-title")).toHaveCount(1);
  });

  test("tab_switch_keeps_overlay_state", async ({ page }) => {
    await page.goto("/", { waitUntil: "networkidle" });
    await openOverlay(page);

    const overlay = page.locator(".simulator-overlay");
    const simulation = overlay.locator(".accordion-group", {
      has: page.locator(".accordion-title", { hasText: "Simulation" }),
    });
    await simulation.locator(".accordion-header").click();
    await expect(simulation).toHaveClass(/\bopen\b/);

    const controls = overlay.locator(".overlay-controls");
    await controls.evaluate((el) => {
      el.scrollTop = 200;
    });
    expect(await controls.evaluate((el) => el.scrollTop)).toBeGreaterThanOrEqual(150);

    await overlay.getByRole("button", { name: "Heat Eq" }).click();
    await expect(overlay.locator(".overlay-tabs button.active")).toHaveText("Heat Eq");

    await expect(simulation).toHaveClass(/\bopen\b/);
    expect(await controls.evaluate((el) => el.scrollTop)).toBeGreaterThanOrEqual(150);
  });

  test("resize_across_breakpoint", async ({ page }) => {
    await page.goto("/", { waitUntil: "networkidle" });
    await openOverlay(page);
    await expect(page.locator("body")).toHaveClass(/sim-overlay-open/);

    await page.setViewportSize(DESKTOP);
    await expect(page.locator(".simulator-overlay")).toHaveCount(0);
    await expect(page.locator(".desktop-sim")).toBeVisible();
    await expect(page.locator("body")).not.toHaveClass(/sim-overlay-open/);

    await page.setViewportSize(MOBILE);
    await expect(page.locator(".mobile-sim-entry")).toBeVisible();
    await expect(page.locator(".simulator-overlay")).toHaveCount(0);
  });

  test("overlay_is_dialog", async ({ page }) => {
    await page.goto("/", { waitUntil: "networkidle" });
    await openOverlay(page);

    const overlay = page.locator(".simulator-overlay");
    await expect(overlay).toHaveAttribute("role", "dialog");
    await expect(overlay).toHaveAttribute("aria-modal", "true");
    const labelId = await overlay.getAttribute("aria-labelledby");
    expect(labelId).toBeTruthy();
    const label = page.locator(`[id="${labelId}"]`);
    await expect(label).toBeVisible();
    expect((await label.innerText()).trim()).not.toBe("");
    await expect(page.getByRole("dialog", { name: (await label.innerText()).trim() })).toBeVisible();
  });

  test("focus_moves_in_and_is_trapped", async ({ page }) => {
    await page.goto("/", { waitUntil: "networkidle" });
    await openOverlay(page);
    await expect.poll(() => focusIsInOverlay(page)).toBe(true);

    for (let i = 0; i < 40; i++) {
      await page.keyboard.press("Tab");
      expect(await focusIsInOverlay(page), `focus left the overlay after ${i + 1} Tab presses`).toBe(true);
    }

    const focusableCount = await page.evaluate((selector) => {
      const overlay = document.querySelector(".simulator-overlay")!;
      const all = [...overlay.querySelectorAll<HTMLElement>(selector)].filter(
        (el) => el.getClientRects().length > 0,
      );
      all[0].dataset.trapProbe = "first";
      all[all.length - 1].dataset.trapProbe = "last";
      return all.length;
    }, FOCUSABLE);
    expect(focusableCount).toBeGreaterThan(1);

    await page.locator('[data-trap-probe="first"]').focus();
    await page.keyboard.press("Shift+Tab");
    expect(await page.evaluate(() => (document.activeElement as HTMLElement).dataset.trapProbe)).toBe("last");

    await page.keyboard.press("Tab");
    expect(await page.evaluate(() => (document.activeElement as HTMLElement).dataset.trapProbe)).toBe("first");
  });

  test("escape_closes_and_restores_focus", async ({ page }) => {
    await page.goto("/", { waitUntil: "networkidle" });
    await openOverlay(page);
    await expect.poll(() => focusIsInOverlay(page)).toBe(true);

    await page.keyboard.press("Escape");
    await expect(page.locator(".simulator-overlay")).toHaveCount(0);
    await expect(page.getByRole("button", { name: "Open Simulator" })).toBeFocused();
    await expect(page.locator("body")).not.toHaveClass(/sim-overlay-open/);
  });

  test("back_closes_overlay_without_reload", async ({ page }) => {
    await page.goto("/", { waitUntil: "networkidle" });
    const url = page.url();
    await page.evaluate(() => {
      (window as unknown as { __marker: number }).__marker = 1;
    });
    await openOverlay(page);

    await page.goBack();
    await expect(page.locator(".simulator-overlay")).toHaveCount(0);
    expect(await page.evaluate(() => (window as unknown as { __marker?: number }).__marker)).toBe(1);
    expect(page.url()).toBe(url);
    await expect(page.getByRole("button", { name: "Open Simulator" })).toBeFocused();
  });

  test("done_then_back_leaves_page", async ({ page }) => {
    await page.goto("/theory", { waitUntil: "networkidle" });
    await expect(page.locator("h1").first()).toHaveText("Theory");
    await page.goto("/", { waitUntil: "networkidle" });
    await expect(page.locator("h1").first()).toHaveText("Prewarm");

    await openOverlay(page);
    await page.locator(".simulator-overlay").getByRole("button", { name: "Done" }).click();
    await expect(page.locator(".simulator-overlay")).toHaveCount(0);
    // Done pops the overlay's history entry; wait until that has happened.
    await expect.poll(() => page.evaluate(() => history.state?.overlay === true)).toBe(false);

    await page.goBack();
    await expect(page.locator("h1").first()).toHaveText("Theory", { timeout: 15_000 });
    expect(page.url()).toMatch(/\/theory\/?$/);
  });

  test("back_after_resize_pops_leftover_entry", async ({ page }) => {
    await page.goto("/", { waitUntil: "networkidle" });
    const url = page.url();
    await openOverlay(page);

    // Leaving the mobile layout unmounts the overlay but leaves its history
    // entry; the next Back pops that entry without reloading the page.
    await page.setViewportSize(DESKTOP);
    await expect(page.locator(".simulator-overlay")).toHaveCount(0);
    await page.evaluate(() => {
      (window as unknown as { __marker: number }).__marker = 1;
    });

    await page.goBack();
    await expect.poll(() => page.evaluate(() => history.state?.overlay === true)).toBe(false);
    expect(await page.evaluate(() => (window as unknown as { __marker?: number }).__marker)).toBe(1);
    expect(page.url()).toBe(url);
    await expect(page.locator(".desktop-sim")).toBeVisible();
  });
});

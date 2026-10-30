import { test, expect } from "@playwright/test";

// ---------------------------------------------------------------------------
// Controls — one MugControls serves the desktop panels and the mobile overlay.
// ---------------------------------------------------------------------------

test("select_follows_reset", async ({ page }) => {
  await page.goto("/", { waitUntil: "networkidle" });
  const material = page.locator(".desktop-sim #material");
  await expect(material).toBeVisible();

  await material.selectOption("glass");
  await expect(material).toHaveValue("glass");

  await page.locator(".desktop-sim .reset-button").click();
  await expect(material).toHaveValue("ceramic");
});

test.describe("mobile", () => {
  test.use({ viewport: { width: 390, height: 844 } });

  test("overlay_has_advanced_controls_for_heat_eq", async ({ page }) => {
    await page.goto("/?tab=heat_equation", { waitUntil: "networkidle" });
    await page.getByRole("button", { name: "Open Simulator" }).click();

    const overlay = page.locator(".simulator-overlay");
    await expect(overlay.locator(".overlay-tabs button.active")).toHaveText("Heat Eq");

    const advanced = overlay.locator(".accordion-group", {
      has: page.locator(".accordion-title", { hasText: "Advanced" }),
    });
    await advanced.locator(".accordion-header").click();
    await expect(advanced.locator("#m-pde-modes")).toBeVisible();
  });
});

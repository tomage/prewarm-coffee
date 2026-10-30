import { test, expect, type Page } from "@playwright/test";

// ---------------------------------------------------------------------------
// Scroll restoration — each page keeps its own scroll position: browser Back
// returns to where the destination page was left, and following a link opens
// the destination at the top.
// ---------------------------------------------------------------------------

const DESKTOP = { width: 1280, height: 800 };
// Back reloads the page (Leptos 0.6 popstate workaround) and the position is
// restored once the app has rendered, so a few px of layout drift are allowed.
const TOLERANCE = 50;

test.use({ viewport: DESKTOP });

async function waitForRoute(page: Page, h1Text: string) {
  await expect(page.locator("h1").first()).toHaveText(h1Text, { timeout: 15_000 });
}

/** Scroll the window to `y` and resolve once the page has seen the scroll event. */
async function scrollTo(page: Page, y: number) {
  await page.evaluate(
    (target) =>
      new Promise<void>((resolve) => {
        window.addEventListener("scroll", () => resolve(), { once: true });
        window.scrollTo(0, target);
      }),
    y,
  );
  expect(await page.evaluate(() => window.scrollY)).toBe(y);
}

/**
 * Follow an in-page link without Playwright first scrolling it into view
 * (that scroll would move the page being left before it is left).
 */
async function clickLink(page: Page, selector: string) {
  await page.locator(selector).first().evaluate((el) => (el as HTMLElement).click());
}

async function scrollY(page: Page): Promise<number> {
  return page.evaluate(() => window.scrollY);
}

test("back_restores_home_position", async ({ page }) => {
  await page.goto("/", { waitUntil: "networkidle" });
  await waitForRoute(page, "Prewarm");
  await scrollTo(page, 1200);

  await clickLink(page, 'a[href="/theory"]');
  await waitForRoute(page, "Theory");
  expect(await scrollY(page)).toBe(0);

  await page.goBack();
  await waitForRoute(page, "Prewarm");
  await expect.poll(() => scrollY(page)).toBeGreaterThanOrEqual(1200 - TOLERANCE);
  expect(await scrollY(page)).toBeLessThanOrEqual(1200 + TOLERANCE);
});

test("theory_keeps_its_own_position", async ({ page }) => {
  await page.goto("/", { waitUntil: "networkidle" });
  await waitForRoute(page, "Prewarm");
  await scrollTo(page, 1200);
  await clickLink(page, 'a[href="/theory"]');
  await waitForRoute(page, "Theory");

  // The theory page scrolls ~480 px at this viewport, so 400 rather than 600.
  await scrollTo(page, 400);
  await clickLink(page, 'a[href="/theory/lumped"]');
  await waitForRoute(page, "Theory: Lumped Capacitance");
  expect(await scrollY(page)).toBe(0);

  await page.goBack();
  await waitForRoute(page, "Theory");
  await expect.poll(() => scrollY(page)).toBeGreaterThanOrEqual(400 - TOLERANCE);
  expect(await scrollY(page)).toBeLessThanOrEqual(400 + TOLERANCE);
});

test("back_from_a_longer_page_keeps_destination_position", async ({ page }) => {
  await page.goto("/theory", { waitUntil: "networkidle" });
  await waitForRoute(page, "Theory");
  await scrollTo(page, 400);
  await clickLink(page, 'a[href="/theory/lumped"]');
  await waitForRoute(page, "Theory: Lumped Capacitance");

  // Far past theory's maximum scroll (~480 px): Back must not save the
  // clamped position under theory while the page reloads.
  await scrollTo(page, 3000);

  await page.goBack();
  await waitForRoute(page, "Theory");
  await expect.poll(() => scrollY(page)).toBeGreaterThanOrEqual(400 - TOLERANCE);
  expect(await scrollY(page)).toBeLessThanOrEqual(400 + TOLERANCE);
});

test("link_click_opens_at_top", async ({ page }) => {
  await page.goto("/", { waitUntil: "networkidle" });
  await waitForRoute(page, "Prewarm");
  await scrollTo(page, 1200);
  await clickLink(page, 'a[href="/theory"]');
  await waitForRoute(page, "Theory");

  // Home has a saved position (1200), but a link opens it at the top.
  await clickLink(page, 'main a[href="/"]');
  await waitForRoute(page, "Prewarm");
  expect(await scrollY(page)).toBe(0);
  // Not just at the top for a moment: no late restore moves it.
  await page.evaluate(() => new Promise((r) => setTimeout(r, 300)));
  expect(await scrollY(page)).toBe(0);
});

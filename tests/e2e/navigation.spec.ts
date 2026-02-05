import { test, expect } from "@playwright/test";

// ---------------------------------------------------------------------------
// Navigation test — exercises the full site graph using both in-page links
// and the browser back/forward buttons, verifying that client-side routing
// works correctly. Back/forward triggers a full page reload (popstate
// workaround for Leptos 0.6), so we use generous timeouts.
// ---------------------------------------------------------------------------

/** Wait for the page to show the expected h1 text (exact match). */
async function waitForRoute(page: import("@playwright/test").Page, h1Text: string) {
  await expect(page.locator("h1").first()).toHaveText(h1Text, { timeout: 15_000 });
}

test("full site navigation via links and browser back button", async ({ page }) => {
  // 1. Start at the home page
  await page.goto("/", { waitUntil: "networkidle" });
  await waitForRoute(page, "Prewarm");
  expect(page.url()).toMatch(/\/$/);

  // 2. Home → Theory (click link)
  await page.click('a[href="/theory"]');
  await waitForRoute(page, "Theory");
  expect(page.url()).toContain("/theory");

  // 3. Theory → Lumped (click link)
  await page.click('a[href="/theory/lumped"]');
  await waitForRoute(page, "Theory: Lumped Capacitance");
  expect(page.url()).toContain("/theory/lumped");

  // 4. Back button: Lumped → Theory
  await page.goBack();
  await waitForRoute(page, "Theory");
  expect(page.url()).toMatch(/\/theory\/?$/);

  // 5. Theory → Heat Equation (click link)
  //    NOTE: pushState clears forward entries — /theory/lumped is gone from history
  await page.click('a[href="/theory/heat-equation"]');
  await waitForRoute(page, "Theory: Heat Equation (PDE)");
  expect(page.url()).toContain("/theory/heat-equation");

  // 6. Use the in-page "< Back to Theory" link
  await page.click('a[href="/theory"]:visible');
  await waitForRoute(page, "Theory");
  expect(page.url()).toMatch(/\/theory\/?$/);

  // 7. Back button: Theory → Heat Equation
  await page.goBack();
  await waitForRoute(page, "Theory: Heat Equation (PDE)");
  expect(page.url()).toContain("/theory/heat-equation");

  // 8. Back button: Heat Equation → Theory
  await page.goBack();
  await waitForRoute(page, "Theory");
  expect(page.url()).toMatch(/\/theory\/?$/);

  // 9. Back button: Theory → Home
  //    History: [/, /theory, /theory/heat-equation, /theory] — going back 3× lands on /
  await page.goBack();
  await waitForRoute(page, "Prewarm");
  expect(page.url()).toMatch(/\/$/);

  // 10. Home → About (click link)
  await page.click('a[href="/about"]');
  await waitForRoute(page, "Prewarm");
  expect(page.url()).toContain("/about");

  // 11. Use the in-page "< Back to simulation" link
  await page.click('a[href="/"]:visible');
  await waitForRoute(page, "Prewarm");
  expect(page.url()).toMatch(/\/$/);

  // 12. Back button: Home → About
  await page.goBack();
  await waitForRoute(page, "Prewarm");
  expect(page.url()).toContain("/about");

  // 13. Back button: About → Home
  await page.goBack();
  await waitForRoute(page, "Prewarm");
  expect(page.url()).toMatch(/\/$/);

  // 14. Forward button: Home → About
  await page.goForward();
  await waitForRoute(page, "Prewarm");
  expect(page.url()).toContain("/about");
});

test("?tab= query parameter selects the correct tab", async ({ page }) => {
  // 1. Go to home page and select Heat Equation tab
  await page.goto("/", { waitUntil: "networkidle" });
  await waitForRoute(page, "Prewarm");
  await page.click('button:has-text("Heat Equation")');
  await expect(page.locator('button:has-text("Heat Equation")')).toHaveClass(/active/);

  // 2. Navigate to /theory
  await page.click('a[href="/theory"]');
  await waitForRoute(page, "Theory");

  // 3. Click the "Compare tab" link which goes to /?tab=compare
  await page.click('a[href="/?tab=compare"]');
  await waitForRoute(page, "Prewarm");

  // 4. Verify the URL has the query parameter and Compare tab is active
  expect(page.url()).toContain("tab=compare");
  await expect(page.locator('button:has-text("Compare")')).toHaveClass(/active/, { timeout: 5_000 });
  await expect(page.locator('button:has-text("Heat Equation")')).not.toHaveClass(/active/);
});

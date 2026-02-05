import { test, expect, type Page, type ConsoleMessage } from "@playwright/test";

const MAX_LOAD_MS = 5_000;

/** Collect console errors and page crashes while a callback runs. */
async function withErrorTracking(
  page: Page,
  fn: () => Promise<void>
): Promise<string[]> {
  const errors: string[] = [];

  const onConsole = (msg: ConsoleMessage) => {
    if (msg.type() === "error") {
      const text = msg.text();
      // Trunk dev server noise: transient rebuild errors and integrity hash mismatches
      if (text.includes("Build failed")) return;
      if (text.includes("integrity") && text.includes("digest")) return;
      errors.push(`console.error: ${text}`);
    }
  };
  const onPageError = (err: Error) => {
    errors.push(`pageerror: ${err.message}`);
  };

  page.on("console", onConsole);
  page.on("pageerror", onPageError);

  await fn();

  page.off("console", onConsole);
  page.off("pageerror", onPageError);

  return errors;
}

// ---------------------------------------------------------------------------
// Route smoke tests
// ---------------------------------------------------------------------------

const routes: {
  path: string;
  name: string;
  h1: string;
  extraCheck?: (page: Page) => Promise<void>;
}[] = [
  {
    path: "/",
    name: "Home",
    h1: "Prewarm",
    extraCheck: async (page) => {
      await expect(page.locator("svg").first()).toBeVisible();
    },
  },
  {
    path: "/about",
    name: "About",
    h1: "Prewarm",
  },
  {
    path: "/theory",
    name: "Theory",
    h1: "Theory",
  },
  {
    path: "/theory/lumped",
    name: "Theory: Lumped",
    h1: "Theory: Lumped Capacitance",
  },
  {
    path: "/theory/heat-equation",
    name: "Theory: Heat Equation",
    h1: "Theory: Heat Equation (PDE)",
  },
  {
    path: "/does-not-exist",
    name: "404",
    h1: "404",
  },
];

for (const route of routes) {
  test(`${route.name} (${route.path}) loads without errors`, async ({
    page,
  }) => {
    const errors = await withErrorTracking(page, async () => {
      const start = Date.now();
      await page.goto(route.path, { waitUntil: "networkidle" });
      const elapsed = Date.now() - start;

      expect(elapsed, `Page load took ${elapsed}ms (max ${MAX_LOAD_MS}ms)`).toBeLessThan(MAX_LOAD_MS);
    });

    // Heading renders
    await expect(page.locator("h1", { hasText: route.h1 })).toBeVisible();

    // Extra assertions (e.g. SVG chart present)
    if (route.extraCheck) {
      await route.extraCheck(page);
    }

    // No WASM panics or JS errors
    expect(errors, "Unexpected browser errors").toEqual([]);
  });
}

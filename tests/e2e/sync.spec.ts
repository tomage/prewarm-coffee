import { test, expect, type Page } from "@playwright/test";

// ---------------------------------------------------------------------------
// Settings sync — settings persist to localStorage and follow across tabs.
// ---------------------------------------------------------------------------

const SETTINGS_KEY = "prewarm_settings";

/** Wait until the app has mounted the desktop controls. */
async function openHome(page: Page) {
  await page.goto("/", { waitUntil: "networkidle" });
  await expect(page.locator("#volume")).toBeVisible();
}

async function storedSettings(page: Page): Promise<Record<string, unknown>> {
  const json = await page.evaluate((key) => localStorage.getItem(key), SETTINGS_KEY);
  expect(json, "settings are stored").not.toBeNull();
  return JSON.parse(json as string);
}

test("cross_tab_follows", async ({ context }) => {
  const pageA = await context.newPage();
  await openHome(pageA);
  // Seed a non-default field and one from a future version; both must survive.
  await pageA.evaluate((key) => {
    localStorage.setItem(
      key,
      JSON.stringify({ version: 1, room_temp_c: 18, future_field: "kept" })
    );
  }, SETTINGS_KEY);
  await openHome(pageA);

  const pageB = await context.newPage();
  await openHome(pageB);
  const volumeLabelB = pageB.locator('label[for="volume"]');
  await expect(volumeLabelB).toHaveText("Volume: 180 ml");

  await pageA.locator("#volume").focus();
  for (let i = 0; i < 10; i++) {
    await pageA.keyboard.press("ArrowRight");
  }

  await expect(volumeLabelB).toHaveText("Volume: 280 ml", { timeout: 2_000 });
  const stored = await storedSettings(pageB);
  expect(stored.volume_ml).toBe(280);
  expect(stored.room_temp_c).toBe(18);
  expect(stored.future_field).toBe("kept");
  expect(stored.material).toBe("ceramic");
  expect(stored.version).toBe(1);
});

test("drag_writes_storage_rarely", async ({ page }) => {
  // Count settings writes from the first script on the page onwards.
  await page.addInitScript((key) => {
    const w = window as unknown as { __settingsWrites: number };
    w.__settingsWrites = 0;
    const setItem = Storage.prototype.setItem;
    Storage.prototype.setItem = function (k: string, v: string) {
      if (this === window.localStorage && k === key) w.__settingsWrites++;
      return setItem.call(this, k, v);
    };
  }, SETTINGS_KEY);
  const writes = () =>
    page.evaluate(() => (window as unknown as { __settingsWrites: number }).__settingsWrites);

  await openHome(page);
  expect(await writes(), "startup writes nothing").toBe(0);

  await page.locator("#volume").focus();
  for (let i = 0; i < 30; i++) {
    await page.keyboard.press("ArrowRight");
  }

  // 180 ml + 30 steps of 10 ml.
  await expect
    .poll(() =>
      page.evaluate((key) => JSON.parse(localStorage.getItem(key) ?? "{}").volume_ml, SETTINGS_KEY)
    )
    .toBe(480);
  expect(await writes()).toBeGreaterThan(0);
  expect(await writes()).toBeLessThanOrEqual(3);
});

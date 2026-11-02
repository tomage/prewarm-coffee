import { test, expect, type BrowserContext, type Page } from "@playwright/test";

// ---------------------------------------------------------------------------
// Settings sync — settings persist to localStorage and follow across tabs.
// ---------------------------------------------------------------------------

const SETTINGS_KEY = "prewarm_settings";

/** Must outlast the app's 250 ms save debounce (`SAVE_DEBOUNCE_MS`). */
const SAVE_DEBOUNCE_WINDOW_MS = 250 + 350;

/** Count settings writes on every page of `context` from its first script on. */
async function countSettingsWrites(context: BrowserContext) {
  await context.addInitScript((key) => {
    const w = window as unknown as { __settingsWrites: number };
    w.__settingsWrites = 0;
    const setItem = Storage.prototype.setItem;
    Storage.prototype.setItem = function (k: string, v: string) {
      if (this === window.localStorage && k === key) w.__settingsWrites++;
      return setItem.call(this, k, v);
    };
  }, SETTINGS_KEY);
}

const writesOf = (page: Page) =>
  page.evaluate(() => (window as unknown as { __settingsWrites: number }).__settingsWrites);

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

test("no_echo_write", async ({ context }) => {
  await countSettingsWrites(context);
  const pageA = await context.newPage();
  await openHome(pageA);
  const pageB = await context.newPage();
  await openHome(pageB);
  const volumeLabelB = pageB.locator('label[for="volume"]');

  await pageA.locator("#volume").focus();
  for (let i = 0; i < 10; i++) {
    await pageA.keyboard.press("ArrowRight");
  }
  await expect(volumeLabelB).toHaveText("Volume: 280 ml", { timeout: 2_000 });
  // A's debounced save has landed (B only follows what A stored).
  expect((await storedSettings(pageA)).volume_ml).toBe(280);

  // Negative assertion: an echo write from B would be armed by the sync that
  // just applied and fire one debounce after it. Waiting one debounce plus
  // margin from A's landed write bounds that.
  await pageB.waitForTimeout(SAVE_DEBOUNCE_WINDOW_MS);
  expect(await writesOf(pageB), "B writes nothing after following A").toBe(0);

  // The focus re-read applies stored settings again; it must not echo either.
  await pageB.evaluate(() => window.dispatchEvent(new Event("focus")));
  await expect(volumeLabelB).toHaveText("Volume: 280 ml");
  await pageB.waitForTimeout(SAVE_DEBOUNCE_WINDOW_MS);
  expect(await writesOf(pageB), "B writes nothing after a focus re-read").toBe(0);
});

test("pagehide_flushes_pending_save", async ({ page }) => {
  await openHome(page);
  await page.locator("#volume").focus();
  await page.keyboard.press("ArrowRight");
  // A synthetic pagehide, read back in the same task: the debounce timer
  // cannot have fired in between, and the visibilitychange flush (which also
  // covers real navigations) is not involved, so only the pagehide listener
  // can account for the stored value.
  const { before, after } = await page.evaluate((key) => {
    const volume = () => JSON.parse(localStorage.getItem(key) ?? "{}").volume_ml;
    const before = volume();
    window.dispatchEvent(new Event("pagehide"));
    return { before, after: volume() };
  }, SETTINGS_KEY);
  expect(before, "the edit is still unsaved before pagehide").toBeUndefined();
  expect(after).toBe(190);
});

test("navigating_away_keeps_pending_save", async ({ page }) => {
  await openHome(page);
  await page.locator("#volume").focus();
  await page.keyboard.press("ArrowRight");
  // Leave before the debounce fires.
  await page.goto("/about", { waitUntil: "networkidle" });
  expect((await storedSettings(page)).volume_ml).toBe(190);
});

test("pending_save_wins_over_sync", async ({ context }) => {
  const pageA = await context.newPage();
  await openHome(pageA);
  const pageB = await context.newPage();
  await openHome(pageB);
  const volumeLabelB = pageB.locator('label[for="volume"]');

  // B edits (190 ml, unsaved), then within the debounce window another tab
  // stores a different volume: B must save its own edit, not follow.
  await pageB.locator("#volume").focus();
  await pageB.keyboard.press("ArrowRight");
  await pageA.evaluate(
    (key) => localStorage.setItem(key, JSON.stringify({ version: 1, volume_ml: 400 })),
    SETTINGS_KEY
  );

  // The storage event reaches B asynchronously and flushes B's edit over 400.
  await expect(volumeLabelB).toHaveText("Volume: 190 ml");
  await expect.poll(async () => (await storedSettings(pageA)).volume_ml).toBe(190);
  // Focus is dispatched only after that settles: sent earlier it could run
  // before the storage event, and the then-stale event would apply 400.
  await pageB.evaluate(() => window.dispatchEvent(new Event("focus")));
  // Past the debounce window nothing has reverted B to the other tab's value.
  await pageB.waitForTimeout(SAVE_DEBOUNCE_WINDOW_MS);
  await expect(volumeLabelB).toHaveText("Volume: 190 ml");
  expect((await storedSettings(pageA)).volume_ml).toBe(190);
});

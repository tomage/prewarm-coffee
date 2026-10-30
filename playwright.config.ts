import { defineConfig } from "@playwright/test";

const TEST_PORT = 8089;
const externalURL = process.env.BASE_URL;

export default defineConfig({
  testDir: "./tests/e2e",
  timeout: 30_000,
  expect: { timeout: 10_000 },
  fullyParallel: true,
  retries: 0,
  reporter: "list",
  use: {
    baseURL: externalURL || `http://localhost:${TEST_PORT}`,
    headless: true,
  },
  projects: [
    {
      name: "smoke",
      testMatch: /smoke\.spec\.ts/,
      use: { browserName: "chromium" },
    },
    {
      name: "navigation",
      testMatch: /navigation\.spec\.ts/,
      dependencies: ["smoke"],
      use: { browserName: "chromium" },
    },
    {
      // Behaviour specs (settings sync, reactivity, controls, overlay, scroll):
      // one project, so `make test-e2e` runs them all with one flag.
      name: "behaviour",
      testMatch: [/sync\.spec\.ts/, /reactivity\.spec\.ts/, /controls\.spec\.ts/, /overlay\.spec\.ts/, /scroll\.spec\.ts/],
      dependencies: ["smoke"],
      use: { browserName: "chromium" },
    },
    {
      name: "perf",
      testMatch: /perf\.spec\.ts/,
      dependencies: ["smoke"],
      timeout: 120_000,
      use: { browserName: "chromium" },
    },
  ],
  // Skip local dev server when testing against an external URL
  ...(!externalURL && {
    webServer: {
      // The run writes target-test/, dist-test/ and test-results/. --ignore keeps
      // trunk from rebuilding on them: a rebuild replaces the served files while
      // a test loads a page, and the app then fails to mount. trunk refuses to
      // ignore a path that does not exist, hence the mkdir. --no-autoreload: a
      // live reload would navigate a page mid-test and destroy its state (e.g.
      // the sync spec's write counter).
      command:
        `mkdir -p target-test dist-test test-results && ` +
        `CARGO_TARGET_DIR=target-test trunk serve --port ${TEST_PORT} --dist dist-test --no-autoreload ` +
        `--ignore ./target-test --ignore ./dist-test --ignore ./test-results`,
      url: `http://localhost:${TEST_PORT}`,
      reuseExistingServer: false,
      timeout: 120_000,
    },
  }),
});

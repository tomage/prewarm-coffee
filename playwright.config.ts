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
      command: `CARGO_TARGET_DIR=target-test trunk serve --port ${TEST_PORT} --dist dist-test`,
      url: `http://localhost:${TEST_PORT}`,
      reuseExistingServer: false,
      timeout: 120_000,
    },
  }),
});

// @ts-check

/** @type {import('@playwright/test').PlaywrightTestConfig} */
const config = {
  testDir: '.',
  timeout: 60000,
  retries: 1,
  reporter: 'list',
  use: {
    baseURL: 'http://localhost:3000',
    headless: false,
    viewport: { width: 1280, height: 720 },
    ignoreHTTPSErrors: true,
    video: 'on-first-retry',
    screenshot: 'on',
    trace: 'on-first-retry',
  },
  // Let's not try to start the server automatically since it's causing issues
  // Instead, we'll document that users should start the server manually
  webServer: null,
  projects: [
    {
      name: 'chromium',
      use: {
        browserName: 'chromium',
      },
    },
  ],
};

module.exports = config;

const { chromium } = require('playwright');
const fs = require('fs');

(async () => {
  const browser = await chromium.launch();
  const page = await browser.newPage();
  
  try {
    console.log('Navigating to https://example.com');
    await page.goto('https://example.com');
    console.log('Getting accessibility snapshot');
    const snapshot = await page.accessibility.snapshot();
    console.log(JSON.stringify({ success: true, data: { snapshot } }));
  } catch (error) {
    console.error(error);
    console.log(JSON.stringify({ success: false, error: error.message }));
  } finally {
    await browser.close();
  }
})();


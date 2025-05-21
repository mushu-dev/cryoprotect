const { chromium } = require('playwright');
const fs = require('fs');

(async () => {
  const browser = await chromium.launch();
  const page = await browser.newPage();
  
  try {
    await page.goto('https://example.com');
    const snapshot = await page.accessibility.snapshot();
    const result = { success: true, data: { snapshot } };
    fs.writeFileSync('/app/shared/snapshot_result.json', JSON.stringify(result));
  } catch (error) {
    const result = { success: false, error: error.message };
    fs.writeFileSync('/app/shared/snapshot_result.json', JSON.stringify(result));
  } finally {
    await browser.close();
  }
})();

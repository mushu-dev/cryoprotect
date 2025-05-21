const { chromium } = require('playwright');
const fs = require('fs');

(async () => {
  const browser = await chromium.launch();
  const page = await browser.newPage();
  
  try {
    await page.goto('https://example.com');
    await page.screenshot({ path: '/app/shared/screenshot.png' });
    const result = { success: true, data: { filename: 'screenshot.png' } };
    fs.writeFileSync('/app/shared/screenshot_result.json', JSON.stringify(result));
  } catch (error) {
    const result = { success: false, error: error.message };
    fs.writeFileSync('/app/shared/screenshot_result.json', JSON.stringify(result));
  } finally {
    await browser.close();
  }
})();

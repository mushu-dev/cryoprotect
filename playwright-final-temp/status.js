// status.js
const { chromium } = require('playwright');

(async () => {
  try {
    const browser = await chromium.launch();
    await browser.close();
    console.log(JSON.stringify({ success: true, status: "Playwright is working correctly" }));
  } catch (error) {
    console.error('Error:', error.message);
    console.log(JSON.stringify({ success: false, error: error.message }));
  }
})();

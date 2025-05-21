
const { chromium } = require('playwright');
const fs = require('fs');

(async () => {
  const browser = await chromium.launch();
  const page = await browser.newPage();
  
  try {
    console.log('Navigating to https://example.com');
    await page.goto('https://example.com');
    const title = await page.title();
    console.log(JSON.stringify({ success: true, data: { title, url: page.url() } }));
  } catch (error) {
    console.error(error);
    console.log(JSON.stringify({ success: false, error: error.message }));
  } finally {
    await browser.close();
  }
})();


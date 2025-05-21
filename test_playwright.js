// test_playwright.js
const { chromium } = require('playwright');

(async () => {
  console.log('Starting Playwright test...');
  try {
    const browser = await chromium.launch();
    console.log('Browser launched successfully!');
    
    const page = await browser.newPage();
    console.log('Page created successfully!');
    
    await page.goto('https://example.com');
    console.log('Navigated to example.com successfully!');
    
    await page.screenshot({ path: 'example.png' });
    console.log('Screenshot taken successfully!');
    
    await browser.close();
    console.log('Test completed successfully!');
  } catch (error) {
    console.error('Error running test:', error);
  }
})();
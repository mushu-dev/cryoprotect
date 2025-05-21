// Simple Playwright test
const { chromium } = require('playwright');

(async () => {
  console.log('Starting simple Playwright test...');
  
  try {
    console.log('Launching browser...');
    const browser = await chromium.launch();
    
    console.log('Creating page...');
    const page = await browser.newPage();
    
    console.log('Navigating to example.com...');
    await page.goto('https://example.com');
    
    console.log('Getting page title...');
    const title = await page.title();
    console.log(`Page title: ${title}`);
    
    console.log('Closing browser...');
    await browser.close();
    
    console.log('Test completed successfully!');
  } catch (error) {
    console.error('Error:', error);
  }
})();
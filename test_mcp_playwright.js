// test_mcp_playwright.js
// This file tests the MCP Playwright integration

console.log('MCP Playwright Test');
console.log('LD_LIBRARY_PATH:', process.env.LD_LIBRARY_PATH);

// In a real MCP integration, you would use their API
// This is just to confirm our environment is working
const { chromium } = require('playwright');

(async () => {
  try {
    console.log('Launching browser...');
    const browser = await chromium.launch();
    
    console.log('Creating page...');
    const page = await browser.newPage();
    
    console.log('Navigating to website...');
    await page.goto('https://example.com');
    
    console.log('Taking screenshot...');
    await page.screenshot({ path: 'mcp_example.png' });
    
    console.log('Closing browser...');
    await browser.close();
    
    console.log('Test successful!');
  } catch (error) {
    console.error('Test failed:', error);
  }
})();
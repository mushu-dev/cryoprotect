// navigate.js
const { chromium } = require('playwright');

(async () => {
  const browser = await chromium.launch();
  const page = await browser.newPage();
  
  try {
    console.log('Navigating to URL: https://cryoprotect.app/molecules/1');
    await page.goto('https://cryoprotect.app/molecules/1');
    const title = await page.title();
    console.log(JSON.stringify({ 
      success: true, 
      data: { 
        title: title, 
        url: page.url() 
      } 
    }));
  } catch (error) {
    console.error('Error:', error.message);
    console.log(JSON.stringify({ 
      success: false, 
      error: error.message 
    }));
  } finally {
    await browser.close();
  }
})();

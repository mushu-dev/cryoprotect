
// MCP Playwright operation: browser_take_screenshot
// Generated at: 2025-05-17T04:07:00.098Z
// Operation ID: 7fb3f45c1a0c6758

const fs = require('fs');
const path = require('path');
const { chromium, firefox, webkit } = require('playwright');

// Operation parameters
const params = {
  "url": "https://example.com",
  "filename": "example-screenshot.png"
};

// Result file
const resultPath = path.join(__dirname, 'result-7fb3f45c1a0c6758.json');

// Log function
function log(message) {
  const logPath = path.join(__dirname, 'operation-7fb3f45c1a0c6758.log');
  fs.appendFileSync(logPath, message + '\n');
}

// Main operation function
async function performOperation() {
  let browser;
  let page;
  let result = { success: false, error: null, data: null };
  
  try {
    log('Starting operation browser_take_screenshot');
    
    // Default to chromium if not specified
    const browserType = params.browserType || 'chromium';
    log(`Using browser: ${browserType}`);
    
    // Launch browser based on type
    switch (browserType) {
      case 'firefox':
        browser = await firefox.launch({ headless: true });
        break;
      case 'webkit':
        browser = await webkit.launch({ headless: true });
        break;
      case 'chromium':
      default:
        browser = await chromium.launch({ headless: true });
        break;
    }
    
    log('Browser launched');
    page = await browser.newPage();
    log('Page created');

    // Take screenshot
    const screenshotPath = path.join(__dirname, 'screenshot-7fb3f45c1a0c6758.png');
    log(`Taking screenshot: ${screenshotPath}`);
    
    // Navigate to URL if provided
    if (params.url) {
      log(`Navigating to: ${params.url}`);
      await page.goto(params.url, { waitUntil: 'domcontentloaded' });
    }
    
    // Take the screenshot
    await page.screenshot({ path: screenshotPath, fullPage: true });
    result.success = true;
    result.data = { 
      screenshotPath: screenshotPath,
      relativePath: `screenshot-7fb3f45c1a0c6758.png`
    };

  } catch (error) {
    log(`Error during operation: ${error.message}`);
    result.success = false;
    result.error = error.message;
  } finally {
    if (page) {
      log('Closing page');
      await page.close();
    }
    if (browser) {
      log('Closing browser');
      await browser.close();
    }
    
    log('Writing result file');
    fs.writeFileSync(resultPath, JSON.stringify(result, null, 2));
    
    log('Operation completed');
    return result;
  }
}

// Execute operation
performOperation()
  .then(result => {
    log(`Operation result: ${result.success ? 'SUCCESS' : 'FAILURE'}`);
    if (!result.success) {
      log(`Error: ${result.error}`);
    }
  })
  .catch(error => {
    log(`Unexpected error: ${error.message}`);
    const failResult = { success: false, error: error.message, data: null };
    fs.writeFileSync(resultPath, JSON.stringify(failResult, null, 2));
  });

#!/usr/bin/env node

/**
 * MCP Playwright Bridge
 * 
 * This script provides a bridge between MCP Playwright commands and the containerized Playwright environment.
 * It intercepts MCP Playwright commands, translates them to containerized operations, and returns the results.
 */

const fs = require('fs');
const path = require('path');
const { execSync, spawn } = require('child_process');
const crypto = require('crypto');

// Constants
const SCRIPT_DIR = __dirname;
const SHARED_DIR = path.join(SCRIPT_DIR, 'playwright-shared');
const LOG_FILE = path.join(SCRIPT_DIR, 'mcp-playwright-bridge.log');

// Ensure shared directory exists
if (!fs.existsSync(SHARED_DIR)) {
  fs.mkdirSync(SHARED_DIR, { recursive: true });
}

// Logging function
function log(level, message) {
  const timestamp = new Date().toISOString();
  const logEntry = `[${timestamp}] [${level}] ${message}\n`;
  fs.appendFileSync(LOG_FILE, logEntry);
  
  // Also print to console
  const colors = {
    INFO: '\x1b[32m',  // Green
    WARN: '\x1b[33m',  // Yellow
    ERROR: '\x1b[31m', // Red
    DEBUG: '\x1b[34m'  // Blue
  };
  const color = colors[level] || '\x1b[0m';
  console.log(`${color}[${level}] ${message}\x1b[0m`);
}

// Run a command and return the result
function runCommand(command) {
  try {
    log('INFO', `Executing command: ${command}`);
    const result = execSync(command, { encoding: 'utf8' });
    return { success: true, output: result.trim() };
  } catch (error) {
    log('ERROR', `Command execution failed: ${error.message}`);
    return { success: false, error: error.message };
  }
}

// Check container status
function checkContainerStatus() {
  const result = runCommand(`${SCRIPT_DIR}/playwright-container.sh status`);
  return result.success;
}

// Start container if needed
function ensureContainerRunning() {
  if (!checkContainerStatus()) {
    log('INFO', 'Starting Playwright container');
    return runCommand(`${SCRIPT_DIR}/playwright-container.sh start`).success;
  }
  return true;
}

// Generate a unique ID for the operation
function generateOperationId() {
  return crypto.randomBytes(8).toString('hex');
}

// Create a Playwright script for the MCP operation
function createPlaywrightScript(operation, params) {
  const operationId = generateOperationId();
  const scriptPath = path.join(SHARED_DIR, `mcp-operation-${operationId}.js`);
  
  let scriptContent = `
// MCP Playwright operation: ${operation}
// Generated at: ${new Date().toISOString()}
// Operation ID: ${operationId}

const fs = require('fs');
const path = require('path');
const { chromium, firefox, webkit } = require('playwright');

// Operation parameters
const params = ${JSON.stringify(params, null, 2)};

// Result file
const resultPath = path.join(__dirname, 'result-${operationId}.json');

// Log function
function log(message) {
  const logPath = path.join(__dirname, 'operation-${operationId}.log');
  fs.appendFileSync(logPath, message + '\\n');
}

// Main operation function
async function performOperation() {
  let browser;
  let page;
  let result = { success: false, error: null, data: null };
  
  try {
    log('Starting operation ${operation}');
    
    // Default to chromium if not specified
    const browserType = params.browserType || 'chromium';
    log(\`Using browser: \${browserType}\`);
    
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
`;

  // Add operation-specific code
  switch (operation) {
    case 'browser_navigate':
      scriptContent += `
    // Navigate to URL
    log(\`Navigating to: \${params.url}\`);
    await page.goto(params.url, { waitUntil: 'domcontentloaded' });
    result.success = true;
    result.data = { url: page.url(), title: await page.title() };
`;
      break;
      
    case 'browser_take_screenshot':
      scriptContent += `
    // Take screenshot
    const screenshotPath = path.join(__dirname, 'screenshot-${operationId}.png');
    log(\`Taking screenshot: \${screenshotPath}\`);
    
    // Navigate to URL if provided
    if (params.url) {
      log(\`Navigating to: \${params.url}\`);
      await page.goto(params.url, { waitUntil: 'domcontentloaded' });
    }
    
    // Take the screenshot
    await page.screenshot({ path: screenshotPath, fullPage: true });
    result.success = true;
    result.data = { 
      screenshotPath: screenshotPath,
      relativePath: \`screenshot-${operationId}.png\`
    };
`;
      break;
      
    case 'browser_snapshot':
      scriptContent += `
    // Get page snapshot (accessibility tree)
    log('Capturing accessibility snapshot');
    
    // Navigate to URL if provided
    if (params.url) {
      log(\`Navigating to: \${params.url}\`);
      await page.goto(params.url, { waitUntil: 'domcontentloaded' });
    }
    
    // Get snapshot
    const snapshot = await page.accessibility.snapshot();
    result.success = true;
    result.data = { snapshot };
`;
      break;
      
    case 'browser_click':
      scriptContent += `
    // Click on element
    log(\`Clicking element: \${params.element || params.ref}\`);
    
    // Navigate to URL if provided
    if (params.url) {
      log(\`Navigating to: \${params.url}\`);
      await page.goto(params.url, { waitUntil: 'domcontentloaded' });
    }
    
    // Click the element
    if (params.ref) {
      await page.click(\`[data-testid="\${params.ref}"]\`);
    } else if (params.element) {
      await page.click(params.element);
    } else {
      throw new Error('No element or ref provided for click');
    }
    
    result.success = true;
    result.data = { clicked: true };
`;
      break;
      
    case 'browser_type':
      scriptContent += `
    // Type text into element
    log(\`Typing text in element: \${params.element || params.ref}\`);
    
    // Navigate to URL if provided
    if (params.url) {
      log(\`Navigating to: \${params.url}\`);
      await page.goto(params.url, { waitUntil: 'domcontentloaded' });
    }
    
    // Type into the element
    if (params.ref) {
      await page.fill(\`[data-testid="\${params.ref}"]\`, params.text);
    } else if (params.element) {
      await page.fill(params.element, params.text);
    } else {
      throw new Error('No element or ref provided for typing');
    }
    
    // Submit if requested
    if (params.submit) {
      log('Submitting form');
      await page.press(params.element || \`[data-testid="\${params.ref}"]\`, 'Enter');
    }
    
    result.success = true;
    result.data = { typed: true, text: params.text };
`;
      break;
      
    default:
      scriptContent += `
    // Unknown operation
    throw new Error(\`Unsupported operation: ${operation}\`);
`;
      break;
  }

  // Close resources and return result
  scriptContent += `
  } catch (error) {
    log(\`Error during operation: \${error.message}\`);
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
    log(\`Operation result: \${result.success ? 'SUCCESS' : 'FAILURE'}\`);
    if (!result.success) {
      log(\`Error: \${result.error}\`);
    }
  })
  .catch(error => {
    log(\`Unexpected error: \${error.message}\`);
    const failResult = { success: false, error: error.message, data: null };
    fs.writeFileSync(resultPath, JSON.stringify(failResult, null, 2));
  });
`;

  // Write the script to the file
  fs.writeFileSync(scriptPath, scriptContent);
  log('INFO', `Created operation script: ${scriptPath}`);
  return { scriptPath, operationId };
}

// Execute a Playwright operation through the container
async function executePlaywrightOperation(operation, params) {
  log('INFO', `Executing Playwright operation: ${operation}`);
  log('DEBUG', `Operation parameters: ${JSON.stringify(params)}`);
  
  // Ensure container is running
  if (!ensureContainerRunning()) {
    log('ERROR', 'Failed to start Playwright container');
    return { success: false, error: 'Failed to start Playwright container' };
  }
  
  // Create script for the operation
  const { scriptPath, operationId } = createPlaywrightScript(operation, params);
  
  // Execute the script in the container
  const scriptCmd = `${SCRIPT_DIR}/playwright-container.sh run-script ${scriptPath} 2>&1 | grep -v '\\[INFO\\]' | grep -v '\\[DEBUG\\]'`;
  const scriptResult = runCommand(scriptCmd);
  
  if (!scriptResult.success) {
    log('ERROR', `Failed to execute operation script: ${scriptResult.error}`);
    return { success: false, error: `Script execution failed: ${scriptResult.error}` };
  }
  
  // Read the result file
  const resultPath = path.join(SHARED_DIR, `result-${operationId}.json`);
  if (!fs.existsSync(resultPath)) {
    log('ERROR', `Result file not found: ${resultPath}`);
    return { success: false, error: 'Result file not found' };
  }
  
  try {
    const resultContent = fs.readFileSync(resultPath, 'utf8');
    const result = JSON.parse(resultContent);
    
    // Handle screenshot paths if present
    if (result.success && result.data && result.data.relativePath) {
      result.data.fullPath = path.join(SHARED_DIR, result.data.relativePath);
      
      // Copy screenshot to original requested location if specified
      if (params.filename && result.data.screenshotPath) {
        fs.copyFileSync(result.data.screenshotPath, params.filename);
        log('INFO', `Copied screenshot to requested location: ${params.filename}`);
        result.data.userPath = params.filename;
      }
    }
    
    log('INFO', `Operation completed with result: ${result.success ? 'SUCCESS' : 'FAILURE'}`);
    return result;
  } catch (error) {
    log('ERROR', `Failed to parse result: ${error.message}`);
    return { success: false, error: `Failed to parse result: ${error.message}` };
  }
}

// Main function
async function main() {
  const args = process.argv.slice(2);
  
  if (args.length < 1) {
    console.log('MCP Playwright Bridge');
    console.log('Usage: node mcp-playwright-bridge.js <operation> [params]');
    console.log('');
    console.log('Operations:');
    console.log('  browser_navigate         Navigate to a URL');
    console.log('  browser_take_screenshot  Take a screenshot');
    console.log('  browser_snapshot         Get accessibility snapshot');
    console.log('  browser_click            Click on an element');
    console.log('  browser_type             Type text into an element');
    console.log('');
    process.exit(0);
  }
  
  const operation = args[0];
  let params = {};
  
  if (args.length > 1) {
    try {
      params = JSON.parse(args[1]);
    } catch (error) {
      log('ERROR', `Failed to parse parameters: ${error.message}`);
      console.error(`Failed to parse parameters: ${error.message}`);
      process.exit(1);
    }
  }
  
  const result = await executePlaywrightOperation(operation, params);
  console.log(JSON.stringify(result, null, 2));
}

// Run the main function
main().catch(error => {
  log('ERROR', `Unhandled error: ${error.message}`);
  console.error(`Unhandled error: ${error.message}`);
  process.exit(1);
});
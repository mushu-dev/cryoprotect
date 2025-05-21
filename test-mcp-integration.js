// test-mcp-integration.js
// Tests if our MCP Playwright integration works correctly

const { execSync } = require('child_process');
const fs = require('fs');

console.log('MCP Playwright Integration Test');
console.log('===============================');

// Function to run a Playwright command
function runCommand(command) {
  console.log(`\nRunning: ${command}`);
  try {
    const output = execSync(command, { encoding: 'utf8' });
    console.log('Output:');
    console.log(output.slice(0, 500)); // Show first 500 chars if output is long
    return { success: true, output };
  } catch (error) {
    console.error('Error:', error.message);
    return { success: false, error: error.message };
  }
}

// Test navigation
console.log('\nTest: Browser navigation');
const navResult = runCommand('./mcp-playwright-final.sh browser_navigate https://example.com');

// Test screenshot
console.log('\nTest: Taking a screenshot');
const screenshotResult = runCommand('./mcp-playwright-final.sh browser_take_screenshot https://example.com mcp-test-screenshot.png');

// Check if screenshot was created
if (fs.existsSync('mcp-test-screenshot.png')) {
  console.log('✅ Screenshot created successfully');
} else {
  console.log('❌ Screenshot was not created');
}

// Test snapshot
console.log('\nTest: Getting accessibility snapshot');
const snapshotResult = runCommand('./mcp-playwright-final.sh browser_snapshot https://example.com');

// Summary
console.log('\nTest Summary');
console.log('===============================');
console.log(`Navigation test: ${navResult.success ? '✅ Passed' : '❌ Failed'}`);
console.log(`Screenshot test: ${screenshotResult.success && fs.existsSync('mcp-test-screenshot.png') ? '✅ Passed' : '❌ Failed'}`);
console.log(`Snapshot test: ${snapshotResult.success ? '✅ Passed' : '❌ Failed'}`);

// Clean up
if (fs.existsSync('mcp-test-screenshot.png')) {
  fs.unlinkSync('mcp-test-screenshot.png');
}

console.log('\nTest completed!');
// test-mcp-playwright-container.js
// Test script for the containerized MCP Playwright solution

console.log('Starting MCP Playwright Container Test');
console.log('========================================');

const { execSync } = require('child_process');
const fs = require('fs');
const path = require('path');

// Test cases
const tests = [
  {
    name: 'Container Startup',
    command: './playwright-container.sh start',
    validate: (output) => output.includes('Playwright container is ready')
  },
  {
    name: 'Container Status',
    command: './playwright-container.sh status',
    validate: (output) => output.includes('Container is running') && output.includes('Playwright is functioning correctly')
  },
  {
    name: 'Browser Navigate',
    command: './mcp-playwright.sh browser_navigate \'{"url": "https://example.com"}\'',
    validate: (output) => {
      const result = JSON.parse(output);
      return result.success && result.data && result.data.title.includes('Example');
    }
  },
  {
    name: 'Take Screenshot',
    command: './mcp-playwright.sh browser_take_screenshot \'{"url": "https://example.com", "filename": "example-screenshot.png"}\'',
    validate: (output) => {
      const result = JSON.parse(output);
      return result.success && fs.existsSync('example-screenshot.png');
    }
  },
  {
    name: 'Get Accessibility Snapshot',
    command: './mcp-playwright.sh browser_snapshot \'{"url": "https://example.com"}\'',
    validate: (output) => {
      const result = JSON.parse(output);
      return result.success && result.data && result.data.snapshot;
    }
  }
];

// Run tests
async function runTests() {
  let passedCount = 0;
  let failedCount = 0;
  
  for (const test of tests) {
    process.stdout.write(`Testing: ${test.name}... `);
    
    try {
      const output = execSync(test.command, { encoding: 'utf8' });
      const passed = test.validate(output);
      
      if (passed) {
        console.log('\x1b[32mPASSED\x1b[0m');
        passedCount++;
      } else {
        console.log('\x1b[31mFAILED\x1b[0m');
        console.log(`Output: ${output}`);
        failedCount++;
      }
    } catch (error) {
      console.log('\x1b[31mERROR\x1b[0m');
      console.log(`Error: ${error.message}`);
      failedCount++;
    }
  }
  
  // Print summary
  console.log('========================================');
  console.log(`Test summary: ${passedCount} passed, ${failedCount} failed`);
  
  if (failedCount === 0) {
    console.log('\x1b[32mAll tests passed! The containerized MCP Playwright solution is working correctly.\x1b[0m');
  } else {
    console.log('\x1b[31mSome tests failed. Please review the errors above.\x1b[0m');
  }
  
  // Clean up
  try {
    if (fs.existsSync('example-screenshot.png')) {
      fs.unlinkSync('example-screenshot.png');
    }
  } catch (error) {
    console.log(`Error during cleanup: ${error.message}`);
  }
  
  return failedCount === 0;
}

// Run the tests
runTests()
  .then(success => {
    process.exit(success ? 0 : 1);
  })
  .catch(error => {
    console.error(`Unhandled error: ${error.message}`);
    process.exit(1);
  });
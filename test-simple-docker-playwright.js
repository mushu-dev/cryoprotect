// test-simple-docker-playwright.js
// Test script for the simplified Playwright Docker MCP solution

console.log('Starting Playwright Docker MCP Test');
console.log('=======================================');

const { execSync } = require('child_process');
const fs = require('fs');
const path = require('path');

// Test cases
const tests = [
  {
    name: 'Container Status',
    command: './playwright-docker-mcp.sh status',
    validate: (output) => {
      try {
        const result = JSON.parse(output);
        return result.success === true;
      } catch (e) {
        console.error('Failed to parse JSON output:', e.message);
        return false;
      }
    }
  },
  {
    name: 'Browser Navigate',
    command: './playwright-docker-mcp.sh browser_navigate https://example.com',
    validate: (output) => {
      try {
        const result = JSON.parse(output);
        return result.success === true && result.data && result.data.title.includes('Example');
      } catch (e) {
        console.error('Failed to parse JSON output:', e.message);
        return false;
      }
    }
  },
  {
    name: 'Take Screenshot',
    command: './playwright-docker-mcp.sh browser_take_screenshot https://example.com test-screenshot.png',
    validate: (output) => {
      try {
        const result = JSON.parse(output);
        return result.success === true && fs.existsSync('test-screenshot.png');
      } catch (e) {
        console.error('Failed to parse JSON output:', e.message);
        return false;
      }
    }
  },
  {
    name: 'Get Accessibility Snapshot',
    command: './playwright-docker-mcp.sh browser_snapshot https://example.com',
    validate: (output) => {
      try {
        const result = JSON.parse(output);
        return result.success === true && result.data && result.data.snapshot;
      } catch (e) {
        console.error('Failed to parse JSON output:', e.message);
        return false;
      }
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
  console.log('=======================================');
  console.log(`Test summary: ${passedCount} passed, ${failedCount} failed`);
  
  if (failedCount === 0) {
    console.log('\x1b[32mAll tests passed! The Playwright Docker MCP solution is working correctly.\x1b[0m');
  } else {
    console.log('\x1b[31mSome tests failed. Please review the errors above.\x1b[0m');
  }
  
  // Clean up
  try {
    if (fs.existsSync('test-screenshot.png')) {
      fs.unlinkSync('test-screenshot.png');
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
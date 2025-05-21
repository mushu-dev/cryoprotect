// test-direct-playwright.js
// Test script for the direct Playwright solution

console.log('Starting Direct Playwright Test');
console.log('=======================================');

const { execSync } = require('child_process');
const fs = require('fs');
const path = require('path');

// Helper function to extract JSON from output
function extractJson(output) {
  const jsonRegex = /\{[\s\S]*\}/;
  const match = output.match(jsonRegex);
  if (match) {
    try {
      return JSON.parse(match[0]);
    } catch (e) {
      console.error('Failed to parse JSON from output:', match[0].substring(0, 100) + '...');
      return null;
    }
  }
  return null;
}

// Test cases
const tests = [
  {
    name: 'Playwright Status',
    command: './mcp-playwright-direct.sh status',
    validate: (output) => {
      const json = extractJson(output);
      return json && json.success === true;
    }
  },
  {
    name: 'Browser Navigate',
    command: './mcp-playwright-direct.sh browser_navigate https://example.com',
    validate: (output) => {
      const json = extractJson(output);
      return json && json.success === true && json.data && json.data.title.includes('Example');
    }
  },
  {
    name: 'Take Screenshot',
    command: './mcp-playwright-direct.sh browser_take_screenshot https://example.com direct-test.png',
    validate: (output) => {
      const json = extractJson(output);
      return json && json.success === true && fs.existsSync('direct-test.png');
    }
  },
  {
    name: 'Get Accessibility Snapshot',
    command: './mcp-playwright-direct.sh browser_snapshot https://example.com',
    validate: (output) => {
      const json = extractJson(output);
      return json && json.success === true && json.data && json.data.snapshot;
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
    console.log('\x1b[32mAll tests passed! The Direct Playwright solution is working correctly.\x1b[0m');
  } else {
    console.log('\x1b[31mSome tests failed. Please review the errors above.\x1b[0m');
  }
  
  // Clean up
  try {
    if (fs.existsSync('direct-test.png')) {
      fs.unlinkSync('direct-test.png');
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
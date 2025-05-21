// test_mcp_services.js
// A comprehensive test script for all MCP services
// Run with: ./mcp-cli.sh node test_mcp_services.js

const fs = require('fs');
const path = require('path');
const { chromium } = require('playwright');

// Results object to store test outcomes
const results = {
  timestamp: new Date().toISOString(),
  tests: {},
  summary: {
    passed: 0,
    failed: 0,
    skipped: 0
  }
};

// Test runner function
async function runTest(name, testFn) {
  console.log(`\n📋 Running test: ${name}`);
  console.log('--------------------------------------------------');
  
  try {
    await testFn();
    results.tests[name] = { status: 'passed', error: null };
    results.summary.passed++;
    console.log(`✅ Test passed: ${name}`);
  } catch (error) {
    results.tests[name] = { status: 'failed', error: error.message };
    results.summary.failed++;
    console.error(`❌ Test failed: ${name}`);
    console.error(`   Error: ${error.message}`);
  }
}

// Skip a test and record it
function skipTest(name, reason) {
  console.log(`\n📋 Skipping test: ${name}`);
  console.log(`   Reason: ${reason}`);
  results.tests[name] = { status: 'skipped', reason };
  results.summary.skipped++;
}

// Main test suite
async function runTests() {
  console.log('🔍 Starting MCP Services Test Suite');
  console.log('==================================================');
  
  // Test 1: Playwright Browser functionality
  await runTest('Playwright Browser', async () => {
    const browser = await chromium.launch();
    const page = await browser.newPage();
    await page.goto('https://example.com');
    const title = await page.title();
    await page.screenshot({ path: 'mcp_test_screenshot.png' });
    await browser.close();
    
    if (!title.includes('Example')) {
      throw new Error(`Expected page title to contain 'Example', got: ${title}`);
    }
    
    // Verify screenshot was created
    if (!fs.existsSync('mcp_test_screenshot.png')) {
      throw new Error('Screenshot file was not created');
    }
    
    console.log('   Browser launched successfully');
    console.log('   Page navigation successful');
    console.log('   Screenshot captured successfully');
  });

  // Test 2: File operations (mock MCP functionality via Node.js)
  await runTest('File Operations', async () => {
    const testFile = 'mcp_test_file.txt';
    const testContent = 'Test content for MCP file operations';
    
    // Write file
    fs.writeFileSync(testFile, testContent);
    console.log('   File created successfully');
    
    // Read file
    const readContent = fs.readFileSync(testFile, 'utf8');
    if (readContent !== testContent) {
      throw new Error('File content does not match what was written');
    }
    console.log('   File read successfully');
    
    // Delete file
    fs.unlinkSync(testFile);
    if (fs.existsSync(testFile)) {
      throw new Error('Failed to delete test file');
    }
    console.log('   File deleted successfully');
  });

  // Test 3: Library compatibility verification
  await runTest('Library Compatibility', async () => {
    // Check if LD_LIBRARY_PATH is set correctly
    const ldLibraryPath = process.env.LD_LIBRARY_PATH || '';
    if (!ldLibraryPath.includes('lib_links')) {
      throw new Error('LD_LIBRARY_PATH does not include lib_links directory');
    }
    console.log('   LD_LIBRARY_PATH includes lib_links directory');
    
    // Check if symlinks exist
    const libLinksDir = path.join(process.cwd(), 'lib_links');
    if (!fs.existsSync(libLinksDir)) {
      throw new Error('lib_links directory does not exist');
    }
    
    const requiredLinks = [
      'libicudata.so.66',
      'libicui18n.so.66',
      'libicuuc.so.66',
      'libjpeg.so.8',
      'libwebp.so.6',
      'libffi.so.7'
    ];
    
    for (const link of requiredLinks) {
      const linkPath = path.join(libLinksDir, link);
      if (!fs.existsSync(linkPath)) {
        throw new Error(`Symlink ${link} does not exist`);
      }
    }
    console.log('   All required symlinks exist');
  });

  // For other MCP services, we can only verify indirectly
  // since we can't access them directly outside the Claude Code environment
  
  // Test 4: Context7 MCP Server (skipped in standalone test)
  skipTest('Context7 MCP Server', 'Cannot be tested outside Claude Code CLI');
  
  // Test 5: Supabase MCP Server (skipped in standalone test)
  skipTest('Supabase MCP Server', 'Cannot be tested outside Claude Code CLI');

  // Write results to a file
  const resultsFile = 'mcp_test_results.json';
  fs.writeFileSync(resultsFile, JSON.stringify(results, null, 2));
  console.log(`\n📝 Test results written to ${resultsFile}`);
  
  // Print summary
  console.log('\n📊 Test Summary');
  console.log('==================================================');
  console.log(`Passed: ${results.summary.passed}`);
  console.log(`Failed: ${results.summary.failed}`);
  console.log(`Skipped: ${results.summary.skipped}`);
  console.log('==================================================');
  
  if (results.summary.failed > 0) {
    console.log('❌ Some tests failed. See details above.');
    process.exit(1);
  } else {
    console.log('✅ All tests passed or were skipped as expected!');
  }
}

// Run the tests
runTests().catch(error => {
  console.error('Unhandled error in test suite:', error);
  process.exit(1);
});
// Verify CORS configuration between all CryoProtect services
// Usage: node verify-cors-config.js
// This script only checks CORS, not full connections

const axios = require('axios');

// Configuration defaults
const config = {
  timeout: 10000,
  frontend: process.env.FRONTEND_URL || 'https://cryoprotect.netlify.app',
  api: process.env.API_URL || 'https://cryoprotect-8030e4025428.herokuapp.com',
  rdkit: process.env.RDKIT_URL || 'https://rdkit.cryoprotect.app',
  convex: process.env.CONVEX_URL || 'https://dynamic-mink-63.convex.cloud'
};

// Format colors for terminal output
const colors = {
  reset: '\x1b[0m',
  bright: '\x1b[1m',
  dim: '\x1b[2m',
  fg: {
    red: '\x1b[31m',
    green: '\x1b[32m',
    yellow: '\x1b[33m',
    blue: '\x1b[34m',
    cyan: '\x1b[36m',
  },
  bg: {
    blue: '\x1b[44m',
    red: '\x1b[41m',
    green: '\x1b[42m',
  }
};

/**
 * Test if CORS is properly configured
 * @param {string} targetUrl - URL to test CORS against
 * @param {string} originUrl - Origin to use in the request
 * @returns {Promise<boolean>} - True if CORS is correctly configured
 */
async function testCORS(targetUrl, originUrl) {
  const name = `CORS from ${new URL(originUrl).hostname} to ${new URL(targetUrl).hostname}`;
  
  try {
    // First, test with OPTIONS request to check preflight
    const preflightResponse = await axios({
      method: 'OPTIONS',
      url: `${targetUrl}/test-cors`,
      headers: {
        'Origin': originUrl,
        'Access-Control-Request-Method': 'GET',
        'Access-Control-Request-Headers': 'X-Requested-With'
      },
      timeout: config.timeout
    });
    
    const corsHeaders = preflightResponse.headers;
    const allowOrigin = corsHeaders['access-control-allow-origin'];
    const allowMethods = corsHeaders['access-control-allow-methods'];
    
    if (!allowOrigin || !allowMethods) {
      console.log(`${colors.fg.red}❌ ${name}: Missing CORS headers${colors.reset}`);
      return false;
    }
    
    if (allowOrigin !== '*' && allowOrigin !== originUrl) {
      console.log(`${colors.fg.red}❌ ${name}: Origin not allowed (${allowOrigin})${colors.reset}`);
      return false;
    }
    
    // Now test with actual GET request
    const response = await axios.get(`${targetUrl}/test-cors`, {
      headers: {
        'Origin': originUrl
      },
      timeout: config.timeout
    });
    
    console.log(`${colors.fg.green}✅ ${name}: CORS correctly configured${colors.reset}`);
    console.log(`   Allowed Origin: ${allowOrigin}`);
    console.log(`   Allowed Methods: ${allowMethods}`);
    
    return true;
  } catch (error) {
    console.log(`${colors.fg.red}❌ ${name}: CORS test failed${colors.reset}`);
    
    if (error.response) {
      console.log(`   Status: ${error.response.status}`);
      console.log(`   Message: ${error.message}`);
    } else if (error.request) {
      console.log(`   Network error: No response received`);
      console.log(`   Is the ${targetUrl} server running and accessible?`);
    } else {
      console.log(`   Error: ${error.message}`);
    }
    
    return false;
  }
}

/**
 * Test debug CORS endpoint to show actual configuration
 * @param {string} targetUrl - URL to test
 * @returns {Promise<void>}
 */
async function checkCORSConfiguration(targetUrl) {
  const name = `CORS Config for ${new URL(targetUrl).hostname}`;
  
  try {
    const response = await axios.get(`${targetUrl}/debug/cors`, {
      timeout: config.timeout
    });
    
    console.log(`${colors.fg.green}ℹ️ ${name}: Configuration retrieved${colors.reset}`);
    
    const config = response.data;
    console.log(`   Configured Origins: ${JSON.stringify(config.configured_origins)}`);
    console.log(`   Frontend URL: ${config.frontend_url}`);
    
    if (config.rdkit_service_url) {
      console.log(`   RDKit Service URL: ${config.rdkit_service_url}`);
    }
    
    if (config.convex_url) {
      console.log(`   Convex URL: ${config.convex_url}`);
    }
    
    if (config.configured_origins.includes('https://dynamic-mink-63.convex.cloud')) {
      console.log(`   ${colors.fg.green}✅ Convex URL is properly included in CORS configuration${colors.reset}`);
    } else {
      console.log(`   ${colors.fg.red}❌ Convex URL is missing from CORS configuration${colors.reset}`);
    }
    
    return true;
  } catch (error) {
    console.log(`${colors.fg.yellow}⚠️ ${name}: Cannot retrieve CORS configuration${colors.reset}`);
    
    if (error.response) {
      console.log(`   Status: ${error.response.status}`);
      console.log(`   Message: ${error.message}`);
    } else if (error.request) {
      console.log(`   Network error: No response received`);
    } else {
      console.log(`   Error: ${error.message}`);
    }
    
    return false;
  }
}

/**
 * Run all CORS tests
 */
async function runTests() {
  console.log(`${colors.bg.blue}${colors.fg.yellow}${colors.bright} CryoProtect CORS Configuration Verification ${colors.reset}`);
  console.log(`${colors.bright}=============================================${colors.reset}`);
  console.log();
  console.log(`${colors.fg.yellow}Testing CORS between:${colors.reset}`);
  console.log(`- Frontend: ${config.frontend}`);
  console.log(`- API: ${config.api}`);
  console.log(`- RDKit Service: ${config.rdkit}`);
  console.log(`- Convex: ${config.convex}`);
  console.log();
  
  // Track overall status
  let testsPassed = 0;
  let testsFailed = 0;
  
  console.log(`${colors.fg.cyan}${colors.bright}1. CORS Configuration Information${colors.reset}`);
  console.log(`${colors.bright}-------------------------${colors.reset}`);
  
  // Check configuration on the API server
  await checkCORSConfiguration(config.api);
  
  console.log();
  console.log(`${colors.fg.cyan}${colors.bright}2. CORS Communication Tests${colors.reset}`);
  console.log(`${colors.bright}-------------------------${colors.reset}`);
  
  // Test critical paths for CORS
  console.log(`${colors.fg.yellow}Frontend to API CORS Tests:${colors.reset}`);
  if (await testCORS(config.api, config.frontend)) testsPassed++; else testsFailed++;
  
  console.log(`\n${colors.fg.yellow}Frontend to RDKit CORS Tests:${colors.reset}`);
  if (await testCORS(config.rdkit, config.frontend)) testsPassed++; else testsFailed++;
  
  console.log(`\n${colors.fg.yellow}API to RDKit CORS Tests:${colors.reset}`);
  if (await testCORS(config.rdkit, config.api)) testsPassed++; else testsFailed++;
  if (await testCORS(config.api, config.rdkit)) testsPassed++; else testsFailed++;
  
  // Print summary
  console.log();
  console.log(`${colors.bright}=============================================${colors.reset}`);
  console.log(`${colors.bright}Test Summary:${colors.reset}`);
  console.log(`${colors.fg.green}✅ Passed: ${testsPassed} tests${colors.reset}`);
  console.log(`${colors.fg.red}❌ Failed: ${testsFailed} tests${colors.reset}`);
  
  // Overall result
  if (testsFailed === 0) {
    console.log();
    console.log(`${colors.bg.green}${colors.fg.yellow}${colors.bright} ALL CORS TESTS PASSED! ${colors.reset}`);
    console.log(`${colors.fg.green}The CORS configuration is working correctly between all services.${colors.reset}`);
  } else {
    console.log();
    console.log(`${colors.bg.red}${colors.fg.yellow}${colors.bright} CORS CONFIGURATION ISSUES DETECTED ${colors.reset}`);
    console.log(`${colors.fg.yellow}See above for details on the failed tests.${colors.reset}`);
    console.log(`${colors.fg.yellow}To fix:${colors.reset}`);
    console.log(`1. Make sure the CORS configuration is deployed to all services`);
    console.log(`2. Verify environment variables are set correctly on the Heroku app`);
    console.log(`3. Run the deploy-cors-heroku.sh script to update the configuration`);
  }
}

// Run all the tests
runTests();
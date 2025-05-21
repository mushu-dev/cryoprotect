// Test connections between all CryoProtect services
// Usage: node test-connection.js

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
  underscore: '\x1b[4m',
  blink: '\x1b[5m',
  reverse: '\x1b[7m',
  hidden: '\x1b[8m',
  
  fg: {
    black: '\x1b[30m',
    red: '\x1b[31m',
    green: '\x1b[32m',
    yellow: '\x1b[33m',
    blue: '\x1b[34m',
    magenta: '\x1b[35m',
    cyan: '\x1b[36m',
    white: '\x1b[37m',
    crimson: '\x1b[38m'
  },
  
  bg: {
    black: '\x1b[40m',
    red: '\x1b[41m',
    green: '\x1b[42m',
    yellow: '\x1b[43m',
    blue: '\x1b[44m',
    magenta: '\x1b[45m',
    cyan: '\x1b[46m',
    white: '\x1b[47m',
    crimson: '\x1b[48m'
  }
};

/**
 * Test a connection to a service
 * @param {string} name - Service name
 * @param {string} url - Base URL of the service
 * @param {string} endpoint - Endpoint to test (optional)
 * @param {object} options - Axios request options (optional)
 * @returns {Promise<boolean>} - True if connection successful
 */
async function testConnection(name, url, endpoint = '', options = {}) {
  const fullUrl = `${url}${endpoint}`;
  console.log(`${colors.bright}Testing connection to ${colors.fg.yellow}${name}${colors.reset}${colors.bright}: ${colors.fg.cyan}${fullUrl}${colors.reset}`);
  
  try {
    const response = await axios.get(fullUrl, { 
      timeout: config.timeout,
      ...options
    });
    
    const status = response.status;
    const contentType = response.headers['content-type'] || '';
    
    console.log(`${colors.fg.green}✅ ${name}: Connected successfully (${status})${colors.reset}`);
    console.log(`   Content-Type: ${contentType}`);
    
    if (contentType.includes('application/json')) {
      console.log(`   Response data: ${JSON.stringify(response.data).substring(0, 100)}...`);
    }
    
    return true;
  } catch (error) {
    console.log(`${colors.fg.red}❌ ${name}: Connection failed${colors.reset}`);
    
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
      url: targetUrl + '/test-cors',
      headers: {
        'Origin': originUrl,
        'Access-Control-Request-Method': 'GET',
        'Access-Control-Request-Headers': 'X-Requested-With'
      }
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
    const response = await axios.get(targetUrl + '/test-cors', {
      headers: {
        'Origin': originUrl
      }
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
    } else {
      console.log(`   Error: ${error.message}`);
    }
    
    return false;
  }
}

/**
 * Test Netlify redirects
 * @param {string} netlifyUrl - Netlify frontend URL
 * @param {string} redirectPath - Path that should be redirected
 * @param {string} expectedContent - Content that should be in the response
 * @returns {Promise<boolean>} - True if redirect works correctly
 */
async function testRedirect(netlifyUrl, redirectPath, expectedContent) {
  const name = `Netlify Redirect: ${redirectPath}`;
  const fullUrl = `${netlifyUrl}${redirectPath}`;
  
  try {
    const response = await axios.get(fullUrl);
    const responseData = JSON.stringify(response.data).toLowerCase();
    
    if (responseData.includes(expectedContent.toLowerCase())) {
      console.log(`${colors.fg.green}✅ ${name}: Redirect working correctly${colors.reset}`);
      return true;
    } else {
      console.log(`${colors.fg.red}❌ ${name}: Got response but missing expected content${colors.reset}`);
      console.log(`   Expected to find: ${expectedContent}`);
      console.log(`   Response (partial): ${responseData.substring(0, 100)}...`);
      return false;
    }
  } catch (error) {
    console.log(`${colors.fg.red}❌ ${name}: Redirect test failed${colors.reset}`);
    
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
 * Run all connection tests
 */
async function runTests() {
  console.log(`${colors.bg.blue}${colors.fg.white}${colors.bright} CryoProtect Connection Test Tool ${colors.reset}`);
  console.log(`${colors.bright}=============================================${colors.reset}`);
  console.log();
  console.log(`${colors.fg.yellow}Testing configurations:${colors.reset}`);
  console.log(`- Frontend URL: ${config.frontend}`);
  console.log(`- API URL: ${config.api}`);
  console.log(`- RDKit Service URL: ${config.rdkit}`);
  console.log(`- Convex URL: ${config.convex}`);
  console.log();
  
  // Track overall status
  let testsPassed = 0;
  let testsFailed = 0;
  
  console.log(`${colors.fg.cyan}${colors.bright}1. Basic Connection Tests${colors.reset}`);
  console.log(`${colors.bright}-------------------------${colors.reset}`);
  
  // Test direct connections
  if (await testConnection('Frontend', config.frontend)) testsPassed++; else testsFailed++;
  if (await testConnection('API Health', config.api, '/health')) testsPassed++; else testsFailed++;
  if (await testConnection('RDKit Service Health', config.rdkit, '/health')) testsPassed++; else testsFailed++;
  // Convex requires authentication so just check if it's reachable
  if (await testConnection('Convex', config.convex)) testsPassed++; else testsFailed++;
  
  console.log();
  console.log(`${colors.fg.cyan}${colors.bright}2. CORS Configuration Tests${colors.reset}`);
  console.log(`${colors.bright}-------------------------${colors.reset}`);
  
  // Test CORS from frontend to API
  if (await testCORS(config.api, config.frontend)) testsPassed++; else testsFailed++;
  
  // Test CORS from frontend to RDKit
  if (await testCORS(config.rdkit, config.frontend)) testsPassed++; else testsFailed++;
  
  // Test CORS between services
  if (await testCORS(config.api, config.rdkit)) testsPassed++; else testsFailed++;
  if (await testCORS(config.rdkit, config.api)) testsPassed++; else testsFailed++;
  
  console.log();
  console.log(`${colors.fg.cyan}${colors.bright}3. Netlify Redirect Tests${colors.reset}`);
  console.log(`${colors.bright}-------------------------${colors.reset}`);
  
  // Test Netlify redirects
  if (await testRedirect(config.frontend, '/api/health', 'ok')) testsPassed++; else testsFailed++;
  if (await testRedirect(config.frontend, '/api/v1/health/connectivity', 'connectivity')) testsPassed++; else testsFailed++;
  
  // Test Netlify CSP headers
  console.log();
  console.log(`${colors.fg.cyan}${colors.bright}4. Content Security Policy Test${colors.reset}`);
  console.log(`${colors.bright}-------------------------${colors.reset}`);
  
  try {
    const response = await axios.get(config.frontend);
    const cspHeader = response.headers['content-security-policy'];
    
    if (cspHeader) {
      // Check if CSP allows connections to all our services
      const connectSrc = cspHeader.match(/connect-src\s+([^;]+)/)?.[1] || '';
      
      const allowsApi = connectSrc.includes(new URL(config.api).hostname) || connectSrc.includes("'self'") || connectSrc.includes('*');
      const allowsRdkit = connectSrc.includes(new URL(config.rdkit).hostname) || connectSrc.includes("'self'") || connectSrc.includes('*');
      const allowsConvex = connectSrc.includes('convex.cloud') || connectSrc.includes('*');
      
      if (allowsApi && allowsRdkit && allowsConvex) {
        console.log(`${colors.fg.green}✅ Content-Security-Policy: Allows connections to all services${colors.reset}`);
        testsPassed++;
      } else {
        console.log(`${colors.fg.red}❌ Content-Security-Policy: Missing some required domains${colors.reset}`);
        console.log(`   CSP connect-src: ${connectSrc}`);
        if (!allowsApi) console.log(`   Missing API domain: ${new URL(config.api).hostname}`);
        if (!allowsRdkit) console.log(`   Missing RDKit domain: ${new URL(config.rdkit).hostname}`);
        if (!allowsConvex) console.log(`   Missing Convex domains: convex.cloud`);
        testsFailed++;
      }
    } else {
      console.log(`${colors.fg.yellow}⚠️ No Content-Security-Policy header found${colors.reset}`);
      testsFailed++;
    }
  } catch (error) {
    console.log(`${colors.fg.red}❌ Failed to check Content-Security-Policy header${colors.reset}`);
    console.log(`   Error: ${error.message}`);
    testsFailed++;
  }
  
  // Print summary
  console.log();
  console.log(`${colors.bright}=============================================${colors.reset}`);
  console.log(`${colors.bright}Test Summary:${colors.reset}`);
  console.log(`${colors.fg.green}✅ Passed: ${testsPassed} tests${colors.reset}`);
  console.log(`${colors.fg.red}❌ Failed: ${testsFailed} tests${colors.reset}`);
  
  // Overall result
  if (testsFailed === 0) {
    console.log();
    console.log(`${colors.bg.green}${colors.fg.black}${colors.bright} ALL TESTS PASSED - Integration is working correctly! ${colors.reset}`);
  } else {
    console.log();
    console.log(`${colors.bg.red}${colors.fg.white}${colors.bright} SOME TESTS FAILED - Integration issues detected ${colors.reset}`);
    console.log(`${colors.fg.yellow}See above for details on the failed tests.${colors.reset}`);
  }
}

// Run all the tests
runTests();
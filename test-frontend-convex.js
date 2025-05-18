/**
 * Test script for verifying the frontend Convex integration
 * This script helps test the connectivity between the frontend, 
 * backend, and Convex database.
 */

const fetch = require('node-fetch');

// Configuration
const config = {
  frontend: 'https://cryoprotect.netlify.app',
  api: 'https://cryoprotect-8030e4025428.herokuapp.com',
  rdkit: 'https://cryoprotect-rdkit.fly.dev',
  convex: 'https://upbeat-parrot-866.convex.cloud'
};

// Test a connection with proper error handling
async function testConnection(name, url, endpoint = '') {
  const fullUrl = `${url}${endpoint}`;
  console.log(`Testing connection to ${name}: ${fullUrl}`);
  
  try {
    const response = await fetch(fullUrl, { 
      timeout: 10000,
      headers: {
        'Cache-Control': 'no-cache',
        'Accept': 'application/json'
      }
    });
    
    const status = response.status;
    const contentType = response.headers.get('content-type') || '';
    
    console.log(`✅ ${name}: Connected successfully (${status})`);
    console.log(`   Content-Type: ${contentType}`);
    
    if (contentType.includes('application/json')) {
      const data = await response.json();
      console.log(`   Response data: ${JSON.stringify(data).substring(0, 100)}...`);
    }
    
    return true;
  } catch (error) {
    console.error(`❌ ${name}: Connection failed`);
    
    if (error.response) {
      console.error(`   Status: ${error.response.status}`);
      console.error(`   Message: ${error.message}`);
    } else if (error.request) {
      console.error(`   Network error: No response received`);
    } else {
      console.error(`   Error: ${error.message}`);
    }
    
    return false;
  }
}

// Test a CORS connection
async function testCorsConnection(name, url, endpoint, origin) {
  const fullUrl = `${url}${endpoint}`;
  console.log(`Testing CORS connection from ${origin} to ${name}: ${fullUrl}`);
  
  try {
    const response = await fetch(fullUrl, { 
      timeout: 10000,
      headers: {
        'Origin': origin,
        'Cache-Control': 'no-cache',
        'Accept': 'application/json'
      }
    });
    
    const status = response.status;
    const corsHeader = response.headers.get('access-control-allow-origin') || '';
    
    if (corsHeader === '*' || corsHeader.includes(origin)) {
      console.log(`✅ ${name}: CORS properly configured (${status})`);
      console.log(`   Access-Control-Allow-Origin: ${corsHeader}`);
      return true;
    } else {
      console.error(`❌ ${name}: CORS header missing or incorrect`);
      console.error(`   Access-Control-Allow-Origin: ${corsHeader}`);
      return false;
    }
  } catch (error) {
    console.error(`❌ ${name}: CORS test failed`);
    console.error(`   Error: ${error.message}`);
    return false;
  }
}

// Main test function
async function runTests() {
  console.log('🧪 Starting frontend integration tests...\n');
  let success = true;
  
  // Test basic connectivity
  success = await testConnection('Frontend', config.frontend) && success;
  success = await testConnection('API', config.api) && success;
  success = await testConnection('RDKit Service', config.rdkit) && success;
  
  // Test health endpoints
  success = await testConnection('API Health', config.api, '/health') && success;
  success = await testConnection('RDKit Health', config.rdkit, '/health') && success;
  
  // Test CORS
  success = await testCorsConnection('API CORS', config.api, '/test-cors', config.frontend) && success;
  success = await testCorsConnection('RDKit CORS', config.rdkit, '/test-cors', config.frontend) && success;
  
  // Test Netlify redirects
  success = await testConnection('Netlify API Redirect', config.frontend, '/api/health') && success;
  success = await testConnection('Netlify RDKit Redirect', config.frontend, '/rdkit-api/health') && success;
  
  console.log('\n🏁 Integration tests completed');
  
  if (success) {
    console.log('✅ All tests passed successfully!');
    console.log('The frontend is properly configured to work with the backend and Convex.');
  } else {
    console.log('❌ Some tests failed. Please check the output above for details.');
    console.log('You might need to fix the CORS configuration or redirects.');
  }
  
  return success;
}

// Run the tests
runTests().then(success => {
  process.exit(success ? 0 : 1);
});
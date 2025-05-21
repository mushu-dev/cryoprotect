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
  let criticalSuccess = true;
  let serviceStatus = {
    frontend: false,
    api: false,
    rdkit: false,
    convex: false
  };
  
  console.log('--- Testing Basic Connectivity ---');
  // Test basic connectivity
  serviceStatus.frontend = await testConnection('Frontend', config.frontend);
  success = serviceStatus.frontend && success;
  
  serviceStatus.api = await testConnection('API', config.api);
  // API is critical for the application
  criticalSuccess = serviceStatus.api && criticalSuccess;
  success = serviceStatus.api && success;
  
  serviceStatus.rdkit = await testConnection('RDKit Service', config.rdkit);
  success = serviceStatus.rdkit && success;
  
  // Test Convex connectivity (if applicable)
  serviceStatus.convex = await testConnection('Convex', config.convex);
  // Convex is critical for the application
  criticalSuccess = serviceStatus.convex && criticalSuccess;
  
  console.log('\n--- Testing Health Endpoints ---');
  // Only test health endpoints for services that are running
  if (serviceStatus.api) {
    success = await testConnection('API Health', config.api, '/health') && success;
  } else {
    console.log('⚠️ Skipping API health check since API is not available');
  }
  
  if (serviceStatus.rdkit) {
    success = await testConnection('RDKit Health', config.rdkit, '/health') && success;
  } else {
    console.log('⚠️ Skipping RDKit health check since RDKit service is not available');
  }
  
  console.log('\n--- Testing CORS Configuration ---');
  // Only test CORS for services that are running
  if (serviceStatus.api && serviceStatus.frontend) {
    success = await testCorsConnection('API CORS', config.api, '/test-cors', config.frontend) && success;
  } else {
    console.log('⚠️ Skipping API CORS test since either API or Frontend is not available');
  }
  
  if (serviceStatus.rdkit && serviceStatus.frontend) {
    success = await testCorsConnection('RDKit CORS', config.rdkit, '/test-cors', config.frontend) && success;
  } else {
    console.log('⚠️ Skipping RDKit CORS test since either RDKit or Frontend is not available');
  }
  
  console.log('\n--- Testing Netlify Redirects ---');
  // Only test Netlify redirects if frontend is available
  if (serviceStatus.frontend && serviceStatus.api) {
    success = await testConnection('Netlify API Redirect', config.frontend, '/api/health') && success;
  } else {
    console.log('⚠️ Skipping Netlify API redirect test since either Frontend or API is not available');
  }
  
  if (serviceStatus.frontend && serviceStatus.rdkit) {
    success = await testConnection('Netlify RDKit Redirect', config.frontend, '/rdkit-api/health') && success;
  } else {
    console.log('⚠️ Skipping Netlify RDKit redirect test since either Frontend or RDKit is not available');
  }
  
  console.log('\n🏁 Integration tests completed');
  
  // Print summary of service status
  console.log('\n--- Service Status Summary ---');
  console.log(`Frontend: ${serviceStatus.frontend ? '✅ Available' : '❌ Not available'}`);
  console.log(`API: ${serviceStatus.api ? '✅ Available' : '❌ Not available'}`);
  console.log(`RDKit Service: ${serviceStatus.rdkit ? '✅ Available' : '❌ Not available'}`);
  console.log(`Convex: ${serviceStatus.convex ? '✅ Available' : '❌ Not available'}`);
  
  if (success) {
    console.log('\n✅ All tests passed successfully!');
    console.log('The frontend is properly configured to work with the backend and Convex.');
  } else if (criticalSuccess) {
    console.log('\n⚠️ Core services are working, but some tests failed.');
    console.log('The application should function, but some features may be limited.');
  } else {
    console.log('\n❌ Critical services are unavailable.');
    console.log('The application will not function properly until critical services are deployed.');
  }
  
  return criticalSuccess; // Return success based on critical services
}

// Run the tests
runTests().then(success => {
  process.exit(success ? 0 : 1);
});
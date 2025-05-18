/**
 * Comprehensive integration test for the CryoProtect application
 * Tests all components: Frontend, Backend, Convex, and RDKit
 */

const fetch = require('node-fetch');
const { execSync } = require('child_process');
const { setTimeout } = require('timers/promises');

// Configuration
const NETLIFY_URL = 'https://cryoprotect.netlify.app';
const HEROKU_URL = 'https://cryoprotect.herokuapp.com';
const CONVEX_URL = 'https://dynamic-mink-63.convex.cloud';
const RDKIT_URL = 'https://rdkit.cryoprotect.app';

// ANSI color codes for prettier output
const RESET = '\x1b[0m';
const RED = '\x1b[31m';
const GREEN = '\x1b[32m';
const YELLOW = '\x1b[33m';
const BLUE = '\x1b[34m';
const MAGENTA = '\x1b[35m';

async function testNetlify() {
  console.log(`${BLUE}Testing Netlify frontend...${RESET}`);
  try {
    const response = await fetch(NETLIFY_URL);
    if (!response.ok) {
      throw new Error(`HTTP error! Status: ${response.status}`);
    }
    const text = await response.text();
    if (text.includes('CryoProtect')) {
      console.log(`${GREEN}✓ Netlify frontend is up and running${RESET}`);
      return true;
    } else {
      console.log(`${RED}✗ Netlify frontend doesn't contain expected content${RESET}`);
      return false;
    }
  } catch (error) {
    console.log(`${RED}✗ Failed to reach Netlify frontend: ${error.message}${RESET}`);
    return false;
  }
}

async function testHeroku() {
  console.log(`${BLUE}Testing Heroku backend...${RESET}`);
  try {
    const response = await fetch(`${HEROKU_URL}/v1/health`);
    if (!response.ok) {
      throw new Error(`HTTP error! Status: ${response.status}`);
    }
    const data = await response.json();
    if (data.status === 'ok') {
      console.log(`${GREEN}✓ Heroku backend is up and running${RESET}`);
      return true;
    } else {
      console.log(`${RED}✗ Heroku health check failed${RESET}`);
      return false;
    }
  } catch (error) {
    console.log(`${RED}✗ Failed to reach Heroku backend: ${error.message}${RESET}`);
    return false;
  }
}

async function testConvex() {
  console.log(`${BLUE}Testing Convex connection...${RESET}`);
  try {
    // We can't directly test the Convex endpoint without authentication,
    // so we'll test it through our backend API
    const response = await fetch(`${HEROKU_URL}/v1/molecules?limit=1`);
    if (!response.ok) {
      throw new Error(`HTTP error! Status: ${response.status}`);
    }
    const data = await response.json();
    console.log(`${GREEN}✓ Successfully queried molecules through Convex: ${JSON.stringify(data).substring(0, 100)}...${RESET}`);
    return true;
  } catch (error) {
    console.log(`${RED}✗ Failed to reach Convex through backend: ${error.message}${RESET}`);
    return false;
  }
}

async function testRDKit() {
  console.log(`${BLUE}Testing RDKit service...${RESET}`);
  try {
    const response = await fetch(`${RDKIT_URL}/health`);
    if (!response.ok) {
      throw new Error(`HTTP error! Status: ${response.status}`);
    }
    const data = await response.json();
    if (data.status === 'ok') {
      console.log(`${GREEN}✓ RDKit service is up and running${RESET}`);
      return true;
    } else {
      console.log(`${RED}✗ RDKit health check failed${RESET}`);
      return false;
    }
  } catch (error) {
    console.log(`${RED}✗ Failed to reach RDKit service: ${error.message}${RESET}`);
    return false;
  }
}

async function testCORSConfiguration() {
  console.log(`${BLUE}Testing CORS configuration...${RESET}`);
  
  // Test Heroku CORS configuration
  try {
    const response = await fetch(`${HEROKU_URL}/v1/health`, {
      method: 'OPTIONS',
      headers: {
        'Origin': NETLIFY_URL,
        'Access-Control-Request-Method': 'GET',
      }
    });
    
    const allowOrigin = response.headers.get('Access-Control-Allow-Origin');
    const allowMethods = response.headers.get('Access-Control-Allow-Methods');
    
    if (allowOrigin && allowMethods) {
      console.log(`${GREEN}✓ Heroku CORS is properly configured${RESET}`);
      console.log(`  - Allow-Origin: ${allowOrigin}`);
      console.log(`  - Allow-Methods: ${allowMethods}`);
    } else {
      console.log(`${RED}✗ Heroku CORS headers are missing${RESET}`);
      return false;
    }
  } catch (error) {
    console.log(`${RED}✗ Failed to test Heroku CORS: ${error.message}${RESET}`);
    return false;
  }
  
  // Test RDKit CORS configuration
  try {
    const response = await fetch(`${RDKIT_URL}/health`, {
      method: 'OPTIONS',
      headers: {
        'Origin': NETLIFY_URL,
        'Access-Control-Request-Method': 'GET',
      }
    });
    
    const allowOrigin = response.headers.get('Access-Control-Allow-Origin');
    const allowMethods = response.headers.get('Access-Control-Allow-Methods');
    
    if (allowOrigin && allowMethods) {
      console.log(`${GREEN}✓ RDKit CORS is properly configured${RESET}`);
      console.log(`  - Allow-Origin: ${allowOrigin}`);
      console.log(`  - Allow-Methods: ${allowMethods}`);
    } else {
      console.log(`${YELLOW}⚠ RDKit CORS headers are missing${RESET}`);
      console.log(`  - This may need to be configured`);
    }
  } catch (error) {
    console.log(`${YELLOW}⚠ Failed to test RDKit CORS: ${error.message}${RESET}`);
  }
  
  return true;
}

async function testEndToEndFlow() {
  console.log(`${BLUE}Testing end-to-end flow...${RESET}`);
  
  // Test creating a new molecule
  try {
    // First, create a test molecule
    const createResponse = await fetch(`${HEROKU_URL}/v1/molecules`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        name: `Test Molecule ${Date.now()}`,
        formula: 'C2H5OH',
        smiles: 'CCO',
        status: 'test'
      })
    });
    
    if (!createResponse.ok) {
      throw new Error(`HTTP error! Status: ${createResponse.status}`);
    }
    
    const createData = await createResponse.json();
    const moleculeId = createData[0].id;
    
    console.log(`${GREEN}✓ Successfully created a test molecule with ID: ${moleculeId}${RESET}`);
    
    // Wait a moment for the data to propagate
    await setTimeout(2000);
    
    // Now fetch the molecule
    const fetchResponse = await fetch(`${HEROKU_URL}/v1/molecules/${moleculeId}`);
    
    if (!fetchResponse.ok) {
      throw new Error(`HTTP error! Status: ${fetchResponse.status}`);
    }
    
    const fetchData = await fetchResponse.json();
    
    if (fetchData && fetchData.id === moleculeId) {
      console.log(`${GREEN}✓ Successfully fetched the test molecule${RESET}`);
    } else {
      console.log(`${RED}✗ Failed to fetch the test molecule${RESET}`);
      return false;
    }
    
    // Delete the test molecule
    const deleteResponse = await fetch(`${HEROKU_URL}/v1/molecules/${moleculeId}`, {
      method: 'DELETE'
    });
    
    if (!deleteResponse.ok) {
      console.log(`${YELLOW}⚠ Failed to delete the test molecule: ${deleteResponse.status}${RESET}`);
    } else {
      console.log(`${GREEN}✓ Successfully deleted the test molecule${RESET}`);
    }
    
    return true;
  } catch (error) {
    console.log(`${RED}✗ End-to-end test failed: ${error.message}${RESET}`);
    return false;
  }
}

async function main() {
  console.log(`${MAGENTA}==========================================${RESET}`);
  console.log(`${MAGENTA}  CryoProtect Integration Test Suite  ${RESET}`);
  console.log(`${MAGENTA}==========================================${RESET}`);
  
  let success = true;
  
  // Run all tests
  success = await testNetlify() && success;
  success = await testHeroku() && success;
  success = await testConvex() && success;
  success = await testRDKit() && success;
  success = await testCORSConfiguration() && success;
  success = await testEndToEndFlow() && success;
  
  console.log(`${MAGENTA}==========================================${RESET}`);
  if (success) {
    console.log(`${GREEN}All tests passed successfully!${RESET}`);
  } else {
    console.log(`${RED}Some tests failed. Please check the logs above.${RESET}`);
  }
  console.log(`${MAGENTA}==========================================${RESET}`);
  
  // Return exit code based on success/failure
  process.exit(success ? 0 : 1);
}

main().catch(error => {
  console.error(`${RED}Error running tests: ${error}${RESET}`);
  process.exit(1);
});
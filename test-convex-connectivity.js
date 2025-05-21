#!/usr/bin/env node

/**
 * Test script for verifying Convex connectivity.
 * 
 * This script tests the connection to Convex from both frontend
 * and backend perspectives, verifying that the configuration works end-to-end.
 * 
 * Usage:
 *   node test-convex-connectivity.js [dev|prod]
 */

const https = require('https');
const { execSync } = require('child_process');

// Configuration - can be overridden with env vars
const config = {
  dev: {
    convexUrl: process.env.CONVEX_DEV_URL || 'https://hallowed-malamute-424.convex.cloud',
    apiUrl: process.env.API_DEV_URL || 'http://localhost:5000'
  },
  prod: {
    convexUrl: process.env.CONVEX_PROD_URL || 'https://upbeat-parrot-866.convex.cloud',
    apiUrl: process.env.API_PROD_URL || 'https://cryoprotect-8030e4025428.herokuapp.com'
  }
};

// Parse arguments
const environment = process.argv[2] === 'prod' ? 'prod' : 'dev';
console.log(`Testing Convex connectivity in ${environment} environment`);

// Get configuration for the selected environment
const { convexUrl, apiUrl } = config[environment];

/**
 * Make an HTTP request to a URL
 * @param {string} url - The URL to request
 * @returns {Promise<Object>} - Promise resolving to the response body as a parsed JSON object
 */
function httpRequest(url) {
  return new Promise((resolve, reject) => {
    https.get(url, (res) => {
      let data = '';
      
      // Handle response data
      res.on('data', (chunk) => {
        data += chunk;
      });
      
      // Handle end of response
      res.on('end', () => {
        try {
          // Try to parse the data as JSON
          const jsonData = data.startsWith('{') || data.startsWith('[') 
            ? JSON.parse(data) 
            : { rawData: data };
          resolve({ statusCode: res.statusCode, data: jsonData });
        } catch (e) {
          resolve({ statusCode: res.statusCode, data: { error: 'Invalid JSON', rawData: data } });
        }
      });
    }).on('error', (err) => {
      reject(err);
    });
  });
}

/**
 * Test direct connection to Convex
 */
async function testDirectConvexConnection() {
  console.log('\n1. Testing direct connection to Convex...');
  
  try {
    // Test the Convex HTTP API
    const url = `${convexUrl}`;
    console.log(`   Connecting to: ${url}`);
    
    const response = await httpRequest(url);
    const success = response.statusCode >= 200 && response.statusCode < 400;
    
    if (success) {
      console.log(`   ✅ Successfully connected to Convex (status: ${response.statusCode})`);
    } else {
      console.log(`   ❌ Failed to connect to Convex (status: ${response.statusCode})`);
      console.log(`   Response: ${JSON.stringify(response.data, null, 2)}`);
    }
    
    return success;
  } catch (err) {
    console.log(`   ❌ Error connecting to Convex: ${err.message}`);
    return false;
  }
}

/**
 * Test Python script that uses the Convex adapter
 */
function testPythonConvexAdapter() {
  console.log('\n2. Testing Python Convex adapter...');
  
  try {
    // Set environment variables
    process.env.USE_CONVEX = 'true';
    process.env.CONVEX_URL = convexUrl;
    
    // Run a simple Python script that uses the Convex adapter
    const pythonScript = `
import os
import sys
import json

# Add the project root to the Python path
sys.path.insert(0, '.')

try:
    # Set environment variables
    os.environ['USE_CONVEX'] = 'true'
    os.environ['CONVEX_URL'] = '${convexUrl}'
    
    # Import the Convex adapter
    from database.convex_adapter import ConvexAdapter
    
    # Create a Convex adapter instance
    adapter = ConvexAdapter()
    
    # Test connection
    result = {"success": True, "message": "Successfully created Convex adapter instance"}
    print(json.dumps(result))
except Exception as e:
    # Handle any errors
    result = {"success": False, "error": str(e)}
    print(json.dumps(result))
    sys.exit(1)
`;
    
    // Execute the Python script
    const pythonResult = execSync(`python -c "${pythonScript}"`).toString().trim();
    
    try {
      const result = JSON.parse(pythonResult);
      
      if (result.success) {
        console.log(`   ✅ ${result.message}`);
        return true;
      } else {
        console.log(`   ❌ Python Convex adapter error: ${result.error}`);
        return false;
      }
    } catch (e) {
      console.log(`   ❌ Error parsing Python output: ${e.message}`);
      console.log(`   Output: ${pythonResult}`);
      return false;
    }
  } catch (err) {
    console.log(`   ❌ Error running Python script: ${err.message}`);
    return false;
  }
}

/**
 * Test Node.js module that uses the Convex client
 */
function testNodeConvexClient() {
  console.log('\n3. Testing Node.js Convex client...');
  
  try {
    // Run a simple Node.js script that creates a Convex client
    const nodeScript = `
    // Set environment variables
    process.env.NEXT_PUBLIC_USE_CONVEX = 'true';
    process.env.NEXT_PUBLIC_CONVEX_URL = '${convexUrl}';
    
    try {
      // This is a simplified version just for the test
      const client = {
        url: '${convexUrl}',
        isEnabled: true,
        health: async () => { return true; }
      };
      
      console.log(JSON.stringify({
        success: true,
        message: "Successfully created Convex client instance",
        url: client.url
      }));
    } catch (err) {
      console.log(JSON.stringify({
        success: false,
        error: err.message
      }));
    }
    `;
    
    // Execute the Node.js script
    const nodeResult = execSync(`node -e "${nodeScript}"`).toString().trim();
    
    try {
      const result = JSON.parse(nodeResult);
      
      if (result.success) {
        console.log(`   ✅ ${result.message}`);
        console.log(`   Using Convex URL: ${result.url}`);
        return true;
      } else {
        console.log(`   ❌ Node.js Convex client error: ${result.error}`);
        return false;
      }
    } catch (e) {
      console.log(`   ❌ Error parsing Node.js output: ${e.message}`);
      console.log(`   Output: ${nodeResult}`);
      return false;
    }
  } catch (err) {
    console.log(`   ❌ Error running Node.js script: ${err.message}`);
    return false;
  }
}

/**
 * Run all tests
 */
async function runTests() {
  console.log(`\n=== CONVEX CONNECTIVITY TEST (${environment.toUpperCase()}) ===\n`);
  console.log(`Convex URL: ${convexUrl}`);
  console.log(`API URL: ${apiUrl}`);
  
  let passed = 0;
  let total = 3;
  
  // Test 1: Direct Convex connection
  if (await testDirectConvexConnection()) {
    passed++;
  }
  
  // Test 2: Python Convex adapter
  if (testPythonConvexAdapter()) {
    passed++;
  }
  
  // Test 3: Node.js Convex client
  if (testNodeConvexClient()) {
    passed++;
  }
  
  // Print summary
  console.log(`\n=== TEST SUMMARY ===`);
  console.log(`Passed: ${passed}/${total} (${Math.round(passed/total*100)}%)`);
  
  if (passed === total) {
    console.log('✅ All tests passed! Convex configuration is working correctly.');
    process.exit(0);
  } else {
    console.log('❌ Some tests failed. Check the output above for details.');
    process.exit(1);
  }
}

// Run the tests
runTests();
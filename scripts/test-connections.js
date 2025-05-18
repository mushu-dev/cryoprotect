#!/usr/bin/env node

/**
 * Connection Testing Script for CryoProtect Services
 * 
 * This script tests connectivity between all services in the CryoProtect architecture:
 * - Netlify frontend
 * - Heroku backend API
 * - RDKit fly.io service
 * - Convex database
 * 
 * Usage: node test-connections.js
 */

const axios = require('axios');
const { execSync } = require('child_process');

// Configuration
const config = {
  frontend: 'https://cryoprotect.netlify.app',
  api: 'https://api.cryoprotect.app', 
  rdkit: 'https://rdkit.cryoprotect.app',
  convex: 'https://dynamic-mink-63.convex.cloud',
  timeout: 10000
};

async function testConnection(name, url, endpoint = '', options = {}) {
  const fullUrl = `${url}${endpoint}`;
  console.log(`Testing connection to ${name}: ${fullUrl}`);
  
  try {
    const response = await axios.get(fullUrl, { 
      timeout: config.timeout,
      ...options
    });
    
    const status = response.status;
    const contentType = response.headers['content-type'] || '';
    
    console.log(`✅ ${name}: Connected successfully (${status})`);
    console.log(`   Content-Type: ${contentType}`);
    
    if (contentType.includes('application/json')) {
      console.log(`   Response data: ${JSON.stringify(response.data).substring(0, 100)}...`);
    }
    
    return true;
  } catch (error) {
    console.error(`❌ ${name}: Connection failed`);
    
    if (error.response) {
      console.error(`   Status: ${error.response.status}`);
      console.error(`   Message: ${error.message}`);
    } else if (error.request) {
      console.error(`   Network error: ${error.message}`);
    } else {
      console.error(`   Error: ${error.message}`);
    }
    
    return false;
  }
}

async function testCORS(name, url, origin) {
  console.log(`Testing CORS from ${origin} to ${name}: ${url}`);
  
  try {
    const response = await axios.get(`${url}/test-cors`, { 
      timeout: config.timeout,
      headers: {
        'Origin': origin
      }
    });
    
    console.log(`✅ ${name}: CORS test successful from ${origin}`);
    return true;
  } catch (error) {
    console.error(`❌ ${name}: CORS test failed from ${origin}`);
    
    if (error.response) {
      console.error(`   Status: ${error.response.status}`);
      if (error.response.headers['access-control-allow-origin']) {
        console.error(`   Allowed origins: ${error.response.headers['access-control-allow-origin']}`);
      }
    } else {
      console.error(`   Error: ${error.message}`);
    }
    
    return false;
  }
}

function testEnvVariables() {
  console.log('\nChecking Netlify environment variables:');
  
  try {
    // Use curl with netlify CLI to get environment variables
    // Note: This requires netlify CLI to be installed and authenticated
    console.log('Running: netlify env:list');
    const result = execSync('netlify env:list --json').toString();
    const variables = JSON.parse(result);
    
    const requiredVars = [
      'NEXT_PUBLIC_API_URL',
      'NEXT_PUBLIC_RDKIT_API_URL',
      'NEXT_PUBLIC_CONVEX_URL',
      'NEXT_PUBLIC_USE_CONVEX'
    ];
    
    let allFound = true;
    requiredVars.forEach(varName => {
      const found = variables.some(v => v.key === varName);
      if (found) {
        console.log(`✅ Found ${varName}`);
      } else {
        console.log(`❌ Missing ${varName}`);
        allFound = false;
      }
    });
    
    return allFound;
  } catch (error) {
    console.error(`❌ Error checking environment variables: ${error.message}`);
    console.log('   Note: This test requires netlify CLI to be installed and authenticated');
    return false;
  }
}

async function testAllConnections() {
  console.log('=== CryoProtect Connection Testing ===\n');
  
  // Test direct connections to all services
  console.log('\n--- Testing Basic Connectivity ---');
  const frontendOk = await testConnection('Frontend', config.frontend);
  const apiOk = await testConnection('API', config.api, '/health');
  const rdkitOk = await testConnection('RDKit', config.rdkit, '/health');
  const convexOk = await testConnection('Convex', config.convex);
  
  // Test CORS if APIs are reachable
  console.log('\n--- Testing CORS Configuration ---');
  if (apiOk) {
    await testCORS('API', config.api, config.frontend);
  }
  
  if (rdkitOk) {
    await testCORS('RDKit', config.rdkit, config.frontend);
  }
  
  // Check Netlify redirects by testing through the frontend
  console.log('\n--- Testing Netlify Redirects ---');
  const apiRedirectOk = await testConnection('API via Netlify', config.frontend, '/api/health');
  const rdkitRedirectOk = await testConnection('RDKit via Netlify', config.frontend, '/rdkit-api/health');
  
  // Test Netlify environment variables
  console.log('\n--- Testing Netlify Environment Variables ---');
  const envVarsOk = testEnvVariables();
  
  // Summary
  console.log('\n=== Connection Test Summary ===');
  console.log(`Frontend: ${frontendOk ? '✅ Connected' : '❌ Failed'}`);
  console.log(`API: ${apiOk ? '✅ Connected' : '❌ Failed'}`);
  console.log(`RDKit: ${rdkitOk ? '✅ Connected' : '❌ Failed'}`);
  console.log(`Convex: ${convexOk ? '✅ Connected' : '❌ Failed'}`);
  console.log(`API via Netlify: ${apiRedirectOk ? '✅ Connected' : '❌ Failed'}`);
  console.log(`RDKit via Netlify: ${rdkitRedirectOk ? '✅ Connected' : '❌ Failed'}`);
  console.log(`Environment Variables: ${envVarsOk ? '✅ Configured' : '❌ Issues found'}`);
  
  const allOk = frontendOk && apiOk && rdkitOk && convexOk && apiRedirectOk && rdkitRedirectOk && envVarsOk;
  
  console.log(`\nOverall status: ${allOk ? '✅ All connections successful' : '❌ Some connections failed'}`);
  
  if (!allOk) {
    console.log('\nSuggested fixes:');
    if (!frontendOk) console.log('- Check if Netlify site is deployed and accessible');
    if (!apiOk) console.log('- Verify Heroku API is running and the DNS record is correct');
    if (!rdkitOk) console.log('- Ensure RDKit service on fly.io is running and the DNS record is correct');
    if (!convexOk) console.log('- Check Convex deployment URL and ensure it\'s active');
    if (!apiRedirectOk) console.log('- Verify Netlify redirect rules for /api/* are correctly configured');
    if (!rdkitRedirectOk) console.log('- Verify Netlify redirect rules for /rdkit-api/* are correctly configured');
    if (!envVarsOk) console.log('- Update Netlify environment variables as per the integration plan');
  }
}

// Run all tests
testAllConnections();
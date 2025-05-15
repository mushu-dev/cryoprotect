#!/usr/bin/env node
/**
 * Deployment Verification Script
 * 
 * This script performs basic checks on a deployed frontend to verify it's working correctly.
 * It checks for:
 * 1. HTTP status of the main page
 * 2. Presence of expected HTML elements
 * 3. Basic functionality checks
 * 
 * Usage:
 *   node scripts/verify-deployment.js --url=https://example.com [--timeout=30000]
 */

const https = require('https');
const http = require('http');
const { URL } = require('url');
const { spawn } = require('child_process');

// Parse command line arguments
const args = process.argv.slice(2).reduce((acc, arg) => {
  const [key, value] = arg.replace(/^--/, '').split('=');
  acc[key] = value;
  return acc;
}, {});

// Get deployment URL to verify
const deployUrl = args.url;
if (!deployUrl) {
  console.error('Error: No deployment URL provided. Use --url=https://example.com');
  process.exit(1);
}

// Timeout for each check (default: 30 seconds)
const timeout = parseInt(args.timeout) || 30000;

// Validation results
const results = {
  statusCheck: { status: 'pending', message: '', time: 0 },
  contentCheck: { status: 'pending', message: '', time: 0 },
  errorCheck: { status: 'pending', message: '', time: 0 },
};

let exitCode = 0;

/**
 * Fetch a URL and return the response
 */
function fetchUrl(url) {
  return new Promise((resolve, reject) => {
    const requestUrl = new URL(url);
    const client = requestUrl.protocol === 'https:' ? https : http;
    
    const req = client.get(url, (res) => {
      let data = '';
      
      // A chunk of data has been received
      res.on('data', (chunk) => {
        data += chunk;
      });
      
      // The whole response has been received
      res.on('end', () => {
        resolve({
          statusCode: res.statusCode,
          headers: res.headers,
          body: data,
        });
      });
    });
    
    // Set a timeout
    req.setTimeout(timeout, () => {
      req.abort();
      reject(new Error(`Request timed out after ${timeout}ms`));
    });
    
    // Handle errors
    req.on('error', (error) => {
      reject(error);
    });
  });
}

/**
 * Check HTTP status
 */
async function checkStatus() {
  const startTime = Date.now();
  
  try {
    console.log(`Checking HTTP status of ${deployUrl}...`);
    const response = await fetchUrl(deployUrl);
    const endTime = Date.now();
    const duration = endTime - startTime;
    
    if (response.statusCode >= 200 && response.statusCode < 400) {
      results.statusCheck = {
        status: 'success',
        message: `Received status code ${response.statusCode} in ${duration}ms`,
        time: duration,
      };
      return true;
    } else {
      results.statusCheck = {
        status: 'failure',
        message: `Received unexpected status code ${response.statusCode} in ${duration}ms`,
        time: duration,
      };
      return false;
    }
  } catch (error) {
    const endTime = Date.now();
    const duration = endTime - startTime;
    
    results.statusCheck = {
      status: 'failure',
      message: `Error checking status: ${error.message} (${duration}ms)`,
      time: duration,
    };
    return false;
  }
}

/**
 * Check content for expected elements
 */
async function checkContent() {
  const startTime = Date.now();
  
  try {
    console.log(`Checking content of ${deployUrl}...`);
    const response = await fetchUrl(deployUrl);
    const endTime = Date.now();
    const duration = endTime - startTime;
    
    // Check for expected HTML elements
    const expectedElements = [
      '<title>CryoProtect</title>',  // Title
      '<header',                    // Header
      '<main',                      // Main content
      '<footer',                    // Footer
    ];
    
    const missingElements = expectedElements.filter(element => 
      !response.body.includes(element)
    );
    
    if (missingElements.length === 0) {
      results.contentCheck = {
        status: 'success',
        message: `All expected elements found in ${duration}ms`,
        time: duration,
      };
      return true;
    } else {
      results.contentCheck = {
        status: 'failure',
        message: `Missing elements: ${missingElements.join(', ')} (${duration}ms)`,
        time: duration,
      };
      return false;
    }
  } catch (error) {
    const endTime = Date.now();
    const duration = endTime - startTime;
    
    results.contentCheck = {
      status: 'failure',
      message: `Error checking content: ${error.message} (${duration}ms)`,
      time: duration,
    };
    return false;
  }
}

/**
 * Check for client-side errors by running Lighthouse
 * This is a basic implementation - in the future, we could use Puppeteer
 * or other tools for more comprehensive checks
 */
async function checkErrors() {
  // For now, we'll implement a simple check that the page doesn't 
  // contain obvious error messages
  const startTime = Date.now();
  
  try {
    console.log(`Checking for errors in ${deployUrl}...`);
    const response = await fetchUrl(deployUrl);
    const endTime = Date.now();
    const duration = endTime - startTime;
    
    // Check for obvious error messages
    const errorPatterns = [
      'Internal Server Error',
      'Cannot read property',
      'is not defined',
      'is not a function',
      'Error: Failed to fetch',
    ];
    
    const foundErrors = errorPatterns.filter(pattern => 
      response.body.includes(pattern)
    );
    
    if (foundErrors.length === 0) {
      results.errorCheck = {
        status: 'success',
        message: `No obvious errors found in ${duration}ms`,
        time: duration,
      };
      return true;
    } else {
      results.errorCheck = {
        status: 'failure',
        message: `Found error patterns: ${foundErrors.join(', ')} (${duration}ms)`,
        time: duration,
      };
      return false;
    }
  } catch (error) {
    const endTime = Date.now();
    const duration = endTime - startTime;
    
    results.errorCheck = {
      status: 'failure',
      message: `Error checking for errors: ${error.message} (${duration}ms)`,
      time: duration,
    };
    return false;
  }
}

/**
 * Run all checks and print a summary
 */
async function runChecks() {
  console.log(`Verifying deployment at: ${deployUrl}`);
  console.log('-'.repeat(50));
  
  const statusOk = await checkStatus();
  
  // Only continue with other checks if status check passed
  if (statusOk) {
    await Promise.all([
      checkContent(),
      checkErrors(),
    ]);
  }
  
  // Print summary
  console.log('\nDeployment Verification Results:');
  console.log('='.repeat(50));
  
  for (const [checkName, result] of Object.entries(results)) {
    const icon = result.status === 'success' ? '✅' : result.status === 'failure' ? '❌' : '⏳';
    console.log(`${icon} ${checkName}: ${result.status}`);
    console.log(`   ${result.message}`);
  }
  
  // Determine overall status
  const allPassed = Object.values(results).every(result => result.status === 'success');
  
  console.log('='.repeat(50));
  if (allPassed) {
    console.log('✅ All checks passed! Deployment verification successful.');
    exitCode = 0;
  } else {
    console.log('❌ Some checks failed! Deployment verification unsuccessful.');
    exitCode = 1;
  }
  
  // If this were a real implementation, we would also:
  // 1. Write results to a file for CI/CD pipeline to read
  // 2. Take screenshots for comparison
  // 3. Run more thorough tests
  process.exit(exitCode);
}

// Run the verification
runChecks();
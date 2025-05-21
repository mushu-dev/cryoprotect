#!/usr/bin/env node

// Simple script to check if the experimental data enhancement UI is deployed and working

const https = require('https');
const { exec } = require('child_process');

const BASE_URL = 'https://cryoprotect.app';
const PAGES_TO_CHECK = [
  '/',                 // Homepage
  '/experiments',      // Experiments list
  '/protocols',        // Protocols list
  '/experiments/1',    // Experiment detail (dynamic route)
  '/protocols/1',      // Protocol detail (dynamic route)
  '/protocols/create'  // Protocol create page
];

// Function to check if a URL is accessible
function checkUrl(url) {
  return new Promise((resolve, reject) => {
    https.get(url, (res) => {
      const { statusCode } = res;
      resolve({ url, statusCode, success: statusCode === 200 });
    }).on('error', (err) => {
      resolve({ url, success: false, error: err.message });
    });
  });
}

// Function to print a colored message
function colorLog(message, success) {
  const color = success ? '\x1b[32m' : '\x1b[31m'; // Green or Red
  const reset = '\x1b[0m';
  console.log(`${color}${message}${reset}`);
}

// Main function to check all pages
async function checkDeployment() {
  console.log('\nChecking Experimental Data Enhancement UI Deployment');
  console.log('===================================================\n');
  
  let allSuccess = true;
  
  for (const page of PAGES_TO_CHECK) {
    const url = `${BASE_URL}${page}`;
    const result = await checkUrl(url);
    
    if (result.success) {
      colorLog(`✅ ${url} - Status: ${result.statusCode}`, true);
    } else {
      allSuccess = false;
      const errorMsg = result.error || `Status: ${result.statusCode}`;
      colorLog(`❌ ${url} - ${errorMsg}`, false);
    }
  }
  
  console.log('\n===================================================');
  
  if (allSuccess) {
    colorLog('\n✅ SUCCESS: All pages are accessible!', true);
    colorLog('The experimental data enhancement UI is fully deployed and functional.', true);
  } else {
    colorLog('\n❌ WARNING: Some pages are not accessible.', false);
    colorLog('The deployment may need additional fixes or may still be in progress.', false);
    colorLog('Run the deploy-experimental-ui-fix.sh script to fix any issues.', false);
  }
  
  // Open the site in the default browser if available
  if (allSuccess) {
    console.log('\nOpening the site in your browser...');
    try {
      const command = process.platform === 'win32' ? 'start' :
                      process.platform === 'darwin' ? 'open' : 'xdg-open';
      exec(`${command} ${BASE_URL}`);
    } catch (error) {
      console.log('Could not open browser automatically. Please visit:');
      console.log(BASE_URL);
    }
  }
}

// Run the check
checkDeployment();
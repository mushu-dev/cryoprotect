#!/usr/bin/env node
/**
 * Route Verification Script
 * 
 * This script checks if all critical routes are accessible in the deployed application.
 * 
 * Usage:
 *   node scripts/verify-routes.js --url=https://cryoprotect.app
 */

const fetch = require('node-fetch');

// Parse command line arguments
const args = process.argv.slice(2).reduce((acc, arg) => {
  const [key, value] = arg.replace(/^--/, '').split('=');
  acc[key] = value;
  return acc;
}, {});

// Get deployment URL to verify (default to production)
const BASE_URL = args.url || 'https://cryoprotect.app';

// Critical routes that should be accessible
const ROUTES = [
  '/',
  '/molecules',
  '/mixtures',
  '/experiments',
  '/protocols',
  '/properties',
  '/api/health'
];

async function verifyRoutes() {
  console.log(`Starting route verification for ${BASE_URL}...`);
  console.log('-'.repeat(50));
  
  let successCount = 0;
  let failureCount = 0;
  
  for (const route of ROUTES) {
    const url = `${BASE_URL}${route}`;
    try {
      console.log(`Checking ${url}`);
      const response = await fetch(url);
      const status = response.status;
      
      if (status === 200) {
        console.log(`✅ ${route} - OK (${status})`);
        successCount++;
      } else {
        console.error(`❌ ${route} - Failed (${status})`);
        failureCount++;
      }
    } catch (error) {
      console.error(`❌ ${route} - Error: ${error.message}`);
      failureCount++;
    }
  }
  
  console.log('-'.repeat(50));
  console.log(`Route verification complete: ${successCount} passed, ${failureCount} failed`);
  
  if (failureCount > 0) {
    console.log('\n❌ Some routes are inaccessible. Review the Netlify deployment configuration.');
    process.exit(1);
  } else {
    console.log('\n✅ All routes are accessible!');
    process.exit(0);
  }
}

verifyRoutes();
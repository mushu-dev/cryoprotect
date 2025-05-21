// Test the Convex integration directly
// Usage: node test-convex-integration.js

const axios = require('axios');

// Configuration
const config = {
  timeout: 10000,
  api: process.env.API_URL || 'https://cryoprotect-8030e4025428.herokuapp.com',
  convex: process.env.CONVEX_URL || 'https://dynamic-mink-63.convex.cloud'
};

// Format colors for terminal output
const colors = {
  reset: '\x1b[0m',
  bright: '\x1b[1m',
  fg: {
    red: '\x1b[31m',
    green: '\x1b[32m',
    yellow: '\x1b[33m',
    blue: '\x1b[34m',
    cyan: '\x1b[36m',
  },
  bg: {
    blue: '\x1b[44m',
    green: '\x1b[42m',
    red: '\x1b[41m',
  }
};

/**
 * Test the database integration
 */
async function testDatabaseIntegration() {
  console.log(`${colors.bg.blue}${colors.fg.yellow}${colors.bright} CryoProtect Convex Integration Test ${colors.reset}`);
  console.log(`${colors.bright}=============================================${colors.reset}`);
  console.log();
  console.log(`${colors.fg.yellow}Testing configurations:${colors.reset}`);
  console.log(`- API URL: ${config.api}`);
  console.log(`- Convex URL: ${config.convex}`);
  console.log();
  
  // Test database status
  try {
    console.log(`${colors.fg.cyan}${colors.bright}1. Database Status${colors.reset}`);
    console.log(`${colors.bright}-------------------------${colors.reset}`);

    const response = await axios.get(`${config.api}/api/database/status`, {
      timeout: config.timeout
    });
    
    console.log(`${colors.fg.green}✅ Database Status: ${response.data.database_type}${colors.reset}`);
    console.log(`   URL: ${response.data.details.url}`);
    console.log(`   Key: ${response.data.details.deployment_key}`);
    
    if (response.data.database_type === 'Convex') {
      console.log(`${colors.fg.green}✅ Backend is using Convex as the database!${colors.reset}`);
    } else {
      console.log(`${colors.fg.red}❌ Backend is using ${response.data.database_type} as the database, not Convex!${colors.reset}`);
    }
  } catch (error) {
    console.log(`${colors.fg.red}❌ Failed to check database status${colors.reset}`);
    console.log(`   Error: ${error.message}`);
  }
  
  // Test molecules API (which uses the database)
  try {
    console.log();
    console.log(`${colors.fg.cyan}${colors.bright}2. Molecules API${colors.reset}`);
    console.log(`${colors.bright}-------------------------${colors.reset}`);

    const response = await axios.get(`${config.api}/api/molecules?limit=3`, {
      timeout: config.timeout
    });
    
    if (response.data.status === 'success' && Array.isArray(response.data.data)) {
      console.log(`${colors.fg.green}✅ Molecules API: Data retrieved successfully${colors.reset}`);
      console.log(`   Number of molecules: ${response.data.data.length}`);
      
      if (response.data.data.length > 0) {
        const firstMolecule = response.data.data[0];
        console.log(`   First molecule: ${firstMolecule.name || 'Unnamed'} (ID: ${firstMolecule.id})`);
      }
    } else {
      console.log(`${colors.fg.red}❌ Molecules API: Failed to retrieve data${colors.reset}`);
      console.log(`   Response: ${JSON.stringify(response.data)}`);
    }
  } catch (error) {
    console.log(`${colors.fg.red}❌ Failed to check molecules API${colors.reset}`);
    console.log(`   Error: ${error.message}`);
    if (error.response) {
      console.log(`   Status: ${error.response.status}`);
      console.log(`   Data: ${JSON.stringify(error.response.data)}`);
    }
  }
  
  // Test direct Convex connectivity
  try {
    console.log();
    console.log(`${colors.fg.cyan}${colors.bright}3. Direct Convex Connection${colors.reset}`);
    console.log(`${colors.bright}-------------------------${colors.reset}`);

    const response = await axios.get(config.convex, {
      timeout: config.timeout
    });
    
    console.log(`${colors.fg.green}✅ Convex Connection: Connected successfully (${response.status})${colors.reset}`);
    console.log(`   Content Type: ${response.headers['content-type']}`);
  } catch (error) {
    // It's okay if this fails, since direct connection to Convex might be restricted
    console.log(`${colors.fg.yellow}⚠️ Direct Convex Connection: Not accessible directly${colors.reset}`);
    console.log(`   This is expected since Convex might restrict direct access`);
    console.log(`   Error: ${error.message}`);
  }
}

// Run the test
testDatabaseIntegration();
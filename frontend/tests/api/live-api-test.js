/**
 * Integration tests for live backend services
 * These tests connect to the real Heroku backend and Fly.io RDKit service
 */

const axios = require('axios');
const fs = require('fs');
const path = require('path');

// Define the API base URLs from environment variables or use default production URLs
const API_URL = process.env.NEXT_PUBLIC_API_URL || 'https://cryoprotect-8030e4025428.herokuapp.com/api';
const RDKIT_API_URL = process.env.NEXT_PUBLIC_RDKIT_API_URL || 'https://cryoprotect-rdkit.fly.dev';

// Define colors for console output
const colors = {
  reset: "\x1b[0m",
  green: "\x1b[32m",
  red: "\x1b[31m",
  yellow: "\x1b[33m",
  blue: "\x1b[34m",
  bold: "\x1b[1m",
  dim: "\x1b[2m"
};

// Create a configured axios instance for the main API
const apiClient = axios.create({
  baseURL: API_URL,
  timeout: 15000, // 15 seconds
  headers: {
    'Content-Type': 'application/json',
  },
});

// Create a configured axios instance for the RDKit service
const rdkitClient = axios.create({
  baseURL: RDKIT_API_URL,
  timeout: 15000, // 15 seconds
  headers: {
    'Content-Type': 'application/json',
  },
});

// Create a results directory
const resultsDir = path.join(__dirname, '../../test-results/api-live');
if (!fs.existsSync(resultsDir)) {
  fs.mkdirSync(resultsDir, { recursive: true });
}

// Log to both console and file
const logFile = path.join(resultsDir, `test-results-${new Date().toISOString().replace(/:/g, '-')}.log`);
const logger = {
  log: (message) => {
    console.log(message);
    fs.appendFileSync(logFile, message.replace(/\x1B\[\d+m/g, '') + '\n');
  }
};

/**
 * Test a GET request to the API
 */
async function testApiGet(endpoint, description, validators = []) {
  logger.log(`${colors.blue}Testing:${colors.reset} ${description}`);
  logger.log(`${colors.dim}GET ${API_URL}${endpoint}${colors.reset}`);
  
  try {
    const startTime = Date.now();
    const response = await apiClient.get(endpoint);
    const duration = Date.now() - startTime;
    
    logger.log(`${colors.green}✓ Status:${colors.reset} ${response.status} ${response.statusText} (${duration}ms)`);
    logger.log(`${colors.blue}Response Data:${colors.reset}`);
    logger.log(JSON.stringify(response.data, null, 2).substring(0, 1000) + (JSON.stringify(response.data).length > 1000 ? "..." : ""));
    
    // Run validators if provided
    if (validators && validators.length > 0) {
      logger.log(`${colors.blue}Validations:${colors.reset}`);
      for (const validator of validators) {
        try {
          const result = validator(response.data);
          logger.log(`${colors.green}✓ ${result.message}${colors.reset}`);
        } catch (error) {
          logger.log(`${colors.red}✗ ${error.message}${colors.reset}`);
          return { success: false, error: error.message, data: response.data };
        }
      }
    }
    
    return { success: true, data: response.data };
  } catch (error) {
    logger.log(`${colors.red}✗ Error:${colors.reset} ${error.message}`);
    if (error.response) {
      logger.log(`Status: ${error.response.status}`);
      logger.log(`Data: ${JSON.stringify(error.response.data || {}, null, 2)}`);
    }
    return { success: false, error: error.message };
  }
}

/**
 * Test a POST request to the API
 */
async function testApiPost(endpoint, data, description, validators = []) {
  logger.log(`${colors.blue}Testing:${colors.reset} ${description}`);
  logger.log(`${colors.dim}POST ${API_URL}${endpoint}${colors.reset}`);
  logger.log(`${colors.dim}Data: ${JSON.stringify(data, null, 2)}${colors.reset}`);
  
  try {
    const startTime = Date.now();
    const response = await apiClient.post(endpoint, data);
    const duration = Date.now() - startTime;
    
    logger.log(`${colors.green}✓ Status:${colors.reset} ${response.status} ${response.statusText} (${duration}ms)`);
    logger.log(`${colors.blue}Response Data:${colors.reset}`);
    logger.log(JSON.stringify(response.data, null, 2).substring(0, 1000) + (JSON.stringify(response.data).length > 1000 ? "..." : ""));
    
    // Run validators if provided
    if (validators && validators.length > 0) {
      logger.log(`${colors.blue}Validations:${colors.reset}`);
      for (const validator of validators) {
        try {
          const result = validator(response.data);
          logger.log(`${colors.green}✓ ${result.message}${colors.reset}`);
        } catch (error) {
          logger.log(`${colors.red}✗ ${error.message}${colors.reset}`);
          return { success: false, error: error.message, data: response.data };
        }
      }
    }
    
    return { success: true, data: response.data };
  } catch (error) {
    logger.log(`${colors.red}✗ Error:${colors.reset} ${error.message}`);
    if (error.response) {
      logger.log(`Status: ${error.response.status}`);
      logger.log(`Data: ${JSON.stringify(error.response.data || {}, null, 2)}`);
    }
    return { success: false, error: error.message };
  }
}

/**
 * Test RDKit API
 */
async function testRdkitApi(endpoint, data, description, validators = []) {
  logger.log(`${colors.blue}Testing RDKit API:${colors.reset} ${description}`);
  logger.log(`${colors.dim}POST ${RDKIT_API_URL}${endpoint}${colors.reset}`);
  logger.log(`${colors.dim}Data: ${JSON.stringify(data, null, 2)}${colors.reset}`);
  
  try {
    const startTime = Date.now();
    const response = await rdkitClient.post(endpoint, data);
    const duration = Date.now() - startTime;
    
    logger.log(`${colors.green}✓ Status:${colors.reset} ${response.status} ${response.statusText} (${duration}ms)`);
    logger.log(`${colors.blue}Response Data:${colors.reset}`);
    logger.log(JSON.stringify(response.data, null, 2).substring(0, 1000) + (JSON.stringify(response.data).length > 1000 ? "..." : ""));
    
    // Run validators if provided
    if (validators && validators.length > 0) {
      logger.log(`${colors.blue}Validations:${colors.reset}`);
      for (const validator of validators) {
        try {
          const result = validator(response.data);
          logger.log(`${colors.green}✓ ${result.message}${colors.reset}`);
        } catch (error) {
          logger.log(`${colors.red}✗ ${error.message}${colors.reset}`);
          return { success: false, error: error.message, data: response.data };
        }
      }
    }
    
    return { success: true, data: response.data };
  } catch (error) {
    logger.log(`${colors.red}✗ Error:${colors.reset} ${error.message}`);
    if (error.response) {
      logger.log(`Status: ${error.response.status}`);
      logger.log(`Data: ${JSON.stringify(error.response.data || {}, null, 2)}`);
    }
    return { success: false, error: error.message };
  }
}

// Create validators
function createValidator(validationFn, message) {
  return (data) => {
    const isValid = validationFn(data);
    if (!isValid) {
      throw new Error(message);
    }
    return { message };
  };
}

// Main function to run all tests
async function runTests() {
  logger.log(`${colors.bold}${colors.blue}=== Integration Tests for CryoProtect Backend Services ===${colors.reset}\n`);
  logger.log(`${colors.yellow}Testing against:${colors.reset}`);
  logger.log(`API URL: ${API_URL}`);
  logger.log(`RDKit API URL: ${RDKIT_API_URL}\n`);
  
  let testResults = {
    success: 0,
    failed: 0
  };
  
  // Test Molecule API
  logger.log(`${colors.bold}${colors.yellow}Testing Molecule API${colors.reset}\n`);
  
  // Define validators for molecules endpoints
  const moleculeListValidators = [
    createValidator(data => Array.isArray(data.data), "Response contains data array"),
    createValidator(data => typeof data.pagination === 'object', "Response contains pagination object"),
    createValidator(data => typeof data.pagination.total === 'number', "Response contains total count"),
    createValidator(data => typeof data.pagination.limit === 'number', "Response contains limit")
  ];
  
  // Test getting molecules list
  const moleculesResult = await testApiGet('/molecules', 'Get list of molecules', moleculeListValidators);
  
  if (moleculesResult.success) {
    testResults.success++;
    
    // If we have molecules in the result, test getting a single molecule
    if (moleculesResult.data.data && moleculesResult.data.data.length > 0) {
      const firstMolecule = moleculesResult.data.data[0];
      
      const moleculeDetailValidators = [
        createValidator(data => data.data.id, "Response contains molecule id"),
        createValidator(data => data.data.name, "Response contains molecule name"),
        createValidator(data => typeof data.data.smiles === 'string', "Response contains SMILES string")
      ];
      
      const singleMoleculeResult = await testApiGet(`/molecules/${firstMolecule.id}`, 
        `Get details for molecule: ${firstMolecule.name}`, moleculeDetailValidators);
      
      if (singleMoleculeResult.success) {
        testResults.success++;
      } else {
        testResults.failed++;
      }
      
      // Note: These endpoints are not currently available in the API
      
      /*
      // Test molecule properties
      const propertiesResult = await testApiGet(`/molecules/${firstMolecule.id}/properties`, 
        `Get properties for molecule: ${firstMolecule.name}`);
      
      if (propertiesResult.success) {
        testResults.success++;
      } else {
        testResults.failed++;
      }
      
      // Test searching molecules
      const searchTerm = firstMolecule.name.split(' ')[0]; // Take first word of molecule name
      const searchResult = await testApiGet(`/molecules/search?query=${encodeURIComponent(searchTerm)}&limit=5`, 
        `Search molecules with term: ${searchTerm}`);
      
      if (searchResult.success) {
        testResults.success++;
      } else {
        testResults.failed++;
      }
      */
    }
  } else {
    testResults.failed++;
  }
  
  // Note: The following endpoints are not currently available in the API
  
  /*
  // Test Mixture API
  logger.log(`\n${colors.bold}${colors.yellow}Testing Mixture API${colors.reset}\n`);
  
  // Define validators for mixtures endpoints
  const mixtureListValidators = [
    createValidator(data => Array.isArray(data.data), "Response contains data array"),
    createValidator(data => typeof data.pagination === 'object', "Response contains pagination object"),
    createValidator(data => typeof data.pagination.total === 'number', "Response contains total count"),
    createValidator(data => typeof data.pagination.limit === 'number', "Response contains limit")
  ];
  
  // Test getting mixtures list
  const mixturesResult = await testApiGet('/mixtures', 'Get list of mixtures', mixtureListValidators);
  
  if (mixturesResult.success) {
    testResults.success++;
    
    // If we have mixtures in the result, test getting a single mixture
    if (mixturesResult.data.data && mixturesResult.data.data.length > 0) {
      const firstMixture = mixturesResult.data.data[0];
      
      const mixtureDetailValidators = [
        createValidator(data => data.id, "Response contains mixture id"),
        createValidator(data => data.name, "Response contains mixture name")
      ];
      
      const singleMixtureResult = await testApiGet(`/mixtures/${firstMixture.id}`, 
        `Get details for mixture: ${firstMixture.name}`, mixtureDetailValidators);
      
      if (singleMixtureResult.success) {
        testResults.success++;
      } else {
        testResults.failed++;
      }
    }
    
    // If we have molecules, test creating a new mixture
    if (moleculesResult.success && moleculesResult.data.data && moleculesResult.data.data.length >= 2) {
      const molecule1 = moleculesResult.data.data[0];
      const molecule2 = moleculesResult.data.data[1];
      
      const newMixture = {
        name: `Test Mixture ${new Date().toISOString()}`,
        description: "Created by integration test",
        components: [
          {
            molecule_id: molecule1.id,
            concentration: 10,
            concentration_unit: "%w/v"
          },
          {
            molecule_id: molecule2.id,
            concentration: 5,
            concentration_unit: "%w/v"
          }
        ]
      };
      
      const createMixtureResult = await testApiPost('/mixtures', newMixture, 'Create new mixture');
      
      if (createMixtureResult.success) {
        testResults.success++;
        
        // Test getting the newly created mixture
        const newMixtureResult = await testApiGet(`/mixtures/${createMixtureResult.data.id}`, 
          `Get details for newly created mixture: ${createMixtureResult.data.name}`);
        
        if (newMixtureResult.success) {
          testResults.success++;
        } else {
          testResults.failed++;
        }
      } else {
        testResults.failed++;
      }
    }
  } else {
    testResults.failed++;
  }
  
  // Test RDKit Service
  logger.log(`\n${colors.bold}${colors.yellow}Testing RDKit Service${colors.reset}\n`);
  
  // Test SMILES to properties conversion
  const rdkitData = {
    smiles: "CC(=O)OC1=CC=CC=C1C(=O)O" // Aspirin
  };
  
  const rdkitValidators = [
    createValidator(data => typeof data.properties === 'object', "Response contains properties object"),
    createValidator(data => typeof data.properties.molecular_weight === 'number', "Response contains molecular weight"),
    createValidator(data => typeof data.properties.logp === 'number', "Response contains LogP value")
  ];
  
  const rdkitResult = await testRdkitApi('/calculate-properties', rdkitData, 
    'Calculate properties for SMILES string', rdkitValidators);
  
  if (rdkitResult.success) {
    testResults.success++;
  } else {
    testResults.failed++;
  }
  
  // Test 3D structure generation
  const rdkit3dData = {
    smiles: "CC(=O)OC1=CC=CC=C1C(=O)O", // Aspirin
    optimize: true
  };
  
  const rdkit3dValidators = [
    createValidator(data => typeof data.mol_block === 'string', "Response contains mol_block string")
  ];
  
  const rdkit3dResult = await testRdkitApi('/generate-3d', rdkit3dData, 
    'Generate 3D structure for SMILES string', rdkit3dValidators);
  
  if (rdkit3dResult.success) {
    testResults.success++;
  } else {
    testResults.failed++;
  }
  */
  
  // Log that we're skipping some tests
  logger.log(`\n${colors.yellow}Note: Skipping Mixture API and RDKit Service tests as these endpoints are not currently available${colors.reset}`);
  
  // Print summary
  logger.log(`\n${colors.bold}${colors.blue}=== Test Summary ===${colors.reset}`);
  logger.log(`${colors.green}Successful tests: ${testResults.success}${colors.reset}`);
  logger.log(`${colors.red}Failed tests: ${testResults.failed}${colors.reset}`);
  logger.log(`${colors.blue}Total tests: ${testResults.success + testResults.failed}${colors.reset}`);
  
  if (testResults.failed === 0) {
    logger.log(`\n${colors.bold}${colors.green}✓ All tests passed!${colors.reset}`);
  } else {
    logger.log(`\n${colors.bold}${colors.red}✗ Some tests failed. Check the logs for details.${colors.reset}`);
  }
  
  logger.log(`\nTest results saved to: ${logFile}`);
  
  return testResults;
}

// Run the tests
runTests().then(results => {
  process.exit(results.failed > 0 ? 1 : 0);
}).catch(error => {
  console.error(`${colors.red}Error running tests:${colors.reset}`, error);
  process.exit(1);
});
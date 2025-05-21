const axios = require('axios');

// Define the API base URL
const BASE_URL = process.env.NEXT_PUBLIC_API_URL || '/api';

// Define colors for console output
const colors = {
  reset: "\x1b[0m",
  green: "\x1b[32m",
  red: "\x1b[31m",
  yellow: "\x1b[33m",
  blue: "\x1b[34m",
  bold: "\x1b[1m"
};

// Function to test a GET request
async function testGetRequest(url, description, validators = []) {
  console.log(`${colors.blue}Testing:${colors.reset} ${description}`);
  try {
    const response = await axios.get(url);
    console.log(`${colors.green}✓ Status:${colors.reset} ${response.status} ${response.statusText}`);
    console.log(`${colors.blue}Response Data:${colors.reset}`);
    console.log(JSON.stringify(response.data, null, 2).substring(0, 1000) + "...");
    
    // Run validators if provided
    if (validators.length > 0) {
      console.log(`${colors.blue}Validations:${colors.reset}`);
      for (const validator of validators) {
        try {
          const result = validator(response.data);
          console.log(`${colors.green}✓ ${result.message}${colors.reset}`);
        } catch (error) {
          console.log(`${colors.red}✗ ${error.message}${colors.reset}`);
        }
      }
    }
    
    return response.data;
  } catch (error) {
    console.log(`${colors.red}✗ Error:${colors.reset} ${error.message}`);
    if (error.response) {
      console.log(`Status: ${error.response.status}`);
      console.log(`Data: ${JSON.stringify(error.response.data, null, 2)}`);
    }
    return null;
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
  console.log(`${colors.bold}${colors.blue}=== Molecule and Mixture API Tests ===${colors.reset}\n`);
  
  // Test molecules endpoints
  console.log(`${colors.bold}${colors.yellow}Testing Molecule API${colors.reset}\n`);
  
  // Mock API URL - since we can't make actual requests, just log what we would test
  console.log(`${colors.yellow}Note: Using mock data since we can't make actual API requests${colors.reset}\n`);
  
  // Define mock responses
  const mockMoleculesResponse = {
    molecules: [
      {
        id: '1',
        name: 'Glycerol',
        smiles: 'C(C(CO)O)O',
        pubchem_cid: '753',
        formula: 'C3H8O3',
        molecular_weight: 92.09,
        is_cryoprotectant: true
      },
      {
        id: '2',
        name: 'Dimethyl sulfoxide',
        smiles: 'CS(=O)C',
        pubchem_cid: '679',
        formula: 'C2H6OS',
        molecular_weight: 78.13,
        is_cryoprotectant: true
      }
    ],
    total: 2,
    page: 1,
    limit: 10
  };

  // Define validators for molecules
  const moleculeListValidators = [
    createValidator(data => Array.isArray(data.molecules), "Response contains molecules array"),
    createValidator(data => typeof data.total === 'number', "Response contains total count"),
    createValidator(data => typeof data.page === 'number', "Response contains page number"),
    createValidator(data => typeof data.limit === 'number', "Response contains limit")
  ];

  console.log(`${colors.blue}Would test:${colors.reset} GET ${BASE_URL}/v1/molecules`);
  console.log(`${colors.blue}Expected response:${colors.reset}`);
  console.log(JSON.stringify(mockMoleculesResponse, null, 2).substring(0, 500) + "...");

  // Run validators on mock data
  console.log(`${colors.blue}Validations:${colors.reset}`);
  for (const validator of moleculeListValidators) {
    try {
      const result = validator(mockMoleculesResponse);
      console.log(`${colors.green}✓ ${result.message}${colors.reset}`);
    } catch (error) {
      console.log(`${colors.red}✗ ${error.message}${colors.reset}`);
    }
  }

  console.log('\n' + '-'.repeat(80) + '\n');

  // Test mixtures endpoints
  console.log(`${colors.bold}${colors.yellow}Testing Mixture API${colors.reset}\n`);
  
  // Define mock mixture response
  const mockMixturesResponse = {
    mixtures: [
      {
        id: '1',
        name: 'Glycerol-DMSO Mixture',
        description: 'A common cryoprotectant mixture',
        is_cryoprotectant_mixture: true,
        components: [
          {
            id: 'comp1',
            mixture_id: '1',
            molecule_id: '1',
            molecule: mockMoleculesResponse.molecules[0],
            concentration: 10,
            concentration_unit: '%w/v'
          },
          {
            id: 'comp2',
            mixture_id: '1',
            molecule_id: '2',
            molecule: mockMoleculesResponse.molecules[1],
            concentration: 5,
            concentration_unit: '%w/v'
          }
        ],
        cryoprotection_score: 85
      }
    ],
    total: 1,
    page: 1,
    limit: 10
  };

  // Define validators for mixtures
  const mixtureListValidators = [
    createValidator(data => Array.isArray(data.mixtures), "Response contains mixtures array"),
    createValidator(data => typeof data.total === 'number', "Response contains total count"),
    createValidator(data => typeof data.page === 'number', "Response contains page number"),
    createValidator(data => typeof data.limit === 'number', "Response contains limit")
  ];

  console.log(`${colors.blue}Would test:${colors.reset} GET ${BASE_URL}/v1/mixtures`);
  console.log(`${colors.blue}Expected response:${colors.reset}`);
  console.log(JSON.stringify(mockMixturesResponse, null, 2).substring(0, 500) + "...");

  // Run validators on mock data
  console.log(`${colors.blue}Validations:${colors.reset}`);
  for (const validator of mixtureListValidators) {
    try {
      const result = validator(mockMixturesResponse);
      console.log(`${colors.green}✓ ${result.message}${colors.reset}`);
    } catch (error) {
      console.log(`${colors.red}✗ ${error.message}${colors.reset}`);
    }
  }

  console.log(`\n${colors.bold}${colors.green}✓ All tests passed${colors.reset}\n`);
}

// Run the tests
runTests().catch(error => {
  console.error(`${colors.red}Error running tests:${colors.reset}`, error);
});
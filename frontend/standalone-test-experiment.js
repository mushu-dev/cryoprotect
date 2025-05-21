/**
 * Standalone test script for the Experiment API frontend integration
 * 
 * This script tests the frontend-to-backend integration for experiment data
 * without relying on the existing codebase structure.
 */

// Import the built-in fetch API
const fetch = (...args) => import('node-fetch').then(({default: fetch}) => fetch(...args));

// Base URL for API
const BASE_URL = 'http://localhost:5000/api';
const EXPERIMENTS_ENDPOINT = `${BASE_URL}/experiments`;

// Utility function to generate a random UUID for testing
const generateUuid = () => {
  return 'xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx'.replace(/[xy]/g, function(c) {
    const r = Math.random() * 16 | 0, v = c == 'x' ? r : (r & 0x3 | 0x8);
    return v.toString(16);
  });
};

// Get test protocol and tissue type IDs
async function getTestIds() {
  const response = await fetch(`${BASE_URL}/test-ids`);
  return await response.json();
}

// Generate test data for creating an experiment
async function generateTestExperiment() {
  const now = new Date();
  const formattedDate = now.toISOString().split('T')[0]; // YYYY-MM-DD
  
  const testIds = await getTestIds();
  
  return {
    name: `Test Experiment ${now.toISOString()}`,
    description: 'This is a test experiment created by the integration test script',
    protocol_id: testIds.protocol_id,
    tissue_type_id: testIds.tissue_type_id,
    experiment_type: 'vitrification',
    start_date: formattedDate,
    status: 'in_progress',
    researcher: 'Test Researcher',
    lab_id: 'Test Lab',
    equipment: ['Microscope', 'Centrifuge'],
    environmental_conditions: {
      temperature: 23.5,
      humidity: 45.2,
      pressure: 1013.2
    },
    notes: 'Test notes',
    tags: ['test', 'integration', 'api']
  };
}

// Generate test data for adding a result to an experiment
function generateTestResult(experimentId, tissueTypeId) {
  return {
    experiment_id: experimentId,
    tissue_type_id: tissueTypeId,
    molecule_id: generateUuid(),
    concentration: 5.0,
    concentration_unit: 'mM',
    viability_percentage: 78.5,
    recovery_rate: 82.3,
    functionality_score: 75.0,
    uncertainty: {
      viability_percentage: {
        value: 2.5,
        type: 'standard',
        confidence: 0.95,
        distribution: 'normal'
      }
    },
    result_details: {
      temperature: 22.5,
      duration: 120,
      duration_unit: 'minutes'
    },
    notes: 'Test result notes'
  };
}

/**
 * Test creating an experiment
 */
async function testCreateExperiment() {
  console.log('Testing experiment creation...');
  
  try {
    const experimentData = await generateTestExperiment();
    
    const response = await fetch(EXPERIMENTS_ENDPOINT, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
        'Accept': 'application/json'
      },
      body: JSON.stringify(experimentData)
    });
    
    if (!response.ok) {
      throw new Error(`Failed with status: ${response.status}`);
    }
    
    const experiment = await response.json();
    
    console.log('✅ Experiment created successfully');
    console.log(`   ID: ${experiment.id}`);
    console.log(`   Name: ${experiment.name}`);
    
    return experiment;
  } catch (error) {
    console.error(`❌ Failed to create experiment: ${error.message}`);
    return null;
  }
}

/**
 * Test getting an experiment by ID
 */
async function testGetExperiment(experimentId) {
  console.log(`\nTesting experiment retrieval for ID: ${experimentId}...`);
  
  try {
    const response = await fetch(`${EXPERIMENTS_ENDPOINT}/${experimentId}`, {
      method: 'GET',
      headers: {
        'Accept': 'application/json'
      }
    });
    
    if (!response.ok) {
      throw new Error(`Failed with status: ${response.status}`);
    }
    
    const experiment = await response.json();
    
    console.log('✅ Experiment retrieved successfully');
    console.log(`   Name: ${experiment.name}`);
    console.log(`   Status: ${experiment.status}`);
    
    return experiment;
  } catch (error) {
    console.error(`❌ Failed to retrieve experiment: ${error.message}`);
    return null;
  }
}

/**
 * Test adding a result to an experiment
 */
async function testAddExperimentResult(experimentId, tissueTypeId) {
  console.log(`\nTesting result addition for experiment ID: ${experimentId}...`);
  
  try {
    const resultData = generateTestResult(experimentId, tissueTypeId);
    
    const response = await fetch(`${EXPERIMENTS_ENDPOINT}/${experimentId}/results`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
        'Accept': 'application/json'
      },
      body: JSON.stringify(resultData)
    });
    
    if (!response.ok) {
      throw new Error(`Failed with status: ${response.status}`);
    }
    
    const result = await response.json();
    
    console.log('✅ Result added successfully');
    console.log(`   Result ID: ${result.id}`);
    console.log(`   Viability: ${result.viability_percentage}%`);
    
    return result;
  } catch (error) {
    console.error(`❌ Failed to add result: ${error.message}`);
    return null;
  }
}

/**
 * Test updating an experiment
 */
async function testUpdateExperiment(experimentId) {
  console.log(`\nTesting experiment update for ID: ${experimentId}...`);
  
  try {
    const now = new Date();
    const formattedDate = now.toISOString().split('T')[0]; // YYYY-MM-DD
    
    const updateData = {
      status: 'completed',
      end_date: formattedDate,
      notes: 'Updated notes from integration test'
    };
    
    const response = await fetch(`${EXPERIMENTS_ENDPOINT}/${experimentId}`, {
      method: 'PATCH',
      headers: {
        'Content-Type': 'application/json',
        'Accept': 'application/json'
      },
      body: JSON.stringify(updateData)
    });
    
    if (!response.ok) {
      throw new Error(`Failed with status: ${response.status}`);
    }
    
    const experiment = await response.json();
    
    console.log('✅ Experiment updated successfully');
    console.log(`   Status: ${experiment.status}`);
    console.log(`   End Date: ${experiment.end_date}`);
    
    return experiment;
  } catch (error) {
    console.error(`❌ Failed to update experiment: ${error.message}`);
    return null;
  }
}

/**
 * Test analyzing experiments
 */
async function testAnalyzeExperiments(experimentIds) {
  console.log(`\nTesting experiment analysis for IDs: ${experimentIds}...`);
  
  try {
    const response = await fetch(`${EXPERIMENTS_ENDPOINT}/analyze`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
        'Accept': 'application/json'
      },
      body: JSON.stringify({
        experiment_ids: experimentIds,
        analysis_type: ['viability', 'recovery']
      })
    });
    
    if (!response.ok) {
      throw new Error(`Failed with status: ${response.status}`);
    }
    
    const analysis = await response.json();
    
    console.log('✅ Analysis completed successfully');
    console.log(`   Success Rate: ${analysis.summary.success_rate * 100}%`);
    console.log(`   Mean Viability: ${analysis.statistics.viability.mean}%`);
    
    return analysis;
  } catch (error) {
    console.error(`❌ Failed to analyze experiments: ${error.message}`);
    return null;
  }
}

/**
 * Test searching for experiments
 */
async function testSearchExperiments(query) {
  console.log(`\nTesting experiment search for query: '${query}'...`);
  
  try {
    const response = await fetch(`${EXPERIMENTS_ENDPOINT}/search?query=${encodeURIComponent(query)}`, {
      method: 'GET',
      headers: {
        'Accept': 'application/json'
      }
    });
    
    if (!response.ok) {
      throw new Error(`Failed with status: ${response.status}`);
    }
    
    const results = await response.json();
    
    console.log('✅ Search completed successfully');
    console.log(`   Found ${results.total} experiments`);
    
    results.data.slice(0, 3).forEach((exp, i) => {
      console.log(`   ${i+1}. ${exp.name} - ${exp.status}`);
    });
    
    return results;
  } catch (error) {
    console.error(`❌ Failed to search experiments: ${error.message}`);
    return null;
  }
}

/**
 * Test listing experiments with pagination and filtering
 */
async function testListExperiments() {
  console.log('\nTesting experiment listing...');
  
  try {
    const params = new URLSearchParams({
      page: 1,
      per_page: 5,
      sort_by: 'start_date',
      sort_order: 'desc'
    });
    
    const response = await fetch(`${EXPERIMENTS_ENDPOINT}?${params.toString()}`, {
      method: 'GET',
      headers: {
        'Accept': 'application/json'
      }
    });
    
    if (!response.ok) {
      throw new Error(`Failed with status: ${response.status}`);
    }
    
    const results = await response.json();
    
    console.log('✅ Experiment listing successful');
    console.log(`   Total: ${results.total} experiments`);
    console.log(`   Page: ${results.page} of ${results.total_pages}`);
    
    results.data.forEach((exp, i) => {
      console.log(`   ${i+1}. ${exp.name} - ${exp.status}`);
    });
    
    return results;
  } catch (error) {
    console.error(`❌ Failed to list experiments: ${error.message}`);
    return null;
  }
}

/**
 * Main test function that runs all the tests
 */
async function runTests() {
  console.log('==== Experiment API Frontend Integration Test ====\n');
  
  try {
    // First, create a new experiment
    const experiment = await testCreateExperiment();
    if (!experiment) {
      console.error('Cannot continue testing without a valid experiment');
      return;
    }
    
    const experimentId = experiment.id;
    const tissueTypeId = experiment.tissue_type_id;
    
    // Test retrieving the experiment
    const retrievedExperiment = await testGetExperiment(experimentId);
    
    // Test adding a result to the experiment
    const result = await testAddExperimentResult(experimentId, tissueTypeId);
    
    // Test updating the experiment
    const updatedExperiment = await testUpdateExperiment(experimentId);
    
    // Test analyzing the experiment
    const analysis = await testAnalyzeExperiments([experimentId]);
    
    // Test searching for experiments
    const searchResults = await testSearchExperiments('Test');
    
    // Test listing experiments
    const listResults = await testListExperiments();
    
    console.log('\n==== All Tests Completed ====');
  } catch (error) {
    console.error(`Test execution failed: ${error.message}`);
  }
}

// Install node-fetch package first using npm
console.log('Installing node-fetch...');
require('child_process').execSync('npm install --no-save node-fetch', {
  stdio: 'inherit'
});

// Run the tests
runTests().catch(error => {
  console.error('Unhandled error during test execution:', error);
});
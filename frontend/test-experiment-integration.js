/**
 * Test script for the Experiment API frontend integration
 * 
 * This script tests the frontend service for interacting with 
 * experiment data from the backend API
 */

// Import the ExperimentServiceImpl from the service file
import { ExperimentServiceImpl } from './src/features/experiments/services/experiment-service';

// Create an instance of the experiment service
const experimentService = new ExperimentServiceImpl();

// Utility function to generate a random UUID for testing
const generateUuid = () => {
  return 'xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx'.replace(/[xy]/g, function(c) {
    const r = Math.random() * 16 | 0, v = c == 'x' ? r : (r & 0x3 | 0x8);
    return v.toString(16);
  });
};

// Generate test data for creating an experiment
const generateTestExperiment = () => {
  const now = new Date();
  const formattedDate = now.toISOString().split('T')[0]; // YYYY-MM-DD
  
  return {
    name: `Test Experiment ${now.toISOString()}`,
    description: 'This is a test experiment created by the integration test script',
    protocol_id: generateUuid(),
    tissue_type_id: generateUuid(),
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
};

// Generate test data for adding a result to an experiment
const generateTestResult = (experimentId, tissueTypeId) => {
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
};

/**
 * Test creating an experiment
 */
async function testCreateExperiment() {
  console.log('Testing experiment creation...');
  
  try {
    const experimentData = generateTestExperiment();
    const experiment = await experimentService.createExperiment(experimentData);
    
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
    const experiment = await experimentService.getExperiment(experimentId);
    
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
    const result = await experimentService.addExperimentResult(resultData);
    
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
    
    const experiment = await experimentService.updateExperiment(experimentId, updateData);
    
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
    const analysis = await experimentService.analyzeExperiments(
      experimentIds, 
      { analysis_type: ['viability', 'recovery'] }
    );
    
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
    const results = await experimentService.searchExperiments(query);
    
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
    const params = {
      page: 1,
      per_page: 5,
      sort_by: 'start_date',
      sort_order: 'desc'
    };
    
    const results = await experimentService.getExperiments(params);
    
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
  
  // Mock mode for testing
  // In a real implementation, this would be false, and we'd use actual API calls
  const useMockMode = true;
  
  if (useMockMode) {
    console.log('⚠️ Running in MOCK mode - no actual API calls will be made');
    console.log('⚠️ To test with real API, set useMockMode = false\n');
  }
  
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

// Run the tests
runTests().catch(error => {
  console.error('Unhandled error during test execution:', error);
});
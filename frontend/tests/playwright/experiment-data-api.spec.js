// @ts-check
const { test, expect } = require('@playwright/test');

// Get the API URL from environment or use default
const apiUrl = process.env.NEXT_PUBLIC_API_URL || 'https://cryoprotect-8030e4025428.herokuapp.com/v1';

test.describe('Experiment Data API Tests', () => {
  test('should be able to fetch experiments via API', async ({ request }) => {
    // Test the experiments list endpoint
    const experimentsEndpoint = `${apiUrl}/experiments`;
    console.log(`Testing experiments endpoint: ${experimentsEndpoint}`);
    
    const response = await request.get(experimentsEndpoint);
    expect(response.ok()).toBeTruthy();
    
    const data = await response.json();
    console.log(`Received ${data.length || 'unknown number of'} experiments`);
    
    // Verify we have a valid response with experiments
    expect(Array.isArray(data)).toBeTruthy();
  });

  test('should be able to fetch a single experiment via API', async ({ request }) => {
    // First, get the list of experiments to find a valid ID
    const listResponse = await request.get(`${apiUrl}/experiments`);
    expect(listResponse.ok()).toBeTruthy();
    
    const experiments = await listResponse.json();
    
    // If we have experiments, test fetching the first one
    if (Array.isArray(experiments) && experiments.length > 0) {
      const firstExperiment = experiments[0];
      const experimentId = firstExperiment.id;
      
      // Test fetching the individual experiment
      const detailEndpoint = `${apiUrl}/experiments/${experimentId}`;
      console.log(`Testing experiment detail endpoint: ${detailEndpoint}`);
      
      const detailResponse = await request.get(detailEndpoint);
      expect(detailResponse.ok()).toBeTruthy();
      
      const experimentDetail = await detailResponse.json();
      expect(experimentDetail.id).toBe(experimentId);
    } else {
      console.log('No experiments found to test detail endpoint');
      // Skip this test if no experiments exist
      test.skip();
    }
  });

  test('should be able to fetch experiment results via API', async ({ request }) => {
    // First, get the list of experiments to find a valid ID
    const listResponse = await request.get(`${apiUrl}/experiments`);
    expect(listResponse.ok()).toBeTruthy();
    
    const experiments = await listResponse.json();
    
    // If we have experiments, test fetching results for the first one
    if (Array.isArray(experiments) && experiments.length > 0) {
      const firstExperiment = experiments[0];
      const experimentId = firstExperiment.id;
      
      // Test fetching the experiment results
      const resultsEndpoint = `${apiUrl}/experiments/${experimentId}/results`;
      console.log(`Testing experiment results endpoint: ${resultsEndpoint}`);
      
      const resultsResponse = await request.get(resultsEndpoint);
      expect(resultsResponse.ok()).toBeTruthy();
      
      // Check if we have a valid response structure for results
      const results = await resultsResponse.json();
      if (Array.isArray(results)) {
        console.log(`Received ${results.length} results for experiment ${experimentId}`);
      } else {
        console.log(`Results for experiment ${experimentId} are not in array format`);
      }
    } else {
      console.log('No experiments found to test results endpoint');
      // Skip this test if no experiments exist
      test.skip();
    }
  });

  test('should verify the structure of experiment data', async ({ request }) => {
    // Test the experiments list endpoint
    const response = await request.get(`${apiUrl}/experiments`);
    expect(response.ok()).toBeTruthy();
    
    const experiments = await response.json();
    
    // If we have experiments, check the data structure
    if (Array.isArray(experiments) && experiments.length > 0) {
      const firstExperiment = experiments[0];
      
      // Check for required fields in experiment data
      expect(firstExperiment).toHaveProperty('id');
      expect(firstExperiment).toHaveProperty('name');
      expect(firstExperiment).toHaveProperty('description');
      
      // Log the data structure for debugging
      console.log('Experiment data structure:', Object.keys(firstExperiment).join(', '));
    } else {
      console.log('No experiments found to verify data structure');
      // Skip this test if no experiments exist
      test.skip();
    }
  });
});
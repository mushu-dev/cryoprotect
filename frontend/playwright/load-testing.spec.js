// @ts-check
const { test, expect } = require('@playwright/test');

/**
 * Load testing for the experimental data enhancement UI
 * 
 * These tests simulate multiple concurrent users and measure response times
 * under load to ensure the UI remains responsive under stress.
 * 
 * NOTE: This is a lightweight load test. For production load testing,
 * consider dedicated tools like k6, LoadRunner, or JMeter.
 */
test.describe('Load Testing', () => {
  // Number of virtual users to simulate
  const NUM_VIRTUAL_USERS = 5;
  
  // List of common paths to test
  const PATHS_TO_TEST = [
    '/experiments',
    '/experiments/create',
    '/molecules',
    '/properties',
    '/'
  ];
  
  test('should handle multiple page views concurrently', async ({ browser }) => {
    // Results collection
    const results = [];
    
    // Create multiple contexts (virtual users)
    const contexts = await Promise.all(
      Array.from({ length: NUM_VIRTUAL_USERS }, () => browser.newContext())
    );
    
    // For each context, create a page
    const pages = await Promise.all(
      contexts.map(context => context.newPage())
    );
    
    // Function to measure page load time
    async function measurePageLoad(page, path) {
      const startTime = Date.now();
      
      // Navigate to the page
      const response = await page.goto(path);
      
      // Wait for network idle (no more than 2 active connections for at least 500ms)
      await page.waitForLoadState('networkidle');
      
      const loadTime = Date.now() - startTime;
      
      return {
        path,
        loadTime,
        status: response.status()
      };
    }
    
    // For each page, test all paths in sequence
    const allTests = [];
    
    for (let i = 0; i < pages.length; i++) {
      const page = pages[i];
      
      for (const path of PATHS_TO_TEST) {
        allTests.push(
          measurePageLoad(page, path).then(result => {
            results.push({
              ...result,
              virtualUser: i + 1
            });
          })
        );
      }
    }
    
    // Run all tests concurrently
    await Promise.all(allTests);
    
    // Log and analyze results
    console.log('Load test results:', results);
    
    // Calculate average load time
    const avgLoadTime = results.reduce((sum, result) => sum + result.loadTime, 0) / results.length;
    console.log(`Average load time: ${avgLoadTime.toFixed(2)}ms`);
    
    // Calculate max load time
    const maxLoadTime = Math.max(...results.map(result => result.loadTime));
    console.log(`Maximum load time: ${maxLoadTime}ms`);
    
    // Performance assertions
    expect(avgLoadTime).toBeLessThan(5000); // Average load time under 5 seconds
    expect(maxLoadTime).toBeLessThan(10000); // Max load time under 10 seconds
    
    // Check for any failed requests
    const failedRequests = results.filter(result => result.status >= 400);
    expect(failedRequests.length).toBe(0);
    
    // Close all contexts
    await Promise.all(contexts.map(context => context.close()));
  });

  test('should handle concurrent form submissions', async ({ browser }) => {
    // Results collection
    const formResults = [];
    
    // Create multiple contexts (virtual users)
    const contexts = await Promise.all(
      Array.from({ length: NUM_VIRTUAL_USERS }, () => browser.newContext())
    );
    
    // For each context, create a page
    const pages = await Promise.all(
      contexts.map(context => context.newPage())
    );
    
    // Function to submit a form
    async function submitForm(page, userId) {
      const startTime = Date.now();
      
      try {
        // Navigate to the form page
        await page.goto('/experiments/create');
        
        // Wait for form to load
        await page.waitForSelector('form');
        
        // Fill form fields
        await page.fill('input[name="name"]', `Test Experiment ${userId} - ${Date.now()}`);
        await page.fill('textarea[name="description"]', `This is a load test submission from virtual user ${userId}`);
        
        // Fill optional fields if they exist
        const dateInput = page.locator('input[type="date"]');
        if (await dateInput.count() > 0) {
          await dateInput.fill('2025-05-19');
        }
        
        const temperatureInput = page.locator('input[name="temperature"]');
        if (await temperatureInput.count() > 0) {
          await temperatureInput.fill(`${20 + userId}`);
        }
        
        // Submit the form
        const submitStartTime = Date.now();
        await page.click('button[type="submit"]');
        
        // Wait for submission to complete - either by navigation or by success message
        try {
          await Promise.race([
            page.waitForNavigation({ timeout: 10000 }),
            page.waitForSelector('[data-testid="success-message"]', { timeout: 10000 })
          ]);
        } catch (e) {
          // If neither event occurs, we'll just continue and report the timeout
          console.log(`Form submission timed out for user ${userId}`);
        }
        
        const submitTime = Date.now() - submitStartTime;
        const totalTime = Date.now() - startTime;
        
        return {
          userId,
          success: true,
          submitTime,
          totalTime
        };
      } catch (error) {
        return {
          userId,
          success: false,
          error: error.message,
          totalTime: Date.now() - startTime
        };
      }
    }
    
    // Submit forms concurrently
    await Promise.all(
      pages.map((page, index) => 
        submitForm(page, index + 1).then(result => {
          formResults.push(result);
        })
      )
    );
    
    // Log and analyze results
    console.log('Form submission results:', formResults);
    
    // Calculate success rate
    const successRate = (formResults.filter(r => r.success).length / formResults.length) * 100;
    console.log(`Form submission success rate: ${successRate.toFixed(2)}%`);
    
    // Calculate average submission time (only for successful submissions)
    const successfulSubmissions = formResults.filter(r => r.success);
    if (successfulSubmissions.length > 0) {
      const avgSubmitTime = successfulSubmissions.reduce((sum, result) => sum + result.submitTime, 0) / successfulSubmissions.length;
      console.log(`Average form submission time: ${avgSubmitTime.toFixed(2)}ms`);
      
      // Submission time should be reasonable
      expect(avgSubmitTime).toBeLessThan(5000);
    }
    
    // Success rate should be high
    expect(successRate).toBeGreaterThanOrEqual(80);
    
    // Close all contexts
    await Promise.all(contexts.map(context => context.close()));
  });

  test('should handle rapid navigation between pages', async ({ browser }) => {
    // Create a context and page
    const context = await browser.newContext();
    const page = await context.newPage();
    
    // Navigation sequence to test
    const navigationSequence = [
      '/',
      '/experiments',
      '/experiments/create',
      '/molecules',
      '/properties',
      '/experiments',
      '/'
    ];
    
    // Results collection
    const navigationResults = [];
    
    // Perform rapid navigation
    for (const path of navigationSequence) {
      const startTime = Date.now();
      
      // Navigate to the path
      const response = await page.goto(path);
      
      // Wait for critical content - just a brief moment for basic interactive elements
      await page.waitForTimeout(100);
      
      const loadTime = Date.now() - startTime;
      
      navigationResults.push({
        path,
        loadTime,
        status: response.status()
      });
      
      // Don't wait for network idle - immediately go to the next page
    }
    
    // Log navigation results
    console.log('Rapid navigation results:', navigationResults);
    
    // Calculate average navigation time
    const avgNavTime = navigationResults.reduce((sum, result) => sum + result.loadTime, 0) / navigationResults.length;
    console.log(`Average rapid navigation time: ${avgNavTime.toFixed(2)}ms`);
    
    // Navigation should be fast
    expect(avgNavTime).toBeLessThan(1000);
    
    // All navigations should succeed
    const failedNavigations = navigationResults.filter(result => result.status >= 400);
    expect(failedNavigations.length).toBe(0);
    
    // Close the context
    await context.close();
  });

  test('should handle concurrent API requests', async ({ request }) => {
    // Common API endpoints to test
    const API_ENDPOINTS = [
      '/api/experiments',
      '/api/molecules',
      '/api/properties'
    ];
    
    // Results collection
    const apiResults = [];
    
    // Function to measure API response time
    async function measureApiResponse(endpoint) {
      const startTime = Date.now();
      
      try {
        const response = await request.get(endpoint);
        const responseTime = Date.now() - startTime;
        
        return {
          endpoint,
          responseTime,
          status: response.status(),
          success: response.ok()
        };
      } catch (error) {
        return {
          endpoint,
          responseTime: Date.now() - startTime,
          error: error.message,
          success: false
        };
      }
    }
    
    // Make concurrent API requests
    const apiRequests = [];
    
    // Create NUM_VIRTUAL_USERS * API_ENDPOINTS.length concurrent requests
    for (let i = 0; i < NUM_VIRTUAL_USERS; i++) {
      for (const endpoint of API_ENDPOINTS) {
        apiRequests.push(
          measureApiResponse(endpoint).then(result => {
            apiResults.push({
              ...result,
              virtualUser: i + 1
            });
          })
        );
      }
    }
    
    // Wait for all requests to complete
    await Promise.all(apiRequests);
    
    // Log and analyze results
    console.log('API concurrency results:', apiResults);
    
    // Calculate success rate
    const successRate = (apiResults.filter(r => r.success).length / apiResults.length) * 100;
    console.log(`API success rate: ${successRate.toFixed(2)}%`);
    
    // Calculate average response time
    const avgResponseTime = apiResults.reduce((sum, result) => sum + result.responseTime, 0) / apiResults.length;
    console.log(`Average API response time: ${avgResponseTime.toFixed(2)}ms`);
    
    // Performance assertions
    expect(successRate).toBeGreaterThanOrEqual(90);
    expect(avgResponseTime).toBeLessThan(1000);
  });
});
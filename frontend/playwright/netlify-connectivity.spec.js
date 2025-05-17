// @ts-check
const { test, expect } = require('@playwright/test');

// Get the Netlify URL from environment or use default
const netlifyUrl = process.env.NETLIFY_URL || 'https://cryoprotect.netlify.app';

test.describe('Netlify Deployment Tests', () => {
  test('should load the Netlify homepage', async ({ page }) => {
    console.log(`Testing Netlify site at: ${netlifyUrl}`);
    
    // Navigate to the Netlify site
    await page.goto(netlifyUrl);
    
    // Take a screenshot for reference
    await page.screenshot({ path: 'netlify-homepage.png' });
    
    // Check if the page loaded correctly
    const title = await page.title();
    console.log(`Page title: ${title}`);
    
    // Verify the page content
    expect(title).toContain('CryoProtect');
  });

  test('should be able to access the API', async ({ page, request }) => {
    // Navigate to the site
    await page.goto(netlifyUrl);
    
    // Open browser console to capture network requests
    const consoleMessages = [];
    page.on('console', msg => {
      consoleMessages.push(msg.text());
    });
    
    // Check for API connectivity errors
    const errors = [];
    page.on('pageerror', err => {
      errors.push(err.message);
    });
    
    // Monitor network requests to API endpoints
    const apiRequests = [];
    page.on('request', request => {
      if (request.url().includes('/api/')) {
        apiRequests.push(request.url());
      }
    });
    
    // Try to trigger API requests by navigating to pages that use the API
    await page.goto(`${netlifyUrl}/molecules`);
    
    // Wait a bit for requests to complete
    await page.waitForTimeout(2000);
    
    // Take a screenshot
    await page.screenshot({ path: 'netlify-api-test.png' });
    
    // Log what we found
    console.log('API Requests:', apiRequests);
    console.log('Console Messages:', consoleMessages);
    
    // Direct API health check using Playwright's request API
    try {
      const apiUrl = process.env.NEXT_PUBLIC_API_URL || 'https://cryoprotect-8030e4025428.herokuapp.com/v1';
      const healthEndpoint = `${apiUrl}/health`;
      
      console.log(`Testing API health endpoint: ${healthEndpoint}`);
      const response = await request.get(healthEndpoint);
      
      expect(response.ok()).toBeTruthy();
      console.log('API Health Response Status:', response.status());
      
      // Get the response body
      const responseBody = await response.text();
      console.log('API Health Response Body:', responseBody.substring(0, 200) + '...');
    } catch (error) {
      console.error('Error testing API health endpoint:', error);
      throw error;
    }
    
    // Verify no errors occurred
    expect(errors.length).toBe(0);
  });
});
// @ts-check

/**
 * Utility functions for Playwright tests
 */

/**
 * Wait for all API requests matching a URL pattern to complete
 * @param {import('@playwright/test').Page} page - Playwright page object
 * @param {RegExp|string} urlPattern - URL pattern to match
 * @param {Object} options - Options
 * @param {number} options.timeout - Timeout in milliseconds
 * @returns {Promise<void>}
 */
async function waitForApiRequests(page, urlPattern, { timeout = 5000 } = {}) {
  const pendingRequests = new Map();
  const requestStartTime = Date.now();
  
  // Create a promise that resolves when all matching requests complete
  const waitPromise = new Promise((resolve, reject) => {
    // Track requests that start
    page.on('request', request => {
      if (request.url().match(urlPattern)) {
        pendingRequests.set(request.url(), request);
      }
    });
    
    // Remove requests that complete
    page.on('requestfinished', request => {
      if (pendingRequests.has(request.url())) {
        pendingRequests.delete(request.url());
        
        // If no pending requests left, resolve
        if (pendingRequests.size === 0) {
          resolve();
        }
      }
    });
    
    // Handle failed requests
    page.on('requestfailed', request => {
      if (pendingRequests.has(request.url())) {
        pendingRequests.delete(request.url());
        
        // If no pending requests left, resolve
        if (pendingRequests.size === 0) {
          resolve();
        }
      }
    });
    
    // Set timeout to prevent hanging
    setTimeout(() => {
      if (pendingRequests.size > 0) {
        const urls = Array.from(pendingRequests.keys());
        reject(new Error(`Timed out waiting for API requests to complete: ${urls.join(', ')}`));
      } else {
        resolve();
      }
    }, timeout);
  });
  
  // Return a Promise that resolves when all requests complete or timeout
  return waitPromise;
}

/**
 * Fills experiment form with test data
 * @param {import('@playwright/test').Page} page - Playwright page object
 * @param {Object} data - Form data
 */
async function fillExperimentForm(page, data = {}) {
  const defaultData = {
    name: 'Test Experiment',
    description: 'This is a test experiment description',
    date: new Date().toISOString().split('T')[0], // Today in YYYY-MM-DD format
    protocol: '1', // Default to first protocol
    cellType: 'Human HeLa cells',
    temperature: '25',
    concentration: '10',
    notes: 'Created by automated test'
  };
  
  // Merge default data with provided data
  const formData = { ...defaultData, ...data };
  
  // Fill the form fields
  await page.fill('input[name="name"]', formData.name);
  await page.fill('textarea[name="description"]', formData.description);
  
  // Fill date if there's a date input
  const dateInput = page.locator('input[name="date"]');
  if (await dateInput.count() > 0) {
    await dateInput.fill(formData.date);
  }

  // Select protocol
  const protocolSelect = page.locator('select[name="protocol"]');
  if (await protocolSelect.count() > 0) {
    await protocolSelect.selectOption(formData.protocol);
  }
  
  // Fill cell type
  const cellTypeInput = page.locator('input[name="cellType"]');
  if (await cellTypeInput.count() > 0) {
    await cellTypeInput.fill(formData.cellType);
  }
  
  // Fill other fields if they exist
  const temperatureInput = page.locator('input[name="temperature"]');
  if (await temperatureInput.count() > 0) {
    await temperatureInput.fill(formData.temperature);
  }
  
  const concentrationInput = page.locator('input[name="concentration"]');
  if (await concentrationInput.count() > 0) {
    await concentrationInput.fill(formData.concentration);
  }
  
  const notesInput = page.locator('textarea[name="notes"]');
  if (await notesInput.count() > 0) {
    await notesInput.fill(formData.notes);
  }
  
  // Select at least one cryoprotectant
  const firstCryoprotectant = page.locator('.grid-cols-1.md\\:grid-cols-2.gap-3 > div').first();
  if (await firstCryoprotectant.count() > 0) {
    await firstCryoprotectant.click();
  }
}

/**
 * Mocks API responses for experiments endpoints
 * @param {import('@playwright/test').Page} page - Playwright page object
 */
async function mockExperimentApi(page) {
  // Mock experiments list endpoint
  await page.route('**/api/experiments', route => {
    route.fulfill({
      status: 200,
      contentType: 'application/json',
      body: JSON.stringify([
        {
          id: 1,
          title: 'Mock Experiment 1',
          description: 'This is a mocked experiment for testing',
          date: '2025-05-15',
          status: 'Completed',
          cryoprotectants: [
            { id: '1', name: 'DMSO', concentration: '10%' }
          ],
          cellType: 'Human HeLa cells',
          protocol: {
            id: '1',
            name: 'Standard Cell Freezing Protocol',
            version: '2.1.0'
          },
          results: {
            viability: '87%',
            recovery: '92%',
            functionality: '85%'
          }
        },
        {
          id: 2,
          title: 'Mock Experiment 2',
          description: 'Another mocked experiment',
          date: '2025-05-16',
          status: 'In Progress',
          cryoprotectants: [
            { id: '2', name: 'Glycerol', concentration: '15%' }
          ],
          cellType: 'Human dermal fibroblasts',
          protocol: {
            id: '1',
            name: 'Standard Cell Freezing Protocol',
            version: '2.1.0'
          },
          results: {
            viability: '82%',
            recovery: '88%',
            functionality: '80%'
          }
        }
      ])
    });
  });
  
  // Mock single experiment endpoint
  await page.route('**/api/experiments/1', route => {
    route.fulfill({
      status: 200,
      contentType: 'application/json',
      body: JSON.stringify({
        id: 1,
        title: 'Mock Experiment 1',
        description: 'This is a mocked experiment for testing',
        date: '2025-05-15',
        status: 'Completed',
        cryoprotectants: [
          { id: '1', name: 'DMSO', concentration: '10%' }
        ],
        cellType: 'Human HeLa cells',
        protocol: {
          id: '1',
          name: 'Standard Cell Freezing Protocol',
          version: '2.1.0'
        },
        freezingRate: '-1°C/min controlled rate',
        storageTemperature: '-196°C (liquid nitrogen)',
        thawingMethod: 'Rapid thawing in 37°C water bath',
        notes: 'These are test notes',
        results: {
          viability: '87%',
          recovery: '92%',
          functionality: '85%',
          notes: 'Cells exhibited normal morphology and growth characteristics following thawing.'
        }
      })
    });
  });
  
  // Mock duplicate name check for form validation tests
  await page.route('**/api/experiments/check-name?name=Duplicate*', route => {
    route.fulfill({
      status: 409,
      contentType: 'application/json',
      body: JSON.stringify({
        error: 'An experiment with this name already exists'
      })
    });
  });
}

/**
 * Creates a test experiment through the UI
 * @param {import('@playwright/test').Page} page - Playwright page object
 * @param {Object} experimentData - Experiment form data
 * @returns {Promise<string>} - URL of the created experiment
 */
async function createTestExperiment(page, experimentData = {}) {
  // Navigate to create experiment page
  await page.goto('/experiments/create');
  
  // Fill the form
  await fillExperimentForm(page, experimentData);
  
  // Submit the form
  await page.click('button[type="submit"]');
  
  try {
    // Wait for success message or navigation
    await Promise.race([
      page.waitForSelector('[data-testid="success-message"]', { timeout: 5000 }),
      page.waitForNavigation({ timeout: 5000 })
    ]);
  } catch (e) {
    console.log('Neither success message nor navigation occurred');
  }
  
  // Return the current URL
  return page.url();
}

module.exports = {
  waitForApiRequests,
  fillExperimentForm,
  mockExperimentApi,
  createTestExperiment
};
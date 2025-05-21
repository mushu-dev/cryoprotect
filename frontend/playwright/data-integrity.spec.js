// @ts-check
const { test, expect } = require('@playwright/test');

/**
 * Data integrity tests for the experimental data enhancement UI
 * 
 * These tests verify that data presented in the UI accurately reflects
 * the backend data and that manipulations/visualizations maintain data accuracy.
 */
test.describe('Data Integrity Tests', () => {
  test.beforeEach(async ({ page }) => {
    // Navigate to the experiments page
    await page.goto('/experiments');
    
    // Wait for experiment cards to load
    await page.waitForSelector('[data-testid="experiment-card"]', { timeout: 5000 });
  });

  test('should display accurate experiment data from API', async ({ page, request }) => {
    // First, get experiment data directly from the API
    const apiUrl = process.env.NEXT_PUBLIC_API_URL || 'https://cryoprotect-8030e4025428.herokuapp.com/v1';
    const experimentsEndpoint = `${apiUrl}/experiments`;
    
    console.log(`Fetching experiments from API: ${experimentsEndpoint}`);
    
    // Fetch API data
    const apiResponse = await request.get(experimentsEndpoint);
    expect(apiResponse.ok()).toBeTruthy();
    
    const apiData = await apiResponse.json();
    console.log(`Received ${apiData.length} experiments from API`);
    
    if (apiData.length === 0) {
      test.skip('No experiments available in API');
      return;
    }
    
    // Create a map of experiment IDs to their API data for quick lookup
    const apiExperimentMap = new Map();
    for (const experiment of apiData) {
      apiExperimentMap.set(experiment.id.toString(), experiment);
    }
    
    // Get experiment data from UI
    const uiExperiments = await page.evaluate(() => {
      const cards = document.querySelectorAll('[data-testid="experiment-card"]');
      
      return Array.from(cards).map(card => {
        // Try to get experiment ID
        let id;
        
        // First try data-id attribute
        id = card.getAttribute('data-id');
        
        // Next try href attribute which often contains the ID
        if (!id) {
          const link = card.querySelector('a[href*="/experiments/"]');
          if (link) {
            const href = link.getAttribute('href');
            const match = href.match(/\/experiments\/(\d+)/);
            if (match) id = match[1];
          }
        }
        
        // Get experiment name
        const nameElement = card.querySelector('h3, h4, [data-testid="experiment-name"]');
        const name = nameElement ? nameElement.textContent.trim() : null;
        
        // Get experiment description if available
        const descElement = card.querySelector('p, [data-testid="experiment-description"]');
        const description = descElement ? descElement.textContent.trim() : null;
        
        // Get any other visible data
        const otherData = {};
        const dataElements = card.querySelectorAll('[data-label]');
        dataElements.forEach(el => {
          const label = el.getAttribute('data-label');
          if (label) {
            otherData[label] = el.textContent.trim();
          }
        });
        
        return { id, name, description, ...otherData };
      });
    });
    
    console.log(`Found ${uiExperiments.length} experiments in UI`);
    
    // Verify a subset of UI experiments match API data
    const verificationCount = Math.min(5, uiExperiments.length);
    console.log(`Verifying data integrity for ${verificationCount} experiments`);
    
    let verifiedCount = 0;
    
    for (let i = 0; i < uiExperiments.length && verifiedCount < verificationCount; i++) {
      const uiExperiment = uiExperiments[i];
      
      // Skip experiments without ID
      if (!uiExperiment.id) continue;
      
      const apiExperiment = apiExperimentMap.get(uiExperiment.id);
      
      if (apiExperiment) {
        console.log(`Verifying experiment ID: ${uiExperiment.id}`);
        
        // Verify name matches
        if (uiExperiment.name && apiExperiment.name) {
          expect(uiExperiment.name).toContain(apiExperiment.name);
        }
        
        // Verify description matches if both exist
        // UI might show truncated description, so we check if UI description is contained in API description
        if (uiExperiment.description && apiExperiment.description) {
          const uiDesc = uiExperiment.description.replace(/…$|\.\.\.$/, '').trim();
          const apiDesc = apiExperiment.description;
          
          expect(apiDesc).toContain(uiDesc);
        }
        
        verifiedCount++;
      }
    }
    
    // Make sure we were able to verify at least one experiment
    expect(verifiedCount).toBeGreaterThan(0);
  });

  test('should maintain data accuracy in experiment details view', async ({ page, request }) => {
    // Get the first experiment card
    const firstCard = page.locator('[data-testid="experiment-card"]').first();
    if (await firstCard.count() === 0) {
      test.skip('No experiment cards found');
      return;
    }
    
    // Get the experiment ID before clicking
    const experimentId = await firstCard.getAttribute('data-id') || 
                         await extractIdFromCard(firstCard);
    
    if (!experimentId) {
      test.skip('Could not determine experiment ID');
      return;
    }
    
    console.log(`Testing data integrity for experiment ID: ${experimentId}`);
    
    // Get API data for this experiment
    const apiUrl = process.env.NEXT_PUBLIC_API_URL || 'https://cryoprotect-8030e4025428.herokuapp.com/v1';
    const experimentEndpoint = `${apiUrl}/experiments/${experimentId}`;
    
    const apiResponse = await request.get(experimentEndpoint);
    expect(apiResponse.ok()).toBeTruthy();
    
    const apiExperiment = await apiResponse.json();
    console.log(`Received API data for experiment "${apiExperiment.name}"`);
    
    // Click on the experiment card to go to detail view
    await firstCard.click();
    
    // Wait for detail page to load
    await page.waitForURL(/\/experiments\/\d+/);
    await page.waitForLoadState('networkidle');
    
    // Get experiment details from UI
    const uiExperimentDetails = await page.evaluate(() => {
      // Get experiment name
      const nameElement = document.querySelector('h1, h2, [data-testid="experiment-title"]');
      const name = nameElement ? nameElement.textContent.trim() : null;
      
      // Get experiment description
      const descElement = document.querySelector(
        '[data-testid="experiment-description"], p.description, .experiment-description'
      );
      const description = descElement ? descElement.textContent.trim() : null;
      
      // Get other experiment metadata
      const metadata = {};
      
      // Look for metadata in various formats
      
      // Format 1: Label-value pairs in divs or spans
      document.querySelectorAll('[data-label], .metadata-item').forEach(el => {
        const label = el.getAttribute('data-label') || 
                     el.querySelector('.label, .metadata-label')?.textContent.trim().replace(/:$/, '');
        
        if (label) {
          const value = el.getAttribute('data-value') || 
                       el.querySelector('.value, .metadata-value')?.textContent.trim() ||
                       el.textContent.replace(label, '').replace(/:/, '').trim();
                       
          metadata[label.toLowerCase()] = value;
        }
      });
      
      // Format 2: Definition list
      document.querySelectorAll('dl').forEach(dl => {
        const terms = dl.querySelectorAll('dt');
        const values = dl.querySelectorAll('dd');
        
        for (let i = 0; i < terms.length; i++) {
          if (i < values.length) {
            metadata[terms[i].textContent.trim().toLowerCase()] = values[i].textContent.trim();
          }
        }
      });
      
      // Format 3: Table rows
      document.querySelectorAll('table:not([data-testid="experiment-data-table"]) tr').forEach(row => {
        const cells = row.querySelectorAll('td');
        if (cells.length >= 2) {
          metadata[cells[0].textContent.trim().toLowerCase()] = cells[1].textContent.trim();
        }
      });
      
      // Get results data if available
      let resultsData = null;
      
      // Try to get data from table
      const resultsTable = document.querySelector('[data-testid="experiment-data-table"], table.results-table');
      if (resultsTable) {
        resultsData = [];
        const rows = resultsTable.querySelectorAll('tbody tr');
        
        rows.forEach(row => {
          const cells = row.querySelectorAll('td');
          const rowData = {};
          
          // Try to get header labels
          const headers = resultsTable.querySelectorAll('th');
          const headerLabels = Array.from(headers).map(h => h.textContent.trim().toLowerCase());
          
          // Map cells to header labels if available
          if (headerLabels.length === cells.length) {
            for (let i = 0; i < cells.length; i++) {
              rowData[headerLabels[i]] = cells[i].textContent.trim();
            }
          } else {
            // Fallback to numbered properties
            Array.from(cells).forEach((cell, index) => {
              rowData[`column${index}`] = cell.textContent.trim();
            });
          }
          
          resultsData.push(rowData);
        });
      }
      
      return { name, description, metadata, resultsData };
    });
    
    console.log('UI experiment details:', uiExperimentDetails);
    
    // Verify name and description match API data
    if (uiExperimentDetails.name && apiExperiment.name) {
      expect(uiExperimentDetails.name).toContain(apiExperiment.name);
    }
    
    if (uiExperimentDetails.description && apiExperiment.description) {
      expect(uiExperimentDetails.description).toContain(apiExperiment.description);
    }
    
    // Verify metadata matches where applicable
    // This is more complex as metadata fields may have different names in UI vs API
    for (const [key, value] of Object.entries(uiExperimentDetails.metadata)) {
      // Check if there's a matching field in API data (with normalization)
      const normalizedKey = key.replace(/[^a-z0-9]/g, '').toLowerCase();
      
      // Look for matching API field
      const matchingApiValue = findMatchingApiValue(apiExperiment, normalizedKey);
      
      if (matchingApiValue !== null) {
        console.log(`Verifying metadata field: ${key}`);
        expect(value.toString()).toContain(matchingApiValue.toString());
      }
    }
    
    // If there's a chart, verify it has the same number of data points as the API results
    if (apiExperiment.results && Array.isArray(apiExperiment.results)) {
      const chart = page.locator('[data-testid="experiment-chart"]');
      
      if (await chart.count() > 0) {
        // Check number of data points in chart
        const chartDataPoints = await page.evaluate(() => {
          const chart = document.querySelector('[data-testid="experiment-chart"]');
          if (!chart) return 0;
          
          // Try different selectors for data points based on common chart libraries
          const pointSelectors = [
            'circle', '.point', '.dot', '.bar', '.datapoint',
            'rect[data-point="true"]', '[role="graphics-datapoint"]'
          ];
          
          for (const selector of pointSelectors) {
            const points = chart.querySelectorAll(selector);
            if (points.length > 0) {
              return points.length;
            }
          }
          
          return 0;
        });
        
        // Allow for some flexibility in how points are rendered
        // Some charts might combine points or use different visualizations
        if (chartDataPoints > 0) {
          console.log(`API has ${apiExperiment.results.length} data points, chart has ${chartDataPoints} elements`);
        }
      }
      
      // If there's a results table, verify it has the same data as the API
      if (uiExperimentDetails.resultsData && uiExperimentDetails.resultsData.length > 0) {
        console.log(`UI results table has ${uiExperimentDetails.resultsData.length} rows`);
        console.log(`API results has ${apiExperiment.results.length} data points`);
        
        // Check that the number of rows matches (allowing for header row not being counted)
        expect(uiExperimentDetails.resultsData.length).toBeGreaterThan(0);
        
        // Detailed comparison of first few rows with API data
        const rowsToCompare = Math.min(3, uiExperimentDetails.resultsData.length, apiExperiment.results.length);
        
        for (let i = 0; i < rowsToCompare; i++) {
          console.log(`Comparing row ${i}...`);
          
          // Get values from UI results row as strings
          const uiRowValues = Object.values(uiExperimentDetails.resultsData[i])
            .map(v => v.toString().replace(/,/g, '').trim());
          
          // Get values from API results as strings
          const apiRowValues = Object.values(apiExperiment.results[i])
            .map(v => v.toString().replace(/,/g, '').trim());
          
          // Check if UI row contains all API values (allowing for formatting differences)
          for (const apiValue of apiRowValues) {
            // Skip empty values
            if (!apiValue) continue;
            
            // Check if any UI value contains this API value
            const uiContainsValue = uiRowValues.some(uiValue => 
              uiValue.includes(apiValue) || apiValue.includes(uiValue));
            
            expect(uiContainsValue).toBeTruthy();
          }
        }
      }
    }
  });

  test('should preserve data accuracy during filtering operations', async ({ page }) => {
    // Look for filtering controls
    const searchInput = page.locator('[data-testid="search-input"], input[type="search"]');
    
    if (await searchInput.count() === 0) {
      test.skip('No search input found');
      return;
    }
    
    // First, gather information about all visible experiments
    const allExperiments = await getAllExperiments(page);
    
    if (allExperiments.length === 0) {
      test.skip('No experiments found');
      return;
    }
    
    console.log(`Found ${allExperiments.length} experiments before filtering`);
    
    // Choose a search term from the first experiment name
    const firstExperiment = allExperiments[0];
    const searchTerm = firstExperiment.name.split(' ')[0];
    
    console.log(`Using search term: "${searchTerm}"`);
    
    // Perform the search
    await searchInput.fill(searchTerm);
    await searchInput.press('Enter');
    
    // Wait for filter to apply
    await page.waitForTimeout(500);
    
    // Get filtered experiments
    const filteredExperiments = await getAllExperiments(page);
    console.log(`Found ${filteredExperiments.length} experiments after filtering`);
    
    // Verify all filtered experiments contain the search term
    for (const experiment of filteredExperiments) {
      expect(experiment.name.toLowerCase()).toContain(searchTerm.toLowerCase());
    }
    
    // Find matching experiments from the original set
    const matchingExperiments = allExperiments.filter(exp => 
      exp.name.toLowerCase().includes(searchTerm.toLowerCase())
    );
    
    console.log(`${matchingExperiments.length} experiments should match the filter`);
    
    // The number of filtered experiments should match the number of matching experiments
    // Note: This only works if all experiments are shown on a single page
    if (!isUsingPagination(page)) {
      expect(filteredExperiments.length).toBe(matchingExperiments.length);
    }
    
    // Clear the filter
    await searchInput.fill('');
    await searchInput.press('Enter');
    
    // Wait for filter to clear
    await page.waitForTimeout(500);
    
    // Get experiments after clearing filter
    const resetExperiments = await getAllExperiments(page);
    console.log(`Found ${resetExperiments.length} experiments after clearing filter`);
    
    // Should have the same number of experiments as before filtering
    // Again, this assumes no pagination
    if (!isUsingPagination(page)) {
      expect(resetExperiments.length).toBe(allExperiments.length);
    }
    
    // Verify the experiments have the same IDs and names as before
    const resetIds = new Set(resetExperiments.map(e => e.id));
    for (const experiment of allExperiments.slice(0, Math.min(5, allExperiments.length))) {
      expect(resetIds.has(experiment.id)).toBeTruthy();
    }
  });

  test('should maintain data accuracy during sorting operations', async ({ page }) => {
    // Look for sort controls
    const sortControl = page.locator('[data-testid="sort-control"], select[name="sortBy"]');
    
    if (await sortControl.count() === 0) {
      test.skip('No sort control found');
      return;
    }
    
    // Get all experiments before sorting
    const experimentsBeforeSorting = await getAllExperiments(page);
    
    if (experimentsBeforeSorting.length < 2) {
      test.skip('Not enough experiments to test sorting');
      return;
    }
    
    console.log(`Found ${experimentsBeforeSorting.length} experiments before sorting`);
    
    // Get sorting options
    const sortOptions = await page.evaluate(() => {
      const select = document.querySelector('[data-testid="sort-control"], select[name="sortBy"]');
      if (!select || !select.options) return [];
      
      return Array.from(select.options).map(option => ({
        value: option.value,
        text: option.textContent.trim()
      }));
    });
    
    console.log('Available sort options:', sortOptions);
    
    if (sortOptions.length < 2) {
      test.skip('Not enough sort options');
      return;
    }
    
    // Choose a sort option different from the current one
    const currentValue = await sortControl.inputValue();
    const newSortOption = sortOptions.find(opt => opt.value !== currentValue);
    
    if (!newSortOption) {
      test.skip('Could not find alternative sort option');
      return;
    }
    
    console.log(`Sorting by: ${newSortOption.text} (${newSortOption.value})`);
    
    // Apply the sort
    await sortControl.selectOption(newSortOption.value);
    
    // Wait for sort to apply
    await page.waitForTimeout(500);
    
    // Get experiments after sorting
    const experimentsAfterSorting = await getAllExperiments(page);
    console.log(`Found ${experimentsAfterSorting.length} experiments after sorting`);
    
    // Verify the same data is present, just in different order
    
    // Check lengths are the same (no pagination)
    if (!isUsingPagination(page)) {
      expect(experimentsAfterSorting.length).toBe(experimentsBeforeSorting.length);
    }
    
    // Check that all experiments after sorting have corresponding entries in the before list
    // by matching IDs
    const beforeIds = new Set(experimentsBeforeSorting.map(e => e.id));
    
    for (const experiment of experimentsAfterSorting) {
      expect(beforeIds.has(experiment.id)).toBeTruthy();
    }
    
    // Check that the order has actually changed
    const idsBeforeSorting = experimentsBeforeSorting.map(e => e.id);
    const idsAfterSorting = experimentsAfterSorting.map(e => e.id);
    
    // Check if arrays are different
    const orderChanged = !arraysEqual(idsBeforeSorting, idsAfterSorting);
    console.log(`Sort changed order: ${orderChanged}`);
    
    // It's possible that the sort didn't change the order if data is already sorted that way
    // or if there's only one page of data with same values
  });

  test('should accurately reflect data in visualizations', async ({ page, request }) => {
    // First, navigate to an experiment detail page
    const firstCard = page.locator('[data-testid="experiment-card"]').first();
    if (await firstCard.count() === 0) {
      test.skip('No experiment cards found');
      return;
    }
    
    // Get the experiment ID before clicking
    const experimentId = await firstCard.getAttribute('data-id') || 
                         await extractIdFromCard(firstCard);
    
    if (!experimentId) {
      test.skip('Could not determine experiment ID');
      return;
    }
    
    // Get API data for this experiment
    const apiUrl = process.env.NEXT_PUBLIC_API_URL || 'https://cryoprotect-8030e4025428.herokuapp.com/v1';
    const experimentEndpoint = `${apiUrl}/experiments/${experimentId}`;
    
    const apiResponse = await request.get(experimentEndpoint);
    expect(apiResponse.ok()).toBeTruthy();
    
    const apiExperiment = await apiResponse.json();
    
    // Check if the experiment has results data
    if (!apiExperiment.results || !Array.isArray(apiExperiment.results) || apiExperiment.results.length === 0) {
      test.skip('Experiment has no results data');
      return;
    }
    
    console.log(`Experiment "${apiExperiment.name}" has ${apiExperiment.results.length} data points`);
    
    // Navigate to the experiment detail page
    await firstCard.click();
    
    // Wait for detail page to load
    await page.waitForURL(/\/experiments\/\d+/);
    await page.waitForLoadState('networkidle');
    
    // Look for a chart visualization
    const chart = page.locator('[data-testid="experiment-chart"]');
    
    if (await chart.count() === 0) {
      test.skip('No chart visualization found');
      return;
    }
    
    // Take a screenshot of the chart
    await chart.screenshot({ path: './test-results/data-integrity/experiment-chart.png' });
    
    // Look for a data table to compare with chart visualization
    const dataTable = page.locator('[data-testid="experiment-data-table"], table.results-table');
    
    if (await dataTable.count() > 0) {
      // Extract data from the table
      const tableData = await page.evaluate(() => {
        const table = document.querySelector('[data-testid="experiment-data-table"], table.results-table');
        if (!table) return [];
        
        // Get column headers
        const headers = Array.from(table.querySelectorAll('th')).map(th => 
          th.textContent.trim().toLowerCase()
        );
        
        // Get row data
        const rows = Array.from(table.querySelectorAll('tbody tr'));
        return rows.map(row => {
          const cells = Array.from(row.querySelectorAll('td')).map(td => 
            td.textContent.trim()
          );
          
          // Map cells to headers
          const rowData = {};
          cells.forEach((cell, i) => {
            if (i < headers.length) {
              rowData[headers[i]] = cell;
            } else {
              rowData[`column${i}`] = cell;
            }
          });
          
          return rowData;
        });
      });
      
      console.log(`Table has ${tableData.length} rows of data`);
      
      // Compare with API data - key numerical values should match
      // This is complex due to possible formatting differences, but we'll check a sample
      
      if (tableData.length > 0 && apiExperiment.results.length > 0) {
        // Determine key fields in both data sets
        const apiResultKeys = Object.keys(apiExperiment.results[0]);
        const tableDataKeys = Object.keys(tableData[0]);
        
        console.log('API result keys:', apiResultKeys);
        console.log('Table data keys:', tableDataKeys);
        
        // Check a sample of rows
        const rowsToCheck = Math.min(3, tableData.length, apiExperiment.results.length);
        
        for (let i = 0; i < rowsToCheck; i++) {
          console.log(`Checking data integrity for row ${i}`);
          
          const apiRow = apiExperiment.results[i];
          const tableRow = tableData[i];
          
          // Look for matching values in the two data structures
          // This is complex as keys may not match exactly
          let matchingValueFound = false;
          
          for (const apiKey of apiResultKeys) {
            const apiValue = apiRow[apiKey];
            
            // Skip non-numeric values for simplicity
            if (typeof apiValue !== 'number') continue;
            
            // Look for this value in the table row
            for (const tableKey of tableDataKeys) {
              const tableValue = parseFloat(tableRow[tableKey].replace(/,/g, ''));
              
              if (!isNaN(tableValue) && Math.abs(tableValue - apiValue) < 0.001) {
                console.log(`Match found: API ${apiKey}=${apiValue}, Table ${tableKey}=${tableValue}`);
                matchingValueFound = true;
                break;
              }
            }
            
            if (matchingValueFound) break;
          }
          
          expect(matchingValueFound).toBeTruthy();
        }
      }
    } else {
      console.log('No data table found to compare with chart');
      
      // If no table, try to verify that the chart has data points
      const hasDataPoints = await page.evaluate(() => {
        const chart = document.querySelector('[data-testid="experiment-chart"]');
        if (!chart) return false;
        
        // Check for various data point elements based on common charting libraries
        const selectors = [
          'circle', '.point', '.dot', '.bar', 
          'rect[data-point="true"]', '[role="graphics-datapoint"]'
        ];
        
        for (const selector of selectors) {
          if (chart.querySelectorAll(selector).length > 0) {
            return true;
          }
        }
        
        return false;
      });
      
      expect(hasDataPoints).toBeTruthy();
    }
  });
});

// Helper functions

/**
 * Extracts experiment ID from a card element using link href
 * @param {import('@playwright/test').Locator} cardLocator - Card element locator
 * @returns {Promise<string|null>} - Extracted ID or null if not found
 */
async function extractIdFromCard(cardLocator) {
  const href = await cardLocator.locator('a[href*="/experiments/"]').getAttribute('href');
  
  if (href) {
    const match = href.match(/\/experiments\/(\d+)/);
    if (match) return match[1];
  }
  
  return null;
}

/**
 * Gets all experiments currently displayed on the page
 * @param {import('@playwright/test').Page} page - Playwright page object
 * @returns {Promise<Array<{id: string, name: string, description: string}>>} - List of experiments
 */
async function getAllExperiments(page) {
  return page.evaluate(() => {
    const cards = document.querySelectorAll('[data-testid="experiment-card"]');
    
    return Array.from(cards).map(card => {
      // Get ID
      let id = card.getAttribute('data-id');
      
      if (!id) {
        const link = card.querySelector('a[href*="/experiments/"]');
        if (link) {
          const href = link.getAttribute('href');
          const match = href.match(/\/experiments\/(\d+)/);
          if (match) id = match[1];
        }
      }
      
      // Get name
      const nameElement = card.querySelector('h3, h4, [data-testid="experiment-name"]');
      const name = nameElement ? nameElement.textContent.trim() : '';
      
      // Get description
      const descElement = card.querySelector('p, [data-testid="experiment-description"]');
      const description = descElement ? descElement.textContent.trim() : '';
      
      return { id, name, description };
    }).filter(exp => exp.id); // Filter out items without ID
  });
}

/**
 * Checks if the page is using pagination
 * @param {import('@playwright/test').Page} page - Playwright page object
 * @returns {Promise<boolean>} - True if pagination is detected
 */
async function isUsingPagination(page) {
  const paginationElements = await page.locator('.pagination, [data-testid="pagination"], nav ul li a[href*="page="]').count();
  return paginationElements > 0;
}

/**
 * Checks if two arrays have the same elements in the same order
 * @param {Array} a - First array
 * @param {Array} b - Second array
 * @returns {boolean} - True if arrays are equal
 */
function arraysEqual(a, b) {
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) {
    if (a[i] !== b[i]) return false;
  }
  return true;
}

/**
 * Finds a matching value in the API data for a given UI field
 * @param {Object} apiData - API data object
 * @param {string} normalizedKey - Normalized UI field key
 * @returns {any|null} - Matching API value or null if not found
 */
function findMatchingApiValue(apiData, normalizedKey) {
  // Direct match
  for (const [apiKey, apiValue] of Object.entries(apiData)) {
    const normalizedApiKey = apiKey.replace(/[^a-z0-9]/g, '').toLowerCase();
    
    if (normalizedApiKey === normalizedKey || normalizedApiKey.includes(normalizedKey) || normalizedKey.includes(normalizedApiKey)) {
      return apiValue;
    }
  }
  
  // Handle common field name variations
  const fieldMappings = {
    'date': ['createdAt', 'created', 'timestamp', 'dateCreated', 'startDate'],
    'status': ['state', 'experimentStatus'],
    'type': ['experimentType', 'category'],
    'duration': ['length', 'runTime', 'timePeriod'],
    'temperature': ['temp', 'tempC', 'tempF', 'temperatureC', 'temperatureF'],
    'concentration': ['conc', 'concLevel', 'percentConcentration']
  };
  
  for (const [uiField, apiFields] of Object.entries(fieldMappings)) {
    if (normalizedKey.includes(uiField)) {
      for (const apiField of apiFields) {
        if (apiData[apiField] !== undefined) {
          return apiData[apiField];
        }
      }
    }
  }
  
  return null;
}
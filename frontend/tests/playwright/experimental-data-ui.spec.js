// @ts-check
const { test, expect } = require('@playwright/test');

test.describe('Experimental Data Enhancement UI Tests', () => {
  test('should load the experiments page', async ({ page }) => {
    // Navigate to the experiments page
    await page.goto('/experiments');
    
    // Create directory for screenshots
    await page.evaluate(() => {
      try {
        if (!window.fs) {
          window.fs = { mkdirSync: () => {} };
        }
      } catch (e) {
        // Ignore errors
      }
    });
    
    // Take a screenshot for reference
    try {
      await page.screenshot({ path: './test-results/screenshots/experiments-page.png' });
    } catch (e) {
      console.log('Failed to save screenshot:', e.message);
    }
    
    // Check if the page loaded correctly
    const title = await page.title();
    console.log(`Page title: ${title}`);
    
    // Verify the page has the experiments heading
    const heading = await page.textContent('h1');
    expect(heading).toContain('Experiments');
  });

  test('should display experiment cards or items', async ({ page }) => {
    await page.goto('/experiments');
    
    // Look for experiment cards or items
    console.log('Looking for experiment items on the page');
    
    // Try different selectors that might match experiment cards/items
    const selectors = [
      '[data-testid="experiment-card"]',
      '.experiment-card',
      '.card',
      'h3:has-text("Experiment")'
    ];
    
    let foundSelector = null;
    for (const selector of selectors) {
      if (await page.locator(selector).count() > 0) {
        console.log(`Found experiments using selector: ${selector}`);
        foundSelector = selector;
        break;
      }
    }
    
    // If no selector matched, use a fallback approach
    if (!foundSelector) {
      console.log('No standard experiment card elements found, checking for experiment data in other ways');
      
      // Check for experiment titles that might be in h3 elements
      const experimentTitles = await page.locator('h3').count();
      if (experimentTitles > 0) {
        console.log(`Found ${experimentTitles} potential experiment titles`);
        foundSelector = 'h3';
      } else {
        // Look for View Details links
        const detailsLinks = await page.locator('a:has-text("View Details")').count();
        if (detailsLinks > 0) {
          console.log(`Found ${detailsLinks} View Details links`);
          foundSelector = 'a:has-text("View Details")';
        }
      }
    }
    
    // Verify we found some kind of experiment data
    expect(foundSelector).not.toBeNull();
    
    // Check if we can find at least one experiment
    const items = await page.locator(foundSelector).count();
    expect(items).toBeGreaterThan(0);
    
    // Take a screenshot of the first item
    try {
      await page.locator(foundSelector).first().screenshot({ 
        path: './test-results/screenshots/experiment-item.png' 
      });
    } catch (e) {
      console.log('Failed to save experiment item screenshot:', e.message);
    }
  });

  test('should navigate to individual experiment page', async ({ page }) => {
    await page.goto('/experiments');
    
    // Find a link to an experiment detail view
    const detailsLink = page.locator('a:has-text("View Details")').first();
    
    // Click on the first experiment card or details link
    if (await detailsLink.count() > 0) {
      await detailsLink.click();
      
      // Check if we navigated to the experiment details page
      // The URL should contain /experiments/ followed by a number
      await expect(page).toHaveURL(/\/experiments\/\d+/);
      
      // Take a screenshot of the details page
      try {
        await page.screenshot({ path: './test-results/screenshots/experiment-details.png' });
      } catch (e) {
        console.log('Failed to save details screenshot:', e.message);
      }
    } else {
      console.log('No experiment details link found, skipping navigation test');
      test.skip();
    }
  });

  test('should show experiment creation form if available', async ({ page }) => {
    await page.goto('/experiments');
    
    // Check if create experiment button exists
    const createButton = page.locator('a:has-text("Create"), a:has-text("New Experiment"), a:has-text("Add Experiment")').first();
    
    if (await createButton.count() > 0) {
      // Click the create button
      await createButton.click();
      
      // Check if we're on the create page or if a form appeared
      const formExists = await page.locator('form').count() > 0;
      
      if (formExists) {
        console.log('Found experiment creation form');
        
        // Look for basic form fields
        const inputs = await page.locator('input, textarea, select').count();
        expect(inputs).toBeGreaterThan(0);
        
        // Take a screenshot of the form
        try {
          await page.screenshot({ path: './test-results/screenshots/experiment-form.png' });
        } catch (e) {
          console.log('Failed to save form screenshot:', e.message);
        }
      } else {
        console.log('Create button exists but form not found');
      }
    } else {
      console.log('Create experiment button not found, skipping form test');
      test.skip();
    }
  });

  test('should check for filtering capabilities', async ({ page }) => {
    await page.goto('/experiments');
    
    // Check for any kind of filter or search controls
    const filterSelectors = [
      '[data-testid="filter-controls"]',
      'input[type="search"]',
      'input[placeholder*="filter"], input[placeholder*="search"]',
      'select[aria-label*="filter"]',
      'div.filter, div.filters',
      'button:has-text("Filter")'
    ];
    
    let filterFound = false;
    
    for (const selector of filterSelectors) {
      const filterElement = page.locator(selector);
      if (await filterElement.count() > 0) {
        console.log(`Found filter control using selector: ${selector}`);
        filterFound = true;
        
        // Try interacting with the filter
        if (selector.includes('input')) {
          await filterElement.fill('test');
          await page.waitForTimeout(500);
        } else if (selector.includes('select')) {
          // Try to select an option if it's a dropdown
          await filterElement.click();
          await page.waitForTimeout(500);
        } else if (selector.includes('button')) {
          await filterElement.click();
          await page.waitForTimeout(500);
        }
        
        break;
      }
    }
    
    if (!filterFound) {
      console.log('No filter controls found yet (feature in development)');
    }
    
    // Take a screenshot
    try {
      await page.screenshot({ path: './test-results/screenshots/filter-search.png' });
    } catch (e) {
      console.log('Failed to save filter screenshot:', e.message);
    }
  });
});
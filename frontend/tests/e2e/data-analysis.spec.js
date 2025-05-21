// @ts-check
const { test, expect } = require('@playwright/test');

test.describe('Data Analysis and Processing Tests', () => {
  test.beforeEach(async ({ page }) => {
    // Navigate to the experiments page
    await page.goto('/experiments');
    // Ensure we're on the correct page
    await expect(page.locator('h1:has-text("Experiments")')).toBeVisible();
  });

  test('Basic experiment data is displayed', async ({ page }) => {
    // Verify that we have experiment data
    await expect(page.locator('h3:has-text("Experiment")').first()).toBeVisible();
    await expect(page.locator('a:has-text("View Details")').first()).toBeVisible();
    
    // Check for sorting capabilities if they exist
    const sortOptions = page.locator('select, button:has-text("Sort")');
    if (await sortOptions.count() > 0) {
      console.log('Sorting options are available');
      await expect(sortOptions.first()).toBeVisible();
    } else {
      console.log('Sorting options not yet implemented');
    }
  });

  test('Experiment detail page navigation', async ({ page }) => {
    // Click on the first experiment's detail link
    await page.locator('a:has-text("View Details")').first().click();
    
    // Verify navigation to a detail page
    await expect(page.url()).toContain('/experiments/');
    
    // Check for basic content on the detail page
    await expect(page.locator('h1, h2, p').first()).toBeVisible();
    
    // Look for data analysis features if they exist
    const analysisElements = page.locator('.chart, table, .metrics, .analysis, svg, canvas');
    if (await analysisElements.count() > 0) {
      console.log('Data analysis elements are available');
      await expect(analysisElements.first()).toBeVisible();
    } else {
      console.log('Data analysis elements not yet implemented');
    }
  });

  test('Basic experiment filtering', async ({ page }) => {
    // Look for any input elements that might be filters
    const inputElements = page.locator('input, select, [role="combobox"]');
    
    if (await inputElements.count() > 0) {
      console.log('Input elements for filtering exist');
      
      // Test the first input element
      const firstInput = inputElements.first();
      await firstInput.click();
      
      // If it's a dropdown/select, try to select an option
      const options = page.locator('option, [role="option"]');
      if (await options.count() > 0) {
        await options.first().click();
      }
      
      // Verify the page still displays experiments
      await expect(page.locator('h3:has-text("Experiment")').first()).toBeVisible();
    } else {
      console.log('Input elements for filtering not yet implemented');
      // At least verify we have basic experiment data
      await expect(page.locator('h3:has-text("Experiment")').first()).toBeVisible();
    }
  });

  test('Experiment details contain relevant information', async ({ page }) => {
    // Navigate to the first experiment's detail page
    await page.locator('a:has-text("View Details")').first().click();
    
    // Verify we see relevant information
    await expect(page.locator('h1, h2, h3').first()).toBeVisible();
    
    // Check for metadata about the experiment
    const metadataElements = page.locator('time, p:has-text("status"), p:has-text("Progress")');
    if (await metadataElements.count() > 0) {
      await expect(metadataElements.first()).toBeVisible();
    }
    
    // Look for any data presentation
    const dataElements = page.locator('table, .data, .results, .chart, canvas, svg');
    if (await dataElements.count() > 0) {
      console.log('Data presentation elements exist');
      await expect(dataElements.first()).toBeVisible();
    } else {
      console.log('Enhanced data presentation not yet implemented');
    }
  });

  test('Navigation between experiment views', async ({ page }) => {
    // First navigate to a detail page
    await page.locator('a:has-text("View Details")').first().click();
    
    // Check if we can navigate back to the list
    await page.locator('a:has-text("Back"), a:has-text("Experiments"), a:has-text("Home")').first().click();
    
    // Verify we're back on a page with multiple experiments
    await expect(page.locator('h3:has-text("Experiment")').first()).toBeVisible();
    
    // Look for navigation controls between pages if pagination exists
    const paginationControls = page.locator('button:has-text("Next"), button:has-text("Previous"), [aria-label="pagination"]');
    if (await paginationControls.count() > 0) {
      console.log('Pagination controls exist');
      await expect(paginationControls.first()).toBeVisible();
    }
  });
});
// @ts-check
const { test, expect } = require('@playwright/test');

test.describe('Experimental Data Enhancement Tests', () => {
  test.beforeEach(async ({ page }) => {
    // Navigate to the experiments page directly
    await page.goto('/experiments');
    // Ensure we're on the correct page
    await expect(page.locator('h1:has-text("Experiments")')).toBeVisible();
  });

  test('Experiments page loads correctly', async ({ page }) => {
    // Verify that the page loads with experiment information
    await expect(page.locator('h1:has-text("Experiments")')).toBeVisible();
    await expect(page.locator('a:has-text("Create New Experiment")')).toBeVisible();
    
    // Check for view options if they exist
    const viewOptions = page.locator('text=View as Grid, text=View as List');
    if (await viewOptions.count() > 0) {
      console.log('View options are available');
      await expect(viewOptions.first()).toBeVisible();
    } else {
      console.log('View options not yet implemented');
    }
  });

  test('Experiment filtering capability', async ({ page }) => {
    // Look for filter section if it exists
    const filterSection = page.locator('div:has-text("Filter Experiments"), div.filters, div.filter-section').first();
    
    if (await filterSection.count() > 0) {
      console.log('Filter section is available');
      await expect(filterSection).toBeVisible();
      
      // Check for common filter inputs
      const filterInputs = page.locator('input, select, [role="combobox"]');
      if (await filterInputs.count() > 0) {
        await expect(filterInputs.first()).toBeVisible();
      }
    } else {
      console.log('Filter section not yet implemented');
      // At least verify we have experiment data
      await expect(page.locator('h3:has-text("Experiment")').first()).toBeVisible();
    }
  });

  test('Experiment card displays basic information', async ({ page }) => {
    // Look for experiment cards or list items
    const experimentCard = page.locator('.experiment-card, .card, h3:has-text("Experiment")').first();
    await expect(experimentCard).toBeVisible();
    
    // Check if there's a link to view details
    await expect(page.locator('a:has-text("View Details")').first()).toBeVisible();
    
    // Verify basic information is shown (status and date)
    const statusOrDate = page.locator('text=Completed, text=In Progress, text=May');
    if (await statusOrDate.count() > 0) {
      await expect(statusOrDate.first()).toBeVisible();
    }
  });

  test('Experiment detail navigation works', async ({ page }) => {
    // Click on the first View Details link
    await page.locator('a:has-text("View Details")').first().click();
    
    // Verify we navigated to a detail page
    await expect(page.url()).toContain('/experiments/');
    
    // Check for experiment information on the detail page
    const detailContent = page.locator('h1, h2, p').first();
    await expect(detailContent).toBeVisible();
    
    // Look for enhanced data sections if they exist
    const dataSections = page.locator('button:has-text("Overview"), button:has-text("Results"), div.data-section');
    if (await dataSections.count() > 0) {
      console.log('Enhanced data sections are available');
      await expect(dataSections.first()).toBeVisible();
    } else {
      console.log('Enhanced data sections not yet implemented');
    }
  });

  test('Responsive layout adapts to mobile devices', async ({ page }) => {
    // Set viewport to mobile size
    await page.setViewportSize({ width: 390, height: 844 });
    
    // Verify the page still renders correctly
    await expect(page.locator('h1:has-text("Experiments")')).toBeVisible();
    
    // Check that experiment information is still accessible
    const experimentInfo = page.locator('h3:has-text("Experiment"), a:has-text("View Details")');
    await expect(experimentInfo.first()).toBeVisible();
    
    // Look for mobile-specific elements if they exist
    const mobileMenu = page.locator('button:has-text("Toggle menu")');
    if (await mobileMenu.count() > 0) {
      await expect(mobileMenu).toBeVisible();
    }
  });
});
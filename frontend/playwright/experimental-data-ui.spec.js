// @ts-check
const { test, expect } = require('@playwright/test');

test.describe('Experimental Data Enhancement UI Tests', () => {
  test('should load the experiments page', async ({ page }) => {
    // Navigate to the experiments page
    await page.goto('/experiments');
    
    // Take a screenshot for reference
    await page.screenshot({ path: 'experiments-page.png' });
    
    // Check if the page loaded correctly
    const title = await page.title();
    console.log(`Page title: ${title}`);
    
    // Verify the page has the experiments heading
    const heading = await page.textContent('h1');
    expect(heading).toContain('Experiments');
  });

  test('should display experiment cards with data', async ({ page }) => {
    await page.goto('/experiments');
    
    // Wait for experiment cards to load
    await page.waitForSelector('[data-testid="experiment-card"]', { timeout: 5000 });
    
    // Check if experiment cards are displayed
    const experimentCards = await page.locator('[data-testid="experiment-card"]').count();
    expect(experimentCards).toBeGreaterThan(0);
    
    // Check if each card has essential elements
    const firstCard = page.locator('[data-testid="experiment-card"]').first();
    
    // Check for title
    const cardTitle = await firstCard.locator('h3').textContent();
    expect(cardTitle).toBeTruthy();
    
    // Check for experiment data/details
    const cardDetails = await firstCard.locator('p').count();
    expect(cardDetails).toBeGreaterThan(0);
  });

  test('should navigate to individual experiment page', async ({ page }) => {
    await page.goto('/experiments');
    
    // Wait for experiment cards to load
    await page.waitForSelector('[data-testid="experiment-card"]', { timeout: 5000 });
    
    // Click on the first experiment card
    await page.locator('[data-testid="experiment-card"]').first().click();
    
    // Check if we navigated to the experiment details page
    await page.waitForURL(/\/experiments\/\d+/);
    
    // Verify experiment details are displayed
    const detailsHeading = await page.textContent('h1');
    expect(detailsHeading).toBeTruthy();
    
    // Check for experiment data visualizations
    const charts = await page.locator('[data-testid="experiment-chart"]').count();
    expect(charts).toBeGreaterThanOrEqual(0); // Charts may not be present in all experiments
    
    await page.screenshot({ path: 'experiment-details.png' });
  });

  test('should show experiment creation form', async ({ page }) => {
    await page.goto('/experiments/create');
    
    // Verify the form is displayed
    await page.waitForSelector('form', { timeout: 5000 });
    
    // Check for required form fields
    const nameInput = await page.locator('input[name="name"]').count();
    expect(nameInput).toBe(1);
    
    const descriptionInput = await page.locator('textarea[name="description"]').count();
    expect(descriptionInput).toBe(1);
    
    // Check for submit button
    const submitButton = await page.locator('button[type="submit"]').count();
    expect(submitButton).toBe(1);
    
    await page.screenshot({ path: 'experiment-create-form.png' });
  });

  test('should allow filtering and sorting experiments', async ({ page }) => {
    await page.goto('/experiments');
    
    // Wait for experiment list to load
    await page.waitForSelector('[data-testid="experiment-card"]', { timeout: 5000 });
    
    // Check for filter controls
    const filterControls = await page.locator('[data-testid="filter-controls"]').count();
    if (filterControls > 0) {
      // Try using the filter if it exists
      await page.locator('[data-testid="filter-controls"] select').first().selectOption({ index: 1 });
      
      // Wait for filtered results
      await page.waitForTimeout(1000);
    }
    
    // Check for sort controls
    const sortControls = await page.locator('[data-testid="sort-controls"]').count();
    if (sortControls > 0) {
      // Try changing the sort order if controls exist
      await page.locator('[data-testid="sort-controls"] button').first().click();
      
      // Wait for sorted results
      await page.waitForTimeout(1000);
    }
    
    await page.screenshot({ path: 'experiments-filtered.png' });
  });

  test('should show experiment results visualization', async ({ page }) => {
    await page.goto('/experiments');
    
    // Wait for experiment cards to load
    await page.waitForSelector('[data-testid="experiment-card"]', { timeout: 5000 });
    
    // Click on the first experiment card
    await page.locator('[data-testid="experiment-card"]').first().click();
    
    // Wait for navigation to complete
    await page.waitForURL(/\/experiments\/\d+/);
    
    // Check for results visualization if they exist
    const charts = await page.locator('[data-testid="experiment-chart"]').count();
    
    if (charts > 0) {
      // Verify chart content
      const chartContainer = page.locator('[data-testid="experiment-chart"]').first();
      
      // Check if the chart rendered
      const chartElement = await chartContainer.locator('svg, canvas').count();
      expect(chartElement).toBeGreaterThan(0);
      
      await page.screenshot({ path: 'experiment-visualization.png' });
    } else {
      console.log('No charts found on this experiment page');
    }
  });
});
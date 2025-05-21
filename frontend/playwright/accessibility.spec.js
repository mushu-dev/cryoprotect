// @ts-check
const { test, expect } = require('@playwright/test');
const AxeBuilder = require('@axe-core/playwright').default;

/**
 * Accessibility tests for the experimental data enhancement UI
 * 
 * These tests use Axe to check for accessibility issues on key pages
 */
test.describe('Accessibility Tests', () => {
  test('experiments list page should not have accessibility violations', async ({ page }) => {
    await page.goto('/experiments');
    
    // Wait for the page to fully load
    await page.waitForSelector('[data-testid="experiment-card"]', { timeout: 5000 });
    
    // Run axe analysis
    const accessibilityScanResults = await new AxeBuilder({ page })
      .withTags(['wcag2a', 'wcag2aa', 'wcag21a', 'wcag21aa'])
      .analyze();
    
    // Log any issues
    if (accessibilityScanResults.violations.length > 0) {
      console.log('Accessibility violations:', accessibilityScanResults.violations);
    }
    
    // Assert no violations
    expect(accessibilityScanResults.violations).toEqual([]);
  });

  test('experiment details page should not have accessibility violations', async ({ page }) => {
    // First go to experiments list
    await page.goto('/experiments');
    
    // Wait for experiment cards to load
    await page.waitForSelector('[data-testid="experiment-card"]', { timeout: 5000 });
    
    // Click on the first experiment if it exists
    const experimentCard = page.locator('[data-testid="experiment-card"]').first();
    if (await experimentCard.count() > 0) {
      await experimentCard.click();
      
      // Wait for navigation to experiment details
      await page.waitForURL(/\/experiments\/\d+/);
      
      // Run axe analysis
      const accessibilityScanResults = await new AxeBuilder({ page })
        .withTags(['wcag2a', 'wcag2aa', 'wcag21a', 'wcag21aa'])
        .analyze();
      
      // Log any issues
      if (accessibilityScanResults.violations.length > 0) {
        console.log('Accessibility violations:', accessibilityScanResults.violations);
      }
      
      // Assert no violations
      expect(accessibilityScanResults.violations).toEqual([]);
    } else {
      test.skip();
    }
  });

  test('experiment creation form should not have accessibility violations', async ({ page }) => {
    await page.goto('/experiments/create');
    
    // Wait for the form to load
    await page.waitForSelector('form', { timeout: 5000 });
    
    // Run axe analysis
    const accessibilityScanResults = await new AxeBuilder({ page })
      .withTags(['wcag2a', 'wcag2aa', 'wcag21a', 'wcag21aa'])
      .analyze();
    
    // Log any issues
    if (accessibilityScanResults.violations.length > 0) {
      console.log('Accessibility violations:', accessibilityScanResults.violations);
    }
    
    // Assert no violations
    expect(accessibilityScanResults.violations).toEqual([]);
  });

  test('data visualizations should have appropriate accessibility features', async ({ page }) => {
    // Go to experiments list
    await page.goto('/experiments');
    
    // Wait for experiment cards to load
    await page.waitForSelector('[data-testid="experiment-card"]', { timeout: 5000 });
    
    // Click on the first experiment if it exists
    const experimentCard = page.locator('[data-testid="experiment-card"]').first();
    if (await experimentCard.count() > 0) {
      await experimentCard.click();
      
      // Wait for navigation to experiment details
      await page.waitForURL(/\/experiments\/\d+/);
      
      // Check for data visualizations
      const charts = page.locator('[data-testid="experiment-chart"]');
      if (await charts.count() > 0) {
        // Check if the chart has aria attributes
        const hasAriaLabel = await charts.first().getAttribute('aria-label');
        const hasAriaDescription = await charts.first().getAttribute('aria-description');
        
        // Check for either aria-label or aria-description
        expect(hasAriaLabel || hasAriaDescription).toBeTruthy();
        
        // Check if there's an accompanying table with the data
        const dataTable = page.locator('table[data-testid="experiment-data-table"]');
        const hasDataTable = await dataTable.count() > 0;
        
        // Check for either a data table or aria attributes
        expect(hasDataTable || hasAriaLabel || hasAriaDescription).toBeTruthy();
      }
    } else {
      test.skip();
    }
  });

  test('keyboard navigation should work correctly', async ({ page }) => {
    await page.goto('/experiments');
    
    // Focus on the first interactive element
    await page.keyboard.press('Tab');
    
    // Get the active element
    const firstFocusedElement = await page.evaluate(() => document.activeElement.tagName);
    
    // Verify we have focus
    expect(firstFocusedElement).not.toBe('BODY');
    
    // Try tabbing through the page
    const focusableElements = [];
    for (let i = 0; i < 10; i++) {
      // Press Tab to move to the next element
      await page.keyboard.press('Tab');
      
      // Get the currently focused element
      const tagName = await page.evaluate(() => document.activeElement.tagName);
      const ariaLabel = await page.evaluate(() => document.activeElement.getAttribute('aria-label'));
      
      // Add to our list of elements we've focused
      focusableElements.push({ tagName, ariaLabel });
    }
    
    // Verify we were able to focus multiple elements
    const uniqueElements = new Set(focusableElements.map(e => e.tagName));
    expect(uniqueElements.size).toBeGreaterThan(1);
  });
});
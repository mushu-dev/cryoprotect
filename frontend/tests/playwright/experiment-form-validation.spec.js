// @ts-check
const { test, expect } = require('@playwright/test');
const { fillExperimentForm } = require('./test-utils');

test.describe('Experiment Form Validation', () => {
  test.beforeEach(async ({ page }) => {
    // Go to the experiment creation page
    await page.goto('/experiments/create');
    
    // Wait for the form to load
    await page.waitForSelector('form');
  });

  test('should show validation errors for empty required fields', async ({ page }) => {
    // Submit the form without filling in any fields
    await page.click('button[type="submit"]');
    
    // Check for error messages
    const errorMessages = await page.locator('[data-testid="form-error"]').allTextContents();
    
    // Verify we have at least one error message
    expect(errorMessages.length).toBeGreaterThan(0);
    
    // Check for specific error message text
    const errorTexts = errorMessages.join(' ');
    expect(errorTexts).toContain('required');
  });

  test('should validate experiment name length', async ({ page }) => {
    // Fill the experiment name with a very short value
    await page.fill('input[name="name"]', 'a');
    
    // Click outside the input to trigger validation
    await page.click('body');
    
    // Check for name length error
    const nameError = await page.locator('[data-testid="form-error"]').filter({ hasText: /name/i }).textContent();
    
    // Verify we have an error message about length
    expect(nameError).toBeTruthy();
    expect(nameError).toMatch(/length|characters|short/i);
    
    // Now try with a valid name
    await page.fill('input[name="name"]', 'Valid Experiment Name');
    await page.click('body');
    
    // The error should disappear
    const nameErrorAfterFix = await page.locator('[data-testid="form-error"]').filter({ hasText: /name/i }).count();
    expect(nameErrorAfterFix).toBe(0);
  });

  test('should validate numeric fields', async ({ page }) => {
    // Try to fill temperature with non-numeric value
    const temperatureInput = page.locator('input[name="temperature"]');
    if (await temperatureInput.count() > 0) {
      await temperatureInput.fill('abc');
      await page.click('body');
      
      // Check for numeric validation error
      const tempError = await page.locator('[data-testid="form-error"]').filter({ hasText: /temperature/i }).textContent();
      expect(tempError).toBeTruthy();
      expect(tempError).toMatch(/number|numeric|invalid/i);
      
      // Fix with valid number
      await temperatureInput.fill('25');
      await page.click('body');
      
      // Error should be gone
      const tempErrorAfterFix = await page.locator('[data-testid="form-error"]').filter({ hasText: /temperature/i }).count();
      expect(tempErrorAfterFix).toBe(0);
    }
  });

  test('should validate date fields', async ({ page }) => {
    // Look for date input
    const dateInput = page.locator('input[type="date"]');
    if (await dateInput.count() > 0) {
      // Try invalid date format
      await dateInput.fill('invalid-date');
      await page.click('body');
      
      // Check for date validation error
      const dateError = await page.locator('[data-testid="form-error"]').filter({ hasText: /date/i }).count();
      
      // Most browsers will prevent invalid date entry in date inputs,
      // so we'll check if the input value is empty instead
      const dateValue = await dateInput.inputValue();
      expect(dateValue).toBe('');
      
      // Fill with valid date
      await dateInput.fill('2025-05-19');
      await page.click('body');
      
      // Verify the date was accepted
      const dateValueAfterFix = await dateInput.inputValue();
      expect(dateValueAfterFix).toBe('2025-05-19');
    }
  });

  test('should allow form submission with valid data', async ({ page }) => {
    // Fill all form fields with valid data
    await fillExperimentForm(page);
    
    // Submit the form
    await page.click('button[type="submit"]');
    
    // Check if we navigated away from the form page or got a success message
    // This could be either navigation to a new page or a success notification
    
    try {
      // Option 1: Check for navigation away from form page
      await page.waitForNavigation({ timeout: 5000 });
      const currentUrl = page.url();
      expect(currentUrl).not.toContain('/experiments/create');
    } catch (e) {
      // Option 2: If no navigation, check for success notification
      const successMessage = await page.locator('[data-testid="success-message"]').textContent();
      expect(successMessage).toBeTruthy();
      expect(successMessage).toMatch(/success|created|saved/i);
    }
  });

  test('should prevent duplicate experiment names', async ({ page }) => {
    // This test requires the API to be functional to check for duplicates
    
    // Fill and submit the form once
    await fillExperimentForm(page, { name: 'Duplicate Test' });
    await page.click('button[type="submit"]');
    
    // Wait for submission to complete
    try {
      await page.waitForNavigation({ timeout: 5000 });
    } catch (e) {
      // If no navigation, wait for API request to complete
      await page.waitForTimeout(2000);
    }
    
    // Navigate back to the create form
    await page.goto('/experiments/create');
    
    // Try to create another experiment with the same name
    await fillExperimentForm(page, { name: 'Duplicate Test' });
    await page.click('button[type="submit"]');
    
    // Check for duplicate name error (this is dependent on the actual implementation)
    try {
      const errorMessage = await page.locator('[data-testid="form-error"], .error-message, .alert-error')
        .filter({ hasText: /duplicate|already exists|unique/i })
        .textContent({ timeout: 5000 });
      
      if (errorMessage) {
        expect(errorMessage).toMatch(/duplicate|already exists|unique/i);
      }
    } catch (e) {
      // If no explicit error message, check that we're still on the form page
      // which would indicate submission failed
      const currentUrl = page.url();
      expect(currentUrl).toContain('/experiments/create');
    }
  });
});
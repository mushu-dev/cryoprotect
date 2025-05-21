// @ts-check
const { test, expect } = require('@playwright/test');

test.describe('Convex Integration Tests', () => {
  test.beforeEach(async ({ page }) => {
    // Navigate to the Convex status page
    await page.goto('/convex-status');
    // Wait for a reasonable amount of time for Convex to connect
    await page.waitForTimeout(2000);
  });

  test('Convex status page loads correctly', async ({ page }) => {
    // Verify the header is visible
    await expect(page.locator('h1:has-text("Convex Database Status")')).toBeVisible();
    
    // Verify the status card is present
    await expect(page.locator('.card')).toBeVisible();
  });

  test('Convex connection status is shown', async ({ page }) => {
    // Look for status badge in the card
    const statusBadge = page.locator('div.card span:has-text("Connected"), div.card span:has-text("Error"), div.card span:has-text("Connecting")');
    await expect(statusBadge).toBeVisible();
    
    // Check that we have URL information
    await expect(page.locator('text=URL:')).toBeVisible();
    
    // Check that we have environment information
    await expect(page.locator('text=Environment:')).toBeVisible();
  });

  test('Navigation includes Convex status link', async ({ page }) => {
    // Go to the homepage
    await page.goto('/');
    
    // Check for Convex status link in the navigation
    const navLink = page.locator('nav >> text=Convex Status');
    
    // Since it might be conditionally shown based on environment variables
    const linkExists = await navLink.count() > 0;
    
    if (linkExists) {
      // If the link exists, it should be visible and clickable
      await expect(navLink).toBeVisible();
      await navLink.click();
      
      // After clicking, we should be on the Convex status page
      await expect(page.url()).toContain('/convex-status');
    } else {
      // If the link doesn't exist, log a message but don't fail the test
      console.log('Convex Status link not found in navigation - this is expected if NEXT_PUBLIC_USE_CONVEX is not set to true');
    }
  });

  test('Refresh connection button is functional', async ({ page }) => {
    // Find the refresh button
    const refreshButton = page.locator('button:has-text("Refresh Connection")');
    await expect(refreshButton).toBeVisible();
    
    // Click the button
    await refreshButton.click();
    
    // After clicking, page should reload and the button should still be there
    await page.waitForLoadState('networkidle');
    await expect(refreshButton).toBeVisible();
  });

  test('About Convex integration section is present', async ({ page }) => {
    // Check for the About section
    await expect(page.locator('h2:has-text("About Convex Integration")')).toBeVisible();
    
    // Check for key benefits
    await expect(page.locator('text=Real-time data synchronization')).toBeVisible();
    await expect(page.locator('text=TypeScript integration')).toBeVisible();
  });
});
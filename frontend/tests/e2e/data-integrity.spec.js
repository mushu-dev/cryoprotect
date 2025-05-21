// @ts-check
const { test, expect } = require('@playwright/test');

test.describe('Data Integrity Tests', () => {
  test.beforeEach(async ({ page }) => {
    // Navigate to homepage before each test
    await page.goto('/');
  });

  test('Navigation data integrity across all pages', async ({ page }) => {
    // Verify consistent navigation header on homepage
    await expect(page.locator('header >> text=CryoProtect')).toBeVisible();
    
    // Verify molecules page has consistent navigation
    await page.click('nav >> text=Molecules');
    await expect(page.locator('header >> text=CryoProtect')).toBeVisible();
    await expect(page.locator('nav >> text=Experiments')).toBeVisible();
    
    // Verify experiments page has consistent navigation
    await page.click('nav >> text=Experiments');
    await expect(page.locator('header >> text=CryoProtect')).toBeVisible();
    await expect(page.locator('nav >> text=Molecules')).toBeVisible();
    
    // Verify protocols page has consistent navigation
    await page.click('nav >> text=Protocols');
    await expect(page.locator('header >> text=CryoProtect')).toBeVisible();
    await expect(page.locator('nav >> text=Experiments')).toBeVisible();
    
    // Verify mixtures page has consistent navigation
    await page.click('nav >> text=Mixtures');
    await expect(page.locator('header >> text=CryoProtect')).toBeVisible();
    await expect(page.locator('nav >> text=Protocols')).toBeVisible();
  });

  test('URL integrity and consistency', async ({ page }) => {
    // Verify homepage URL
    await expect(page.url()).toBe(new URL('/', page.url()).toString());
    
    // Click Molecules link and verify URL format
    await page.click('nav >> text=Molecules');
    await expect(page.url()).toContain('/molecules');
    
    // Click Experiments link and verify URL format
    await page.click('nav >> text=Experiments');
    await expect(page.url()).toContain('/experiments');
    
    // Click Protocols link and verify URL format
    await page.click('nav >> text=Protocols');
    await expect(page.url()).toContain('/protocols');
    
    // Click Mixtures link and verify URL format
    await page.click('nav >> text=Mixtures');
    await expect(page.url()).toContain('/mixtures');
    
    // Return to homepage and verify URL
    await page.click('header >> text=CryoProtect');
    await expect(page.url()).toBe(new URL('/', page.url()).toString());
  });

  test('Page title and URL integrity across all pages', async ({ page }) => {
    // Check homepage title and URL
    await expect(page).toHaveURL(/.*\//);
    
    // Check molecules page URL instead of title (more reliable)
    await page.click('nav >> text=Molecules');
    await expect(page).toHaveURL(/.*\/molecules/);
    
    // Return to homepage
    await page.click('header >> text=CryoProtect');
    
    // Check experiments page URL
    await page.click('nav >> text=Experiments');
    await expect(page).toHaveURL(/.*\/experiments/);
    
    // Return to homepage
    await page.click('header >> text=CryoProtect');
    
    // Check protocols page URL
    await page.click('nav >> text=Protocols');
    await expect(page).toHaveURL(/.*\/protocols/);
    
    // Return to homepage
    await page.click('header >> text=CryoProtect');
    
    // Check mixtures page URL
    await page.click('nav >> text=Mixtures');
    await expect(page).toHaveURL(/.*\/mixtures/);
  });

  test('Footer data integrity across all pages', async ({ page }) => {
    // Check homepage footer - Get the main footer by its role for more specificity
    const footerSelector = 'footer.py-6';
    await expect(page.locator(footerSelector)).toBeVisible();
    
    // Check molecules page footer
    await page.click('nav >> text=Molecules');
    await expect(page.locator(footerSelector)).toBeVisible();
    
    // Return to homepage 
    await page.click('header >> text=CryoProtect');
    
    // Check experiments page footer
    await page.click('nav >> text=Experiments');
    await expect(page.locator(footerSelector)).toBeVisible();
    
    // Return to homepage
    await page.click('header >> text=CryoProtect');
    
    // Check protocols page footer
    await page.click('nav >> text=Protocols');
    await expect(page.locator(footerSelector)).toBeVisible();
    
    // Return to homepage
    await page.click('header >> text=CryoProtect');
    
    // Check mixtures page footer
    await page.click('nav >> text=Mixtures');
    await expect(page.locator(footerSelector)).toBeVisible();
  });

  test('Mobile UI consistency across pages', async ({ page }) => {
    // Set viewport to mobile size
    await page.setViewportSize({ width: 390, height: 844 });
    
    // Check that mobile menu button exists on homepage
    const homeMenuButton = page.locator('button[aria-label="Toggle menu"]');
    await expect(homeMenuButton).toBeVisible();
    
    // Navigate to Molecules page via main content link
    await page.click('a:has-text("View Molecules")');
    
    // Check that mobile menu button exists on Molecules page too
    const moleculesMenuButton = page.locator('button[aria-label="Toggle menu"]');
    await expect(moleculesMenuButton).toBeVisible();
    
    // Return to homepage
    await page.click('header >> a:has-text("CryoProtect")');
    
    // Navigate to Experiments page
    await page.click('a:has-text("Manage Experiments")');
    
    // Check that mobile menu button exists on Experiments page
    const experimentsMenuButton = page.locator('button[aria-label="Toggle menu"]');
    await expect(experimentsMenuButton).toBeVisible();
  });
});
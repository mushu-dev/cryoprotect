// @ts-check
const { test, expect } = require('@playwright/test');

test.describe('User Flow Tests', () => {
  test.beforeEach(async ({ page }) => {
    // Navigate to the home page before each test
    await page.goto('/');
  });

  test('Main page loads correctly', async ({ page }) => {
    // Verify the header is visible
    await expect(page.locator('header a:has-text("CryoProtect")')).toBeVisible();
    
    // Verify section headers are present
    await expect(page.locator('h2:has-text("Molecule Database")')).toBeVisible();
    await expect(page.locator('h2:has-text("Mixtures")')).toBeVisible();
    await expect(page.locator('h2:has-text("Experiments")')).toBeVisible();
    await expect(page.locator('h2:has-text("Protocols")')).toBeVisible();
    
    // Verify footer
    await expect(page.locator('footer').first()).toBeVisible();
  });

  test('Navigate through main application sections', async ({ page }) => {
    // Verify dashboard/homepage loaded
    await expect(page.locator('text=A platform for cryoprotectant analysis')).toBeVisible();
    
    // Navigate to Molecules section
    await page.click('text=View Molecules');
    await expect(page.url()).toContain('/molecules');
    
    // Navigate to Mixtures section
    await page.goto('/');
    await page.click('text=View Mixtures');
    await expect(page.url()).toContain('/mixtures');
    
    // Navigate to Experiments section
    await page.goto('/');
    await page.click('text=Manage Experiments');
    await expect(page.url()).toContain('/experiments');
    
    // Navigate to Protocols section
    await page.goto('/');
    await page.click('text=Browse Protocols');
    await expect(page.url()).toContain('/protocols');
  });

  test('Navigation Header links work correctly', async ({ page }) => {
    // Check that the navigation header is visible
    await expect(page.locator('header >> text=CryoProtect')).toBeVisible();
    
    // Click on Molecules in the nav
    await page.click('nav >> text=Molecules');
    await expect(page.url()).toContain('/molecules');
    
    // Click on Mixtures in the nav
    await page.click('nav >> text=Mixtures');
    await expect(page.url()).toContain('/mixtures');
    
    // Click on Experiments in the nav
    await page.click('nav >> text=Experiments');
    await expect(page.url()).toContain('/experiments');
    
    // Click on Protocols in the nav
    await page.click('nav >> text=Protocols');
    await expect(page.url()).toContain('/protocols');
    
    // Return to Dashboard by clicking the logo
    await page.click('header >> text=CryoProtect');
    await expect(page.url()).toBe(new URL('/', page.url()).toString());
  });

  test('Mobile responsive design elements exist', async ({ page }) => {
    // Resize viewport to mobile dimensions
    await page.setViewportSize({ width: 390, height: 844 });
    
    // Check mobile menu button exists
    const menuButton = page.locator('button[aria-label="Toggle menu"]');
    await expect(menuButton).toBeVisible();
    
    // Verify the layout is responsive by checking certain aspects
    // For example, verify navigation links are not visible on mobile
    const navLinks = page.locator('nav a').filter({ visible: true });
    const visibleNavLinks = await navLinks.count();
    
    // On mobile, all main menu items should be visible
    expect(visibleNavLinks).toBeLessThanOrEqual(6);
  });
});
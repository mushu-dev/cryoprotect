// @ts-check
const { test, expect } = require('@playwright/test');

test.describe('UI Navigation and Component Tests', () => {
  test('should navigate through pages and test UI components', async ({ page }) => {
    // Navigate to the home page
    await page.goto('/');
    
    // Verify home page title
    await expect(page.locator('h1').first()).toHaveText('CryoProtect');
    
    // Take screenshot of home page
    await page.screenshot({ path: 'test-results/home-page.png' });
    
    // Navigate to the dashboard page
    await page.click('text=Dashboard');
    
    // Verify we're on the dashboard page
    await expect(page).toHaveURL(/.*\/dashboard/);
    await expect(page.locator('h1').first()).toHaveText('Dashboard');
    
    // Check that sidebar navigation is visible
    const sidebar = page.locator('.fixed.inset-y-0.left-0.z-50.w-64');
    await expect(sidebar).toBeVisible();
    
    // Verify sidebar navigation items
    const navItems = ['Dashboard', 'Molecules', 'Mixtures', 'Experiments', 'Protocols', 'Properties'];
    for (const item of navItems) {
      await expect(sidebar.locator(`text=${item}`)).toBeVisible();
    }
    
    // Take screenshot of dashboard
    await page.screenshot({ path: 'test-results/dashboard.png' });
    
    // Navigate to the components page
    await page.goto('/components');
    
    // Verify components page title
    await expect(page.locator('h1').first()).toHaveText('Component Library');
    
    // Test basic components are visible
    await expect(page.locator('text=Buttons')).toBeVisible();
    await expect(page.locator('text=Cards')).toBeVisible();
    await expect(page.locator('text=Alerts')).toBeVisible();
    await expect(page.locator('text=Badges')).toBeVisible();
    
    // Test tabs interaction
    await page.locator('button:has-text("Experiments")').click();
    await expect(page.locator('text=Experiments Tab')).toBeVisible();
    
    // Test dialog interaction
    await page.locator('button:has-text("Open Dialog")').click();
    await expect(page.locator('text=Confirm Action')).toBeVisible();
    
    // Close dialog
    await page.locator('button:has-text("Cancel")').first().click();
    
    // Take screenshot of components page
    await page.screenshot({ path: 'test-results/components.png' });
    
    // Navigate to the experiments page
    await page.goto('/experiments');
    
    // Verify experiments page
    await expect(page).toHaveURL(/.*\/experiments/);
    
    // Take screenshot of experiments page
    await page.screenshot({ path: 'test-results/experiments.png' });
  });
});
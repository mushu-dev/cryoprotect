// @ts-check
const { test, expect } = require('@playwright/test');

/**
 * Basic navigation test for the CryoProtect application
 * 
 * This test checks basic navigation and UI components to ensure
 * the core user interface is functioning correctly.
 */
test.describe('Basic Navigation Test', () => {
  test('should navigate through main pages and check UI components', async ({ page }) => {
    // Step 1: Start the development server is handled by the webServer config in playwright.config.js
    console.log('Step 1: Development server should be started by Playwright config');
    
    // Step 2: Navigate to the home page
    console.log('Step 2: Navigating to home page');
    await page.goto('/');
    
    // Verify we're on the home page
    await expect(page.locator('h1')).toHaveText('CryoProtect');
    await expect(page.locator('p.text-xl')).toContainText('cryoprotectant analysis');
    
    // Take screenshot of home page
    await page.screenshot({ path: './test-results/basic-navigation/home-page.png' });
    
    // Step 3: Navigate to the dashboard page
    console.log('Step 3: Navigating to dashboard page');
    
    // Look for dashboard link - might be in navigation or on home page
    const dashboardLink = page.locator('a[href="/dashboard"]').first();
    await dashboardLink.click();
    
    // Verify we're on the dashboard page
    await page.waitForURL('**/dashboard');
    
    // Wait for dashboard content to load
    await page.waitForSelector('.container');
    
    // Take screenshot of dashboard page
    await page.screenshot({ path: './test-results/basic-navigation/dashboard-page.png' });
    
    // Step 4: Check that the navigation sidebar is visible
    console.log('Step 4: Checking navigation sidebar');
    
    // The sidebar should be visible on desktop sizes
    await expect(page.locator('aside nav')).toBeVisible();
    
    // Check if all main navigation items are present
    const navItems = [
      'Dashboard',
      'Molecules',
      'Mixtures',
      'Experiments',
      'Protocols',
      'Properties'
    ];
    
    for (const item of navItems) {
      await expect(page.locator(`nav a:has-text("${item}")`)).toBeVisible();
    }
    
    // Take a screenshot specifically of the sidebar
    await page.locator('aside').screenshot({ 
      path: './test-results/basic-navigation/sidebar.png'
    });
    
    // Step 5: Navigate to the components page
    console.log('Step 5: Navigating to components page');
    
    // Direct navigation since components may not be in the main nav
    await page.goto('/components');
    
    // Wait for the components page to load
    await page.waitForSelector('h1:has-text("Component Library")');
    
    // Take screenshot of components page
    await page.screenshot({ 
      path: './test-results/basic-navigation/components-page.png',
      fullPage: true 
    });
    
    // Step 6: Check that the components are rendered correctly
    console.log('Step 6: Checking component rendering');
    
    // Verify essential component sections exist
    const componentSections = [
      'Buttons',
      'Cards',
      'Alerts',
      'Badges',
      'Tabs',
      'Form Controls'
    ];
    
    for (const section of componentSections) {
      const sectionElement = page.locator(`h2:has-text("${section}")`);
      await expect(sectionElement).toBeVisible();
      
      // Scroll to the section for visibility
      await sectionElement.scrollIntoViewIfNeeded();
      
      // Take a screenshot of each section
      const sectionCard = page.locator(`h2:has-text("${section}")`)
        .locator('xpath=../..')
        .locator('xpath=..');
        
      await sectionCard.screenshot({ 
        path: `./test-results/basic-navigation/component-${section.toLowerCase().replace(/\s+/g, '-')}.png` 
      });
    }
    
    // Test component interactivity - try clicking a tab
    const tabsSection = page.locator('h2:has-text("Tabs")').locator('xpath=../..').locator('xpath=..');
    
    if (await tabsSection.isVisible()) {
      await tabsSection.locator('button:has-text("Experiments")').click();
      
      // Verify the tab content changed
      await expect(tabsSection.locator('div:has-text("Experiments Tab")')).toBeVisible();
      
      // Take screenshot of tab interaction
      await tabsSection.screenshot({ 
        path: './test-results/basic-navigation/tab-interaction.png' 
      });
    }
    
    // Test button interactivity
    const buttonsSection = page.locator('h2:has-text("Buttons")').locator('xpath=../..').locator('xpath=..');
    
    if (await buttonsSection.isVisible()) {
      // Click a button (this won't cause navigation)
      await buttonsSection.locator('button:has-text("Default Button")').click();
      
      // Take screenshot after clicking
      await buttonsSection.screenshot({ 
        path: './test-results/basic-navigation/button-interaction.png' 
      });
    }
    
    // Try opening a dialog
    const dialogsSection = page.locator('h2:has-text("Dialogs")').locator('xpath=../..').locator('xpath=..');
    
    if (await dialogsSection.isVisible()) {
      await dialogsSection.locator('button:has-text("Open Dialog")').click();
      
      // Verify dialog opens
      await expect(page.locator('div[role="dialog"]')).toBeVisible();
      
      // Take screenshot of open dialog
      await page.screenshot({ 
        path: './test-results/basic-navigation/open-dialog.png' 
      });
      
      // Close the dialog by clicking cancel
      await page.locator('div[role="dialog"]').locator('button:has-text("Cancel")').click();
    }
    
    // Step 7: Navigate back to home page to complete the test
    console.log('Step 7: Returning to home page');
    
    // Navigate home by clicking the logo
    await page.locator('a:has-text("CryoProtect")').first().click();
    
    // Verify we're back on the home page
    await expect(page.locator('h1')).toHaveText('CryoProtect');
    
    // Take final screenshot
    await page.screenshot({ path: './test-results/basic-navigation/final-state.png' });
    
    console.log('Basic navigation test completed successfully');
  });
});
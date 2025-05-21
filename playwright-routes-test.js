// @ts-check
const { test, expect } = require('@playwright/test');

/**
 * Comprehensive test to verify all routes work correctly in the CryoProtect application
 */
test('CryoProtect Routes Verification', async ({ page }) => {
  console.log('Starting route verification test');
  
  // Define all routes to test
  const routes = [
    { path: '/', title: 'CryoProtect - A platform for cryoprotectant analysis', expected: 'A platform for cryoprotectant analysis' },
    { path: '/molecules/', title: 'Molecules - CryoProtect', expected: 'Molecules' },
    { path: '/mixtures/', title: 'Mixtures - CryoProtect', expected: 'Mixtures' },
    { path: '/experiments/', title: 'Experiments | CryoProtect', expected: 'Experiments' },
    { path: '/protocols/', title: 'Protocols | CryoProtect', expected: 'Protocols' },
    { path: '/properties/', title: 'Properties - CryoProtect', expected: 'Properties' },
  ];

  // Visit each route and verify content
  for (const route of routes) {
    console.log(`Testing route: ${route.path}`);
    
    // Navigate to the route
    await page.goto(`https://cryoprotect.app${route.path}`);
    
    // Verify page title
    await expect(page).toHaveTitle(route.title);
    
    // Verify page content contains expected heading
    const heading = page.locator('h1');
    await expect(heading).toContainText(route.expected);
    
    // Verify navigation highlighting is correct
    if (route.path !== '/') {
      const activeLink = page.locator(`nav a[href="${route.path}"]`);
      const classes = await activeLink.getAttribute('class');
      expect(classes).toContain('text-black');
    }
    
    // Take a screenshot for visual verification
    await page.screenshot({ path: `route-test-${route.path.replace(/\//g, '-')}.png` });
    
    console.log(`✅ Route ${route.path} verified successfully`);
  }

  // Test navigation flow
  console.log('Testing navigation flow between pages');
  
  // Start at home page
  await page.goto('https://cryoprotect.app/');
  
  // Click each navigation link and verify we reach the correct page
  for (let i = 1; i < routes.length; i++) {
    const route = routes[i];
    console.log(`Testing navigation to: ${route.path}`);
    
    // Click the navigation link
    await page.click(`nav a[href="${route.path}"]`);
    
    // Verify we reached the correct page
    await expect(page).toHaveURL(new RegExp(`.*${route.path}`));
    await expect(page).toHaveTitle(route.title);
    
    console.log(`✅ Navigation to ${route.path} verified successfully`);
  }

  // Test for 404 handling
  console.log('Testing 404 handling');
  await page.goto('https://cryoprotect.app/non-existent-page');
  await expect(page.locator('h1')).toContainText('404');
  await page.screenshot({ path: 'route-test-404.png' });

  console.log('All routes verified successfully!');
});
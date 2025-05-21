// @ts-check
const { test, expect } = require('@playwright/test');

test('CryoProtect Basic Navigation Test', async ({ page }) => {
  // Go to the homepage
  await page.goto('https://cryoprotect.app/');
  
  // Verify the page title
  const title = await page.title();
  expect(title).toBe('CryoProtect - A platform for cryoprotectant analysis');
  
  // Verify the header is present
  const header = await page.locator('h1:has-text("CryoProtect")');
  await expect(header).toBeVisible();
  
  // Verify navigation links are present
  const navLinks = await page.locator('a[href]').filter({
    hasText: /Dashboard|Molecules|Mixtures|Experiments|Protocols/
  });
  
  const count = await navLinks.count();
  expect(count).toBeGreaterThanOrEqual(5);
  
  // Check the API endpoint
  await page.goto('https://cryoprotect.app/api');
  const apiResponse = await page.locator('body').textContent();
  expect(apiResponse).toContain('CryoProtect API is running');
  
  // Take screenshots
  await page.goto('https://cryoprotect.app/');
  await page.screenshot({ path: 'homepage.png' });
  
  await page.goto('https://cryoprotect.app/api');
  await page.screenshot({ path: 'api-endpoint.png' });
  
  console.log('Test completed successfully!');
});
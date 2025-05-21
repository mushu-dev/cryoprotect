# Testing Best Practices for Experimental Data Enhancement UI

This document provides guidelines and best practices for testing the experimental data enhancement UI. These practices will help ensure a high-quality, reliable, and maintainable testing suite.

## General Testing Principles

1. **Test Pyramid**: Follow the test pyramid principle:
   - Many unit tests
   - Moderate number of integration tests
   - Few end-to-end tests

2. **DRY (Don't Repeat Yourself)**: Use utility functions in `test-utils.js` to avoid code duplication.

3. **Independent Tests**: Tests should be independent and not rely on the state from other tests.

4. **Descriptive Names**: Use descriptive test names that clearly state what is being tested.

5. **Arrange-Act-Assert**: Structure tests following the AAA pattern:
   - Arrange: Set up test data and conditions
   - Act: Perform the action being tested
   - Assert: Verify the expected result

## End-to-End Testing with Playwright

### Page Object Pattern

Consider using the Page Object pattern for complex UIs to encapsulate page interactions:

```javascript
// experimentPage.js
class ExperimentPage {
  constructor(page) {
    this.page = page;
    this.cardSelector = '[data-testid="experiment-card"]';
    this.titleSelector = 'h1';
  }

  async navigate() {
    await this.page.goto('/experiments');
    await this.page.waitForSelector(this.cardSelector);
  }

  async getCardCount() {
    return await this.page.locator(this.cardSelector).count();
  }

  async clickFirstCard() {
    await this.page.locator(this.cardSelector).first().click();
  }
}

module.exports = ExperimentPage;
```

### Selecting Elements

1. **Prefer data-testid**: Use data-testid attributes for selecting elements:
   ```html
   <div data-testid="experiment-card">...</div>
   ```

2. **Avoid brittle selectors**: Avoid selectors that could easily change:
   ```javascript
   // Avoid
   await page.locator('.card-container > div:nth-child(2)').click();
   
   // Prefer
   await page.locator('[data-testid="experiment-card"]').click();
   ```

3. **Wait for elements properly**:
   ```javascript
   // Wait for element to be visible
   await page.waitForSelector('[data-testid="experiment-card"]', { state: 'visible' });
   
   // Wait for network to be idle
   await page.waitForLoadState('networkidle');
   ```

### Test Data Management

1. **Mock API responses** for predictable tests:
   ```javascript
   await page.route('**/api/experiments', route => {
     route.fulfill({
       status: 200,
       contentType: 'application/json',
       body: JSON.stringify(mockExperiments)
     });
   });
   ```

2. **Reset state** before each test:
   ```javascript
   test.beforeEach(async ({ page }) => {
     // Reset application state
     await page.evaluate(() => localStorage.clear());
     await page.evaluate(() => sessionStorage.clear());
   });
   ```

### Visual Testing

1. **Ensure consistent environments**: Use consistent browser dimensions, themes, etc.

2. **Version control screenshots** when they represent the expected state.

3. **Focus on critical UI components** rather than capturing entire pages.

## Performance Testing

### Measuring Performance

1. **Define clear metrics**:
   - Time to first paint
   - Time to interactive
   - DOM content loaded time
   - Full page load time

2. **Set explicit performance budgets**:
   ```javascript
   expect(loadTime).toBeLessThan(3000); // Page should load in less than 3 seconds
   ```

3. **Test with realistic data volumes** to catch scale-related issues.

## Security Testing

1. **Input validation**: Test with malicious inputs like XSS payloads.

2. **Check HTTP headers** for security best practices.

3. **Verify authentication flows** with both valid and invalid credentials.

4. **Monitor for information leakage** in responses and error messages.

## Accessibility Testing

1. **Use axe-core for automated tests**:
   ```javascript
   const accessibilityScanResults = await new AxeBuilder({ page }).analyze();
   expect(accessibilityScanResults.violations).toEqual([]);
   ```

2. **Test keyboard navigation** for all interactive elements.

3. **Check proper ARIA attributes** on custom components.

4. **Verify color contrast** meets WCAG standards.

## Continuous Integration

1. **Run critical tests on every commit**.

2. **Parallelize tests** to reduce CI time.

3. **Preserve artifacts** (screenshots, reports) from failed tests.

4. **Fail the build** on critical test failures.

## Debugging Tips

1. **Use Playwright UI mode** for interactive debugging:
   ```bash
   npm run test:e2e:ui
   ```

2. **Use debug mode** to step through tests:
   ```bash
   npm run test:e2e:debug
   ```

3. **Take screenshots** at critical points:
   ```javascript
   await page.screenshot({ path: 'debug-screenshot.png' });
   ```

4. **Check browser console logs**:
   ```javascript
   page.on('console', msg => console.log(`Browser console: ${msg.text()}`));
   ```

## Adding New Tests

1. **Start with a test plan**:
   - What needs to be tested?
   - What assertions should be made?
   - What data is needed?

2. **Write tests from the user's perspective**.

3. **Focus on behavior, not implementation details**.

4. **Add new utility functions** to `test-utils.js` for reusable functionality.

5. **Update this documentation** with new best practices as they emerge.

## Common Pitfalls to Avoid

1. **Flaky tests** due to:
   - Inadequate waiting for elements or network
   - Random test data
   - Race conditions

2. **Slow tests** due to:
   - Unnecessary browser operations
   - Not using test isolation
   - Running tests sequentially that could be parallel

3. **Overlapping tests** that create conflicts with each other

4. **Brittle selectors** that break with minor UI changes

5. **Missing assertions** that don't actually verify functionality

## Resources

- [Playwright Documentation](https://playwright.dev/docs/intro)
- [Web Accessibility Initiative (WAI)](https://www.w3.org/WAI/)
- [OWASP Web Security Testing Guide](https://owasp.org/www-project-web-security-testing-guide/)
- [Web Vitals](https://web.dev/vitals/) - Essential metrics for web performance
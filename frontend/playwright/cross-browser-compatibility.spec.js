// @ts-check
const { test, expect } = require('@playwright/test');

/**
 * Cross-browser compatibility tests for the experimental data enhancement UI
 * 
 * These tests verify that the UI works correctly across different browsers
 * and platforms. Use the --project=chromium,firefox,webkit flag to run these
 * tests on all three browser engines.
 */

// Run these tests on all configured browsers
test.describe('Cross-Browser Compatibility', () => {
  test('should render the experiments page correctly', async ({ page, browserName }) => {
    test.info().annotations.push({
      type: 'browser',
      description: browserName
    });
    
    await page.goto('/experiments');
    
    // Wait for experiments to load
    await page.waitForSelector('[data-testid="experiment-card"]', { timeout: 5000 });
    
    // Take a screenshot with browser name in the filename
    await page.screenshot({ 
      path: `./test-results/cross-browser/${browserName}-experiments-page.png`,
      fullPage: true 
    });
    
    // Basic UI verification
    expect(await page.locator('h1').textContent()).toContain('Experiment');
    expect(await page.locator('[data-testid="experiment-card"]').count()).toBeGreaterThan(0);
    
    // Check for browser-specific styling issues
    const hasLayoutIssues = await page.evaluate(() => {
      // Check for common layout problems
      const problems = [];
      
      // Check for overlapping elements
      document.querySelectorAll('[data-testid="experiment-card"]').forEach(card => {
        const rect = card.getBoundingClientRect();
        const overlaps = Array.from(document.querySelectorAll('[data-testid="experiment-card"]'))
          .filter(otherCard => card !== otherCard)
          .some(otherCard => {
            const otherRect = otherCard.getBoundingClientRect();
            return !(rect.right < otherRect.left || 
                     rect.left > otherRect.right || 
                     rect.bottom < otherRect.top || 
                     rect.top > otherRect.bottom);
          });
        
        if (overlaps) {
          problems.push('Overlapping cards detected');
        }
      });
      
      // Check for unexpected scrollbars
      const hasHorizontalScrollbar = document.body.scrollWidth > window.innerWidth;
      if (hasHorizontalScrollbar) {
        problems.push('Unexpected horizontal scrollbar');
      }
      
      // Check for elements that extend beyond viewport
      document.querySelectorAll('button, a, input, select').forEach(element => {
        const rect = element.getBoundingClientRect();
        if (rect.left < 0 || rect.right > window.innerWidth) {
          problems.push('Element extends beyond viewport width');
        }
      });
      
      return problems;
    });
    
    console.log(`Layout issues on ${browserName}:`, hasLayoutIssues);
    expect(hasLayoutIssues.length).toBe(0);
  });

  test('should render data visualizations correctly', async ({ page, browserName }) => {
    test.info().annotations.push({
      type: 'browser',
      description: browserName
    });
    
    await page.goto('/experiments');
    
    // Wait for experiments to load
    await page.waitForSelector('[data-testid="experiment-card"]', { timeout: 5000 });
    
    // Click on the first experiment card
    const experimentCard = page.locator('[data-testid="experiment-card"]').first();
    if (await experimentCard.count() > 0) {
      await experimentCard.click();
      
      // Wait for navigation
      await page.waitForURL(/\/experiments\/\d+/);
      
      // Check for data visualizations
      const charts = page.locator('[data-testid="experiment-chart"]');
      if (await charts.count() > 0) {
        // Take screenshot of visualization
        await charts.first().screenshot({ 
          path: `./test-results/cross-browser/${browserName}-experiment-chart.png` 
        });
        
        // Check if chart renders correctly
        const chartRenderedCorrectly = await page.evaluate(() => {
          const chartElement = document.querySelector('[data-testid="experiment-chart"]');
          if (!chartElement) return false;
          
          // Check if chart has minimum dimensions
          const rect = chartElement.getBoundingClientRect();
          if (rect.width < 50 || rect.height < 50) return false;
          
          // Check if chart has SVG or Canvas element
          const hasSvg = chartElement.querySelector('svg') !== null;
          const hasCanvas = chartElement.querySelector('canvas') !== null;
          
          return hasSvg || hasCanvas;
        });
        
        expect(chartRenderedCorrectly).toBeTruthy();
      }
    }
  });

  test('should handle forms correctly', async ({ page, browserName }) => {
    test.info().annotations.push({
      type: 'browser',
      description: browserName
    });
    
    await page.goto('/experiments/create');
    
    // Wait for form to load
    await page.waitForSelector('form', { timeout: 5000 });
    
    // Take screenshot of empty form
    await page.screenshot({ 
      path: `./test-results/cross-browser/${browserName}-experiment-form.png` 
    });
    
    // Fill form fields
    await page.fill('input[name="name"]', `Test in ${browserName}`);
    await page.fill('textarea[name="description"]', `Testing form compatibility in ${browserName}`);
    
    // Fill date if present
    const dateInput = page.locator('input[type="date"]');
    if (await dateInput.count() > 0) {
      await dateInput.fill('2025-05-19');
    }
    
    // Fill numeric fields if present
    const numericInputs = page.locator('input[type="number"]');
    for (let i = 0; i < await numericInputs.count(); i++) {
      await numericInputs.nth(i).fill('25');
    }
    
    // Take screenshot of filled form
    await page.screenshot({ 
      path: `./test-results/cross-browser/${browserName}-experiment-form-filled.png` 
    });
    
    // Check form validation
    const submitButton = page.locator('button[type="submit"]');
    expect(await submitButton.isEnabled()).toBeTruthy();
    
    // Check browser-specific form behavior
    const formElementsAccessible = await page.evaluate(() => {
      const formElements = document.querySelectorAll('form input, form textarea, form select, form button');
      
      // Check that all form elements are properly keyboard-focusable
      const allFocusable = Array.from(formElements).every(el => {
        // Try focusing the element
        el.focus();
        return document.activeElement === el;
      });
      
      return allFocusable;
    });
    
    expect(formElementsAccessible).toBeTruthy();
  });

  test('should handle responsive layouts correctly', async ({ page, browserName }) => {
    test.info().annotations.push({
      type: 'browser',
      description: browserName
    });
    
    // Test mobile viewport
    await page.setViewportSize({ width:
375, height: 667 }); // iPhone SE size
    await page.goto('/experiments');
    await page.waitForLoadState('networkidle');
    
    // Take screenshot of mobile view
    await page.screenshot({ 
      path: `./test-results/cross-browser/${browserName}-mobile.png`,
      fullPage: true
    });
    
    // Check mobile-specific UI elements
    const hasMobileNav = await page.evaluate(() => {
      // Check for hamburger menu or other mobile navigation indicators
      return document.querySelector('[data-testid="mobile-nav"]') !== null ||
             document.querySelector('.hamburger-menu') !== null ||
             document.querySelector('[aria-label="Menu"]') !== null;
    });
    
    // Not all UIs will have a mobile-specific nav, so just log this
    console.log(`Mobile navigation in ${browserName}:`, hasMobileNav ? 'Present' : 'Not found');
    
    // Check for layout issues in mobile viewport
    const mobileLayoutIssues = await page.evaluate(() => {
      const problems = [];
      
      // Check for horizontal overflow
      if (document.body.scrollWidth > window.innerWidth + 5) {
        problems.push('Horizontal overflow detected in mobile view');
      }
      
      // Check for tiny tap targets (accessibility issue)
      const smallTapTargets = Array.from(document.querySelectorAll('button, a, [role="button"]'))
        .filter(el => {
          const rect = el.getBoundingClientRect();
          // Minimum recommended tap target size is 44x44px
          return (rect.width < 44 || rect.height < 44) && rect.width > 0 && rect.height > 0;
        }).length;
      
      if (smallTapTargets > 0) {
        problems.push(`${smallTapTargets} small tap targets detected in mobile view`);
      }
      
      return problems;
    });
    
    console.log(`Mobile layout issues in ${browserName}:`, mobileLayoutIssues);
    
    // Test tablet viewport
    await page.setViewportSize({ width: 768, height: 1024 }); // iPad size
    await page.goto('/experiments');
    await page.waitForLoadState('networkidle');
    
    // Take screenshot of tablet view
    await page.screenshot({ 
      path: `./test-results/cross-browser/${browserName}-tablet.png`,
      fullPage: true
    });
    
    // Test desktop viewport
    await page.setViewportSize({ width: 1440, height: 900 }); // Desktop size
    await page.goto('/experiments');
    await page.waitForLoadState('networkidle');
    
    // Take screenshot of desktop view
    await page.screenshot({ 
      path: `./test-results/cross-browser/${browserName}-desktop.png`,
      fullPage: true
    });
  });

  test('should handle browser-specific features correctly', async ({ page, browserName }) => {
    test.info().annotations.push({
      type: 'browser',
      description: browserName
    });
    
    await page.goto('/experiments');
    
    // Test CSS features with browser-specific behavior
    const cssCompatibility = await page.evaluate(() => {
      const testElement = document.createElement('div');
      document.body.appendChild(testElement);
      
      // Test flexbox
      testElement.style.display = 'flex';
      const flexSupported = getComputedStyle(testElement).display === 'flex';
      
      // Test grid
      testElement.style.display = 'grid';
      const gridSupported = getComputedStyle(testElement).display === 'grid';
      
      // Test custom properties (CSS variables)
      testElement.style.setProperty('--test-var', 'red');
      testElement.style.color = 'var(--test-var)';
      const customPropertiesSupported = getComputedStyle(testElement).color !== '';
      
      // Clean up
      document.body.removeChild(testElement);
      
      return {
        flexSupported,
        gridSupported,
        customPropertiesSupported
      };
    });
    
    console.log(`CSS compatibility in ${browserName}:`, cssCompatibility);
    
    // All modern browsers should support these features
    expect(cssCompatibility.flexSupported).toBeTruthy();
    expect(cssCompatibility.gridSupported).toBeTruthy();
    expect(cssCompatibility.customPropertiesSupported).toBeTruthy();
    
    // Test browser-specific JavaScript APIs
    const jsCompatibility = await page.evaluate(() => {
      return {
        // DOM APIs
        intersectionObserver: 'IntersectionObserver' in window,
        resizeObserver: 'ResizeObserver' in window,
        mutationObserver: 'MutationObserver' in window,
        
        // ES6+ features
        promise: typeof Promise !== 'undefined',
        asyncAwait: (function() {
          try {
            eval('async function test() { await Promise.resolve(); }');
            return true;
          } catch (e) {
            return false;
          }
        })(),
        
        // Web APIs
        localStorage: typeof localStorage !== 'undefined',
        sessionStorage: typeof sessionStorage !== 'undefined',
        fetch: typeof fetch !== 'undefined',
      };
    });
    
    console.log(`JS compatibility in ${browserName}:`, jsCompatibility);
    
    // All modern browsers should support these features
    expect(jsCompatibility.promise).toBeTruthy();
    expect(jsCompatibility.localStorage).toBeTruthy();
    expect(jsCompatibility.fetch).toBeTruthy();
  });
});
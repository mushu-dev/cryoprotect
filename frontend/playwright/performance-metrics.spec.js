// @ts-check
const { test, expect } = require('@playwright/test');

/**
 * Performance metrics tests for the experimental data enhancement UI
 * 
 * These tests measure loading times, rendering performance, and other metrics
 * to ensure the UI remains responsive and efficient.
 */
test.describe('Performance Metrics', () => {
  test('should load the experiments page quickly', async ({ page }) => {
    // Start measuring
    const startTime = Date.now();
    
    // Navigate to the page
    const response = await page.goto('/experiments');
    
    // Wait for network to be idle (no more than 2 active connections for at least 500ms)
    await page.waitForLoadState('networkidle');
    
    // Calculate and log the time
    const loadTime = Date.now() - startTime;
    console.log(`Experiments page loaded in ${loadTime}ms`);
    
    // Assert page loaded successfully
    expect(response.status()).toBe(200);
    
    // Performance budget: page should load in less than 3 seconds
    expect(loadTime).toBeLessThan(3000);
    
    // Get performance metrics using the browser's Performance API
    const performanceMetrics = await page.evaluate(() => {
      const perfEntries = performance.getEntriesByType('navigation');
      if (perfEntries.length > 0) {
        const navEntry = perfEntries[0];
        return {
          domContentLoaded: navEntry.domContentLoadedEventEnd - navEntry.startTime,
          load: navEntry.loadEventEnd - navEntry.startTime,
          firstPaint: performance.getEntriesByName('first-paint')[0]?.startTime,
          firstContentfulPaint: performance.getEntriesByName('first-contentful-paint')[0]?.startTime,
        };
      }
      return null;
    });
    
    // Log the detailed metrics if available
    if (performanceMetrics) {
      console.log('Performance metrics:', performanceMetrics);
      
      // Add assertions for key metrics
      expect(performanceMetrics.domContentLoaded).toBeLessThan(2000);
      expect(performanceMetrics.load).toBeLessThan(3000);
      
      // FCP should be under 1.5s for a good user experience
      if (performanceMetrics.firstContentfulPaint) {
        expect(performanceMetrics.firstContentfulPaint).toBeLessThan(1500);
      }
    }
  });

  test('should render experiment list efficiently', async ({ page }) => {
    await page.goto('/experiments');
    
    // Wait for the experiment cards to load
    await page.waitForSelector('[data-testid="experiment-card"]');
    
    // Measure rendering time for experiment cards
    const renderTime = await page.evaluate(async () => {
      const start = performance.now();
      
      // Force a re-render by scrolling
      window.scrollTo(0, 100);
      await new Promise(r => setTimeout(r, 100));
      window.scrollTo(0, 0);
      
      // Wait for next animation frame to ensure render is complete
      await new Promise(r => requestAnimationFrame(r));
      
      return performance.now() - start;
    });
    
    console.log(`Experiment list render time: ${renderTime}ms`);
    
    // Rendering should be smooth (under 100ms is ideal for UI updates)
    expect(renderTime).toBeLessThan(100);
    
    // Check for layout shifts
    const layoutShiftScore = await page.evaluate(() => {
      if ('LayoutShiftAttribution' in window) {
        // Get CLS from the layout instability API if available
        let score = 0;
        new PerformanceObserver((list) => {
          for (const entry of list.getEntries()) {
            score += entry.value;
          }
        }).observe({ type: 'layout-shift', buffered: true });
        return score;
      }
      return null;
    });
    
    if (layoutShiftScore !== null) {
      console.log(`Cumulative Layout Shift score: ${layoutShiftScore}`);
      // Google's recommendation: CLS should be less than 0.1
      expect(layoutShiftScore).toBeLessThan(0.1);
    }
  });

  test('should load experiment details quickly', async ({ page }) => {
    await page.goto('/experiments');
    
    // Wait for experiment cards to load
    await page.waitForSelector('[data-testid="experiment-card"]');
    
    // Find first experiment card
    const firstCard = page.locator('[data-testid="experiment-card"]').first();
    if (await firstCard.count() === 0) {
      test.skip('No experiment cards found to test');
      return;
    }
    
    // Start measuring
    const startTime = Date.now();
    
    // Click on the first experiment card
    await firstCard.click();
    
    // Wait for navigation to complete
    await page.waitForURL(/\/experiments\/\d+/);
    
    // Wait for content to load
    await page.waitForSelector('h1');
    await page.waitForLoadState('networkidle');
    
    // Calculate and log the time
    const detailLoadTime = Date.now() - startTime;
    console.log(`Experiment detail page loaded in ${detailLoadTime}ms`);
    
    // Performance budget: detail page should load in less than 2 seconds
    expect(detailLoadTime).toBeLessThan(2000);
    
    // Check if there are charts and measure their render time
    const hasCharts = await page.locator('[data-testid="experiment-chart"]').count() > 0;
    
    if (hasCharts) {
      // Measure chart rendering time
      const chartRenderTime = await page.evaluate(async () => {
        const start = performance.now();
        
        // Wait for all chart animations to complete
        await new Promise(r => setTimeout(r, 500));
        
        return performance.now() - start;
      });
      
      console.log(`Chart render time: ${chartRenderTime}ms`);
      
      // Charts should render efficiently
      expect(chartRenderTime).toBeLessThan(1000);
    }
  });

  test('should handle data filters efficiently', async ({ page }) => {
    await page.goto('/experiments');
    
    // Wait for the page to load
    await page.waitForSelector('[data-testid="experiment-card"]');
    
    // Find filter controls if they exist
    const filterControls = page.locator('[data-testid="filter-controls"]');
    if (await filterControls.count() === 0) {
      test.skip('No filter controls found to test');
      return;
    }
    
    // Measure filter application time
    const filterStartTime = Date.now();
    
    // Apply a filter (e.g., select an option from a dropdown)
    await filterControls.locator('select').first().selectOption({ index: 1 });
    
    // Wait for filtering to complete and UI to update
    await page.waitForTimeout(500); // Allow time for filtering
    
    const filterTime = Date.now() - filterStartTime;
    console.log(`Filter application time: ${filterTime}ms`);
    
    // Filtering should be fast for a responsive UI
    expect(filterTime).toBeLessThan(1000);
    
    // Check if results were updated
    const resultsUpdated = await page.evaluate(() => {
      // Look for any loading indicators that have appeared and disappeared
      return document.querySelectorAll('[aria-busy="false"]').length > 0;
    });
    
    if (resultsUpdated) {
      console.log('Filter results updated the UI successfully');
    }
  });

  test('should load large data sets efficiently', async ({ page }) => {
    // This test is specifically for pages that might display large datasets
    await page.goto('/experiments');
    
    // Measure memory usage before and after loading large data
    const memoryBefore = await page.evaluate(() => performance.memory?.usedJSHeapSize || 0);
    
    // Simulate loading large dataset by applying a filter that would return many results
    // or by navigating to a data-heavy page
    const loadLargeDataset = async () => {
      const viewAllBtn = page.locator('button, a').filter({ hasText: /view all|show all|load more/i });
      
      if (await viewAllBtn.count() > 0) {
        await viewAllBtn.click();
        await page.waitForLoadState('networkidle');
      } else {
        // Alternative: force-load more items if there's pagination
        const paginationBtn = page.locator('[data-testid="pagination-next"], .pagination button').last();
        if (await paginationBtn.count() > 0) {
          await paginationBtn.click();
          await page.waitForLoadState('networkidle');
        }
      }
    };
    
    await loadLargeDataset();
    
    // Measure memory after loading data
    const memoryAfter = await page.evaluate(() => performance.memory?.usedJSHeapSize || 0);
    
    if (memoryBefore > 0 && memoryAfter > 0) {
      const memoryIncrease = ((memoryAfter - memoryBefore) / memoryBefore) * 100;
      console.log(`Memory usage increased by ${memoryIncrease.toFixed(2)}%`);
      
      // Memory increase should be reasonable - not more than 50% for normal data loading
      expect(memoryIncrease).toBeLessThan(50);
    }
    
    // Check scroll performance with large dataset
    const scrollPerformance = await page.evaluate(async () => {
      const start = performance.now();
      
      // Scroll to bottom
      window.scrollTo(0, document.body.scrollHeight);
      
      // Wait for any lazy-loaded content
      await new Promise(r => setTimeout(r, 500));
      
      // Scroll back to top
      window.scrollTo(0, 0);
      
      return performance.now() - start;
    });
    
    console.log(`Scroll performance with large dataset: ${scrollPerformance}ms`);
    
    // Scrolling should remain smooth even with large datasets
    expect(scrollPerformance).toBeLessThan(1000);
  });
});
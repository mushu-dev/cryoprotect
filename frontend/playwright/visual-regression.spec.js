// @ts-check
const { test, expect } = require('@playwright/test');

/**
 * Visual regression tests for the experimental data enhancement UI
 * 
 * These tests capture screenshots of key UI elements and pages for comparison
 * between builds to detect unintended visual changes.
 */
test.describe('Visual Regression Tests', () => {
  // Desktop viewport for consistent screenshots
  test.use({ viewport: { width: 1280, height: 720 } });

  test('experiments list page visual', async ({ page }) => {
    await page.goto('/experiments');
    
    // Wait for content to fully load
    await page.waitForSelector('[data-testid="experiment-card"]', { timeout: 5000 });
    
    // Take full page screenshot
    await page.screenshot({ 
      path: './test-results/visual/experiments-list.png',
      fullPage: true 
    });
    
    // Take screenshot of just the experiment cards container
    const cardsContainer = page.locator('[data-testid="experiments-list"]');
    if (await cardsContainer.count() > 0) {
      await cardsContainer.screenshot({ 
        path: './test-results/visual/experiment-cards.png' 
      });
    }
  });

  test('experiment details page visual', async ({ page }) => {
    await page.goto('/experiments');
    
    // Wait for content to load
    await page.waitForSelector('[data-testid="experiment-card"]', { timeout: 5000 });
    
    // Navigate to first experiment
    const firstCard = page.locator('[data-testid="experiment-card"]').first();
    if (await firstCard.count() > 0) {
      await firstCard.click();
      
      // Wait for navigation and content to load
      await page.waitForURL(/\/experiments\/\d+/);
      await page.waitForLoadState('networkidle');
      
      // Take full page screenshot
      await page.screenshot({ 
        path: './test-results/visual/experiment-details.png',
        fullPage: true 
      });
      
      // Take screenshot of experiment data visualization
      const dataViz = page.locator('[data-testid="experiment-chart"]');
      if (await dataViz.count() > 0) {
        await dataViz.screenshot({ 
          path: './test-results/visual/experiment-visualization.png' 
        });
      }
    }
  });

  test('experiment creation form visual', async ({ page }) => {
    await page.goto('/experiments/create');
    
    // Wait for form to load
    await page.waitForSelector('form', { timeout: 5000 });
    
    // Take screenshot of the form
    await page.locator('form').screenshot({ 
      path: './test-results/visual/experiment-form.png' 
    });
    
    // Fill in some test data
    await page.fill('input[name="name"]', 'Test Experiment');
    await page.fill('textarea[name="description"]', 'This is a test experiment for visual regression testing');
    
    // Take screenshot of filled form
    await page.locator('form').screenshot({ 
      path: './test-results/visual/experiment-form-filled.png' 
    });
  });

  test('responsive design breakpoints', async ({ page }) => {
    // Test mobile viewport
    await page.setViewportSize({ width: 375, height: 667 }); // iPhone SE
    await page.goto('/experiments');
    await page.waitForLoadState('networkidle');
    await page.screenshot({ 
      path: './test-results/visual/experiments-mobile.png',
      fullPage: true 
    });
    
    // Test tablet viewport
    await page.setViewportSize({ width: 768, height: 1024 }); // iPad
    await page.goto('/experiments');
    await page.waitForLoadState('networkidle');
    await page.screenshot({ 
      path: './test-results/visual/experiments-tablet.png',
      fullPage: true 
    });
    
    // Test desktop viewport
    await page.setViewportSize({ width: 1440, height: 900 }); // Desktop
    await page.goto('/experiments');
    await page.waitForLoadState('networkidle');
    await page.screenshot({ 
      path: './test-results/visual/experiments-desktop.png',
      fullPage: true 
    });
  });

  test('dark mode visual', async ({ page }) => {
    // Test light mode (default)
    await page.goto('/experiments');
    await page.waitForLoadState('networkidle');
    await page.screenshot({ 
      path: './test-results/visual/experiments-light.png',
      fullPage: true 
    });
    
    // Test dark mode if theme toggle exists
    const themeToggle = page.locator('[data-testid="theme-toggle"]');
    if (await themeToggle.count() > 0) {
      await themeToggle.click();
      await page.waitForTimeout(500); // Wait for theme transition
      await page.screenshot({ 
        path: './test-results/visual/experiments-dark.png',
        fullPage: true 
      });
    }
  });
});
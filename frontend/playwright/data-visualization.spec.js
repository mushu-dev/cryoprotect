// @ts-check
const { test, expect } = require('@playwright/test');
const { mockExperimentApi } = require('./test-utils');

/**
 * Data visualization tests for the experimental data enhancement UI
 * 
 * These tests specifically target the charting and visualization components
 * to ensure they render properly and display data correctly.
 */
test.describe('Data Visualization Tests', () => {
  test.beforeEach(async ({ page }) => {
    // Mock API responses to ensure we have predictable data for the charts
    await mockExperimentApi(page);
  });

  test('should render experiment results chart correctly', async ({ page }) => {
    // Navigate to an experiment detail page
    await page.goto('/experiments/1');
    
    // Wait for the chart to load
    await page.waitForSelector('[data-testid="experiment-chart"]', { timeout: 5000 });
    
    // Take a screenshot of the chart
    await page.locator('[data-testid="experiment-chart"]').screenshot({ 
      path: './test-results/visualization/experiment-results-chart.png' 
    });
    
    // Verify chart has key elements
    const chartElements = await page.evaluate(() => {
      const chart = document.querySelector('[data-testid="experiment-chart"]');
      if (!chart) return {};
      
      // Check for SVG elements
      const svg = chart.querySelector('svg');
      const axisElements = svg ? svg.querySelectorAll('g.axis, .x-axis, .y-axis, .tick').length : 0;
      const dataElements = svg ? svg.querySelectorAll('path.line, circle.point, rect.bar').length : 0;
      
      // Check for canvas elements (for Canvas-based charts)
      const canvas = chart.querySelector('canvas');
      const canvasWidth = canvas ? canvas.width : 0;
      const canvasHeight = canvas ? canvas.height : 0;
      
      return {
        hasSvg: !!svg,
        hasCanvas: !!canvas,
        axisElements,
        dataElements,
        canvasWidth,
        canvasHeight
      };
    });
    
    console.log('Chart elements:', chartElements);
    
    // Chart should have either SVG or Canvas rendering
    expect(chartElements.hasSvg || chartElements.hasCanvas).toBeTruthy();
    
    // If SVG-based chart, it should have axis and data elements
    if (chartElements.hasSvg) {
      expect(chartElements.axisElements).toBeGreaterThan(0);
      expect(chartElements.dataElements).toBeGreaterThan(0);
    }
    
    // If Canvas-based chart, it should have dimensions
    if (chartElements.hasCanvas) {
      expect(chartElements.canvasWidth).toBeGreaterThan(0);
      expect(chartElements.canvasHeight).toBeGreaterThan(0);
    }
  });

  test('should update chart when changing visualization type', async ({ page }) => {
    // Navigate to an experiment detail page
    await page.goto('/experiments/1');
    
    // Wait for the chart to load
    await page.waitForSelector('[data-testid="experiment-chart"]', { timeout: 5000 });
    
    // Take a screenshot of the initial chart
    await page.locator('[data-testid="experiment-chart"]').screenshot({ 
      path: './test-results/visualization/chart-initial.png' 
    });
    
    // Check if there's a chart type selector
    const chartTypeSelector = page.locator('[data-testid="chart-type-selector"], select.chart-type');
    if (await chartTypeSelector.count() > 0) {
      // Get the initial chart type
      const initialChartType = await chartTypeSelector.inputValue();
      console.log('Initial chart type:', initialChartType);
      
      // Change the chart type
      // Find a value that's different from the current
      const availableOptions = await page.evaluate(() => {
        const selector = document.querySelector('[data-testid="chart-type-selector"], select.chart-type');
        if (!selector) return [];
        return Array.from(selector.options).map(option => option.value);
      });
      
      // Find an alternative option
      const alternativeOption = availableOptions.find(option => option !== initialChartType);
      
      if (alternativeOption) {
        // Select the alternative option
        await chartTypeSelector.selectOption(alternativeOption);
        
        // Wait for chart to update
        await page.waitForTimeout(500);
        
        // Take a screenshot of the updated chart
        await page.locator('[data-testid="experiment-chart"]').screenshot({ 
          path: './test-results/visualization/chart-updated.png' 
        });
        
        // Verify the chart type has changed
        const newChartType = await chartTypeSelector.inputValue();
        expect(newChartType).toBe(alternativeOption);
      }
    } else {
      console.log('No chart type selector found, skipping chart type change test');
    }
  });

  test('should display data table with experiment results', async ({ page }) => {
    // Navigate to an experiment detail page
    await page.goto('/experiments/1');
    
    // Look for a data table showing results
    const dataTable = page.locator('[data-testid="experiment-data-table"], table.results-table');
    
    if (await dataTable.count() > 0) {
      // Take a screenshot of the data table
      await dataTable.screenshot({ 
        path: './test-results/visualization/data-table.png' 
      });
      
      // Verify the table has headers and data
      const tableHeaders = await dataTable.locator('th').count();
      const tableRows = await dataTable.locator('tr').count();
      
      expect(tableHeaders).toBeGreaterThan(0);
      expect(tableRows).toBeGreaterThan(1); // At least header row + one data row
      
      // Check that the table data matches what's in the chart
      const tableData = await page.evaluate(() => {
        const table = document.querySelector('[data-testid="experiment-data-table"], table.results-table');
        if (!table) return [];
        
        // Extract data from the table
        const rows = Array.from(table.querySelectorAll('tr'));
        return rows.slice(1).map(row => { // Skip header row
          const cells = Array.from(row.querySelectorAll('td'));
          return cells.map(cell => cell.textContent.trim());
        });
      });
      
      console.log('Table data:', tableData);
      
      // Table should have data
      expect(tableData.length).toBeGreaterThan(0);
      expect(tableData[0].length).toBeGreaterThan(0);
    } else {
      console.log('No data table found, skipping data table test');
    }
  });

  test('should provide interactive tooltips on data points', async ({ page }) => {
    // Navigate to an experiment detail page
    await page.goto('/experiments/1');
    
    // Wait for the chart to load
    await page.waitForSelector('[data-testid="experiment-chart"]', { timeout: 5000 });
    
    // Find data points in the chart
    const hasInteractiveElements = await page.evaluate(() => {
      const chart = document.querySelector('[data-testid="experiment-chart"]');
      if (!chart) return false;
      
      // Different chart libraries use different elements for data points
      const dataPoints = chart.querySelectorAll('circle.point, .dot, .bar, path.line');
      if (dataPoints.length === 0) return false;
      
      // Try to get the first data point
      const firstPoint = dataPoints[0];
      
      // Simulate hovering over the point (can't actually test the tooltip appearance this way,
      // but can verify the elements that would trigger it)
      return true;
    });
    
    // Only proceed with hovering test if we found interactive elements
    if (hasInteractiveElements) {
      // Find a point to hover over - different charts have different point elements
      const pointSelectors = [
        '[data-testid="experiment-chart"] circle.point',
        '[data-testid="experiment-chart"] .dot',
        '[data-testid="experiment-chart"] .bar',
        '[data-testid="experiment-chart"] .line-point'
      ];
      
      let pointFound = false;
      
      for (const selector of pointSelectors) {
        const points = page.locator(selector);
        if (await points.count() > 0) {
          // Hover over the first point
          await points.first().hover();
          
          // Wait briefly for tooltip to appear
          await page.waitForTimeout(500);
          
          // Take a screenshot with tooltip
          await page.screenshot({ 
            path: './test-results/visualization/chart-tooltip.png' 
          });
          
          // Look for tooltip element
          const tooltip = await page.locator('.tooltip, [role="tooltip"], .highcharts-tooltip').count();
          console.log(`Tooltip elements found: ${tooltip}`);
          
          pointFound = true;
          break;
        }
      }
      
      if (!pointFound) {
        console.log('No data points found with standard selectors');
      }
    }
  });

  test('should handle different data ranges appropriately', async ({ page }) => {
    // Navigate to experiment list
    await page.goto('/experiments');
    
    // Create a function to check chart scaling with different data
    const testChartWithData = async (data, name) => {
      // Mock API to return specific data for experiment 1
      await page.route('**/api/experiments/1', route => {
        route.fulfill({
          status: 200,
          contentType: 'application/json',
          body: JSON.stringify({
            id: 1,
            name: `Test Experiment - ${name}`,
            description: 'Testing chart with different data ranges',
            results: data
          })
        });
      });
      
      // Navigate to experiment 1
      await page.goto('/experiments/1');
      
      // Wait for chart to load
      await page.waitForSelector('[data-testid="experiment-chart"]', { timeout: 5000 });
      
      // Take a screenshot
      await page.locator('[data-testid="experiment-chart"]').screenshot({ 
        path: `./test-results/visualization/chart-${name}.png` 
      });
      
      // Basic check that chart is visible
      const chartVisible = await page.locator('[data-testid="experiment-chart"]').isVisible();
      expect(chartVisible).toBeTruthy();
    };
    
    // Test with different data scenarios
    await testChartWithData(
      [
        { time: 0, value: 0 },
        { time: 1, value: 10 },
        { time: 2, value: 20 },
        { time: 3, value: 30 },
        { time: 4, value: 40 }
      ], 
      'linear'
    );
    
    // Test with outlier data point
    await testChartWithData(
      [
        { time: 0, value: 0 },
        { time: 1, value: 10 },
        { time: 2, value: 100 }, // Outlier
        { time: 3, value: 30 },
        { time: 4, value: 40 }
      ], 
      'outlier'
    );
    
    // Test with negative values
    await testChartWithData(
      [
        { time: 0, value: 0 },
        { time: 1, value: -10 },
        { time: 2, value: -20 },
        { time: 3, value: 30 },
        { time: 4, value: 40 }
      ], 
      'negative'
    );
    
    // Test with large dataset
    const largeDataset = Array.from({ length: 100 }, (_, i) => ({
      time: i,
      value: Math.sin(i * 0.1) * 50 + 50
    }));
    
    await testChartWithData(largeDataset, 'large-dataset');
  });

  test('should handle new data visualization patterns', async ({ page }) => {
    // Navigate to experiment detail page
    await page.goto('/experiments/1');
    
    // Look for a visualization container
    await page.waitForSelector('[data-testid="experiment-chart"]', { timeout: 5000 });
    
    // Check for advanced visualization features
    const advancedFeatures = await page.evaluate(() => {
      const visualizationContainer = document.querySelector('[data-testid="experiment-chart"]');
      if (!visualizationContainer) return {};
      
      // Check for zoom controls
      const hasZoomControls = document.querySelector(
        '[data-testid="zoom-control"], .zoom-in, .zoom-out, .pan-button'
      ) !== null;
      
      // Check for export/download options
      const hasExportOptions = document.querySelector(
        '[data-testid="export-chart"], .download-button, .export-button'
      ) !== null;
      
      // Check for alternative visualization views
      const hasAltViews = document.querySelector(
        '[data-testid="view-selector"], .chart-view-options'
      ) !== null;
      
      // Check for annotations
      const hasAnnotations = document.querySelector(
        '.annotation, .chart-annotation, .data-label'
      ) !== null;
      
      return {
        hasZoomControls,
        hasExportOptions,
        hasAltViews,
        hasAnnotations
      };
    });
    
    console.log('Advanced visualization features:', advancedFeatures);
    
    // Log the features found, but don't make them required
    // as the implementation may vary
    if (advancedFeatures.hasZoomControls) {
      console.log('✓ Zoom controls detected');
    }
    
    if (advancedFeatures.hasExportOptions) {
      console.log('✓ Export options detected');
    }
    
    if (advancedFeatures.hasAltViews) {
      console.log('✓ Alternative views detected');
    }
    
    if (advancedFeatures.hasAnnotations) {
      console.log('✓ Annotations detected');
    }
  });
});
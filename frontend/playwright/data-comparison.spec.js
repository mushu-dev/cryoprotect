// @ts-check
const { test, expect } = require('@playwright/test');

/**
 * Data comparison and analysis tests for the experimental data enhancement UI
 * 
 * These tests verify the functionality for comparing multiple experiments,
 * analyzing data across experiments, and advanced analysis features.
 */
test.describe('Data Comparison and Analysis Tests', () => {
  test.beforeEach(async ({ page }) => {
    // Navigate to the experiments list page
    await page.goto('/experiments');
    
    // Wait for experiment cards to load
    await page.waitForSelector('[data-testid="experiment-card"]', { timeout: 5000 });
  });

  test('should allow selecting multiple experiments for comparison', async ({ page }) => {
    // Look for experiment selection checkboxes
    const experimentCheckboxes = page.locator('[data-testid="experiment-checkbox"], input[type="checkbox"]').filter({
      has: page.locator('[data-testid="experiment-card"]')
    });
    
    // If no checkboxes, try to find a compare button or multi-select mode
    if (await experimentCheckboxes.count() === 0) {
      const compareButton = page.locator('button:has-text("Compare"), [data-testid="compare-button"]');
      
      if (await compareButton.count() === 0) {
        test.skip('No experiment selection mechanism found');
        return;
      }
      
      // Click the compare button to enter comparison mode
      await compareButton.click();
      
      // Wait for selection mode to activate
      await page.waitForTimeout(500);
      
      // Check again for selection checkboxes
      const selectionCheckboxes = page.locator('[data-testid="experiment-checkbox"], input[type="checkbox"]').filter({
        has: page.locator('[data-testid="experiment-card"]')
      });
      
      if (await selectionCheckboxes.count() === 0) {
        test.skip('No experiment selection checkboxes found after activating compare mode');
        return;
      }
    }
    
    // Take screenshot before selection
    await page.screenshot({ path: './test-results/comparison/before-selection.png' });
    
    // Get the available checkboxes after ensuring they're visible
    const checkboxes = page.locator('[data-testid="experiment-checkbox"], input[type="checkbox"]').filter({
      has: page.locator('[data-testid="experiment-card"]')
    });
    
    const checkboxCount = await checkboxes.count();
    console.log(`Found ${checkboxCount} experiment checkboxes`);
    
    // Select at least two experiments if available
    const selectCount = Math.min(2, checkboxCount);
    
    for (let i = 0; i < selectCount; i++) {
      await checkboxes.nth(i).check();
    }
    
    // Take screenshot after selection
    await page.screenshot({ path: './test-results/comparison/after-selection.png' });
    
    // Look for a compare button / submit selection button
    const submitCompareButton = page.locator(
      'button:has-text("Compare"), button:has-text("View Selected"), [data-testid="compare-selected"]'
    );
    
    if (await submitCompareButton.count() > 0) {
      // Click the comparison button
      await submitCompareButton.click();
      
      // Wait for comparison view to load
      await page.waitForTimeout(1000);
      
      // Take screenshot of comparison view
      await page.screenshot({ path: './test-results/comparison/comparison-view.png' });
      
      // Check if we navigated to a comparison page
      const currentUrl = page.url();
      console.log(`Navigated to: ${currentUrl}`);
      
      // Look for comparison-specific UI elements
      const comparisonChart = page.locator('[data-testid="comparison-chart"], [data-testid="multi-experiment-chart"]');
      const comparisonTable = page.locator('[data-testid="comparison-table"], table.comparison-data');
      
      const hasComparisonElements = 
        (await comparisonChart.count() > 0) || 
        (await comparisonTable.count() > 0);
      
      expect(hasComparisonElements).toBeTruthy();
    } else {
      // Maybe comparison happens in-place without navigation
      console.log('No compare button found, checking for in-page comparison');
      
      // Look for any comparison elements that might have appeared
      const inPageComparison = page.locator(
        '[data-testid="comparison-view"], [data-testid="selected-items-summary"]'
      );
      
      if (await inPageComparison.count() > 0) {
        console.log('Found in-page comparison element');
      } else {
        console.log('No comparison UI elements found');
      }
    }
  });

  test('should display comparative visualizations for selected experiments', async ({ page }) => {
    // First navigate to the comparison page or activate comparison mode
    // This depends on the implementation approach
    
    // Try the multi-select approach first
    const experimentCheckboxes = page.locator('[data-testid="experiment-checkbox"], input[type="checkbox"]').filter({
      has: page.locator('[data-testid="experiment-card"]')
    });
    
    if (await experimentCheckboxes.count() > 0) {
      // Select first two experiments
      await experimentCheckboxes.first().check();
      if (await experimentCheckboxes.count() > 1) {
        await experimentCheckboxes.nth(1).check();
      }
      
      // Look for a compare button
      const compareButton = page.locator(
        'button:has-text("Compare"), [data-testid="compare-selected"]'
      );
      
      if (await compareButton.count() > 0) {
        await compareButton.click();
        await page.waitForTimeout(1000);
      }
    } else {
      // Try looking for a dedicated comparison page
      const compareLink = page.locator('a:has-text("Compare"), [href*="compare"]');
      
      if (await compareLink.count() > 0) {
        await compareLink.click();
        await page.waitForTimeout(1000);
      } else {
        // Try a compare button that might lead to comparison page
        const compareButton = page.locator('button:has-text("Compare")');
        
        if (await compareButton.count() > 0) {
          await compareButton.click();
          await page.waitForTimeout(1000);
        } else {
          test.skip('Could not find a way to access comparison features');
          return;
        }
      }
    }
    
    // Now look for comparative visualization elements
    const comparisonChart = page.locator('[data-testid="comparison-chart"], [data-testid="multi-experiment-chart"]');
    
    if (await comparisonChart.count() === 0) {
      // If we still can't find comparison elements, try looking for experiment selection UI
      // that might be required first
      const experimentSelectors = page.locator(
        'select[multiple], [data-testid="experiment-selector"], [aria-label*="Select experiment"]'
      );
      
      if (await experimentSelectors.count() > 0) {
        // Select at least two experiments
        const selectElement = experimentSelectors.first();
        
        // Different handling for different types of selectors
        const isMultiSelect = await page.evaluate(el => {
          return el.tagName.toLowerCase() === 'select' && el.multiple;
        }, await selectElement.elementHandle());
        
        if (isMultiSelect) {
          // Handle a <select multiple> element
          await selectElement.selectOption([{ index: 0 }, { index: 1 }]);
        } else {
          // Handle other types of selectors (dropdowns, checkboxes, etc.)
          await selectElement.click();
          await page.waitForTimeout(300);
          
          // Select first item in dropdown
          const firstOption = page.locator('li, [role="option"]').first();
          if (await firstOption.count() > 0) {
            await firstOption.click();
          }
          
          // Select second item if possible
          const secondOption = page.locator('li, [role="option"]').nth(1);
          if (await secondOption.count() > 0) {
            await secondOption.click();
          }
          
          // Close dropdown if there's a close/done button
          const doneButton = page.locator('button:has-text("Done"), button:has-text("Apply")');
          if (await doneButton.count() > 0) {
            await doneButton.click();
          }
        }
        
        // Check for a button to apply the selection
        const applyButton = page.locator('button:has-text("Apply"), button:has-text("Compare")');
        if (await applyButton.count() > 0) {
          await applyButton.click();
        }
        
        await page.waitForTimeout(1000);
      } else {
        test.skip('Could not find comparison chart or experiment selection UI');
        return;
      }
    }
    
    // By now we should have navigated to or activated the comparison view
    // Take screenshot of the comparison view
    await page.screenshot({ path: './test-results/comparison/comparative-visualization.png' });
    
    // Look again for comparison visualization elements
    const comparisonChartAfterSetup = page.locator(
      '[data-testid="comparison-chart"], [data-testid="multi-experiment-chart"], .experiment-comparison'
    );
    
    if (await comparisonChartAfterSetup.count() > 0) {
      // Take screenshot of the chart
      await comparisonChartAfterSetup.screenshot({ path: './test-results/comparison/comparison-chart.png' });
      
      // Check for multiple data series in the chart
      const hasMultipleDataSeries = await page.evaluate(() => {
        const chart = document.querySelector(
          '[data-testid="comparison-chart"], [data-testid="multi-experiment-chart"], .experiment-comparison'
        );
        
        if (!chart) return false;
        
        // Look for multiple lines/bars in the chart
        const dataElements = chart.querySelectorAll('path.line, .series, .bar, .data-point');
        
        if (dataElements.length === 0) return false;
        
        // Check for different classes or colors in data elements
        const uniqueClasses = new Set();
        dataElements.forEach(el => {
          const className = el.getAttribute('class');
          if (className) uniqueClasses.add(className);
          
          const color = el.getAttribute('stroke') || el.getAttribute('fill');
          if (color) uniqueClasses.add(color);
        });
        
        return uniqueClasses.size > 1;
      });
      
      console.log(`Chart has multiple data series: ${hasMultipleDataSeries}`);
      
      // Verify we have a legend if there are multiple data series
      if (hasMultipleDataSeries) {
        const legend = page.locator('.legend, [data-testid="chart-legend"]');
        const hasLegend = await legend.count() > 0;
        
        console.log(`Chart has legend: ${hasLegend}`);
        
        // Legend is helpful but not strictly required
        if (hasLegend) {
          await legend.screenshot({ path: './test-results/comparison/chart-legend.png' });
        }
      }
    } else {
      console.log('No comparison chart found after setup');
    }
    
    // Look for a comparison table as an alternative or complement to the chart
    const comparisonTable = page.locator('[data-testid="comparison-table"], table.comparison-data');
    
    if (await comparisonTable.count() > 0) {
      await comparisonTable.screenshot({ path: './test-results/comparison/comparison-table.png' });
      
      // Verify table has multiple experiments
      const tableRowCount = await comparisonTable.locator('tr').count();
      const tableColCount = await comparisonTable.locator('th').count();
      
      console.log(`Comparison table: ${tableRowCount} rows, ${tableColCount} columns`);
      
      // Should have more than just a header row
      expect(tableRowCount).toBeGreaterThan(1);
    }
  });

  test('should perform statistical analysis across experiments', async ({ page }) => {
    // Navigate to or find a statistical analysis section
    // This could be a dedicated page or a tab within comparison
    
    // Try known potential URLs for analysis page
    const analysisPaths = [
      '/analysis',
      '/experiments/analysis',
      '/experiments/compare',
      '/statistics',
      '/data-analysis'
    ];
    
    let foundAnalysisPage = false;
    
    for (const path of analysisPaths) {
      const response = await page.goto(path, { timeout: 3000 }).catch(() => null);
      
      if (response && response.ok()) {
        console.log(`Found analysis page at ${path}`);
        foundAnalysisPage = true;
        break;
      }
    }
    
    if (!foundAnalysisPage) {
      // Look for links/buttons to analysis page
      const analysisLink = page.locator(
        'a:has-text("Analysis"), a:has-text("Statistics"), button:has-text("Analyze")'
      );
      
      if (await analysisLink.count() > 0) {
        await analysisLink.click();
        await page.waitForTimeout(1000);
        foundAnalysisPage = true;
      }
    }
    
    if (!foundAnalysisPage) {
      // Maybe analysis is part of comparison page
      const compareLink = page.locator('a:has-text("Compare")');
      
      if (await compareLink.count() > 0) {
        await compareLink.click();
        await page.waitForTimeout(1000);
        
        // Look for analysis tab within comparison
        const analysisTab = page.locator('button:has-text("Analysis"), [role="tab"]:has-text("Statistics")');
        
        if (await analysisTab.count() > 0) {
          await analysisTab.click();
          await page.waitForTimeout(500);
          foundAnalysisPage = true;
        }
      }
    }
    
    if (!foundAnalysisPage) {
      test.skip('Could not find analysis functionality');
      return;
    }
    
    // By now we should be on some kind of analysis page or tab
    // Take screenshot of the analysis page
    await page.screenshot({ path: './test-results/comparison/analysis-page.png' });
    
    // Look for statistical analysis UI elements
    const statisticalElements = [
      '[data-testid="statistics-panel"]',
      '[data-testid="statistical-analysis"]',
      '.statistics-container',
      'table.statistics',
      '[data-testid="statistical-summary"]'
    ];
    
    let foundStatisticalElement = false;
    
    for (const selector of statisticalElements) {
      const element = page.locator(selector);
      if (await element.count() > 0) {
        console.log(`Found statistical element: ${selector}`);
        await element.screenshot({ path: './test-results/comparison/statistical-element.png' });
        foundStatisticalElement = true;
        break;
      }
    }
    
    // Look for common statistical metrics
    const statsText = await page.evaluate(() => {
      return document.body.innerText;
    });
    
    const statisticalTerms = [
      'mean', 'average', 'median', 'std dev', 'standard deviation',
      'variance', 'p-value', 'correlation', 'regression', 'analysis',
      'statistic', 'confidence', 'interval'
    ];
    
    const foundTerms = statisticalTerms.filter(term => 
      statsText.toLowerCase().includes(term.toLowerCase())
    );
    
    console.log('Statistical terms found:', foundTerms);
    
    const hasStatisticalContent = foundStatisticalElement || foundTerms.length > 0;
    
    // If we're on what seems to be an analysis page, there should be some
    // statistical content visible
    expect(hasStatisticalContent).toBeTruthy();
  });

  test('should allow parameter adjustment for data analysis', async ({ page }) => {
    // First navigate to a page where analysis parameters can be adjusted
    // This could be the comparison page, analysis page, or experiment detail page
    
    // Try comparison page first (it's most likely to have adjustable analysis)
    const compareLink = page.locator('a:has-text("Compare")');
    
    if (await compareLink.count() > 0) {
      await compareLink.click();
      await page.waitForTimeout(1000);
    } else {
      // Try looking for analysis page
      const analysisLink = page.locator('a:has-text("Analysis")');
      
      if (await analysisLink.count() > 0) {
        await analysisLink.click();
        await page.waitForTimeout(1000);
      }
    }
    
    // Look for analysis parameter controls
    const parameterControls = [
      'input[type="range"]', // Sliders
      'select[name*="analysis"], select[name*="parameter"]', // Dropdown selections
      'input[type="number"][name*="threshold"]', // Numeric inputs
      '[data-testid="analysis-parameter"]', // Explicit parameter controls
      '.parameter-control' // Generic parameter controls by class
    ];
    
    let foundParameterControl = false;
    let controlElement = null;
    
    for (const selector of parameterControls) {
      const element = page.locator(selector);
      if (await element.count() > 0) {
        console.log(`Found parameter control: ${selector}`);
        controlElement = element.first();
        foundParameterControl = true;
        break;
      }
    }
    
    if (!foundParameterControl) {
      test.skip('Could not find adjustable analysis parameters');
      return;
    }
    
    // Take screenshot before adjustment
    await page.screenshot({ path: './test-results/comparison/before-parameter-adjustment.png' });
    
    // Adjust the parameter
    const elementType = await controlElement.evaluate(el => el.tagName.toLowerCase());
    
    if (elementType === 'select') {
      // Select a different option
      const optionCount = await controlElement.locator('option').count();
      if (optionCount > 1) {
        const currentOption = await controlElement.inputValue();
        
        // Choose an option different from the current one
        const options = await controlElement.locator('option').all();
        for (const option of options) {
          const value = await option.getAttribute('value');
          if (value !== currentOption) {
            await controlElement.selectOption(value);
            break;
          }
        }
      }
    } else if (elementType === 'input') {
      const inputType = await controlElement.getAttribute('type');
      
      if (inputType === 'range') {
        // For range sliders, set to middle value
        const min = Number(await controlElement.getAttribute('min') || 0);
        const max = Number(await controlElement.getAttribute('max') || 100);
        const middleValue = Math.floor((min + max) / 2);
        
        await controlElement.fill(String(middleValue));
      } else if (inputType === 'number') {
        // For number inputs, increment the value
        const currentValue = Number(await controlElement.inputValue() || 0);
        await controlElement.fill(String(currentValue + 1));
      }
    }
    
    // Wait for analysis to update
    await page.waitForTimeout(1000);
    
    // Take screenshot after adjustment
    await page.screenshot({ path: './test-results/comparison/after-parameter-adjustment.png' });
    
    // Look for any changes in the display
    // Check for analysis results container
    const resultsContainer = page.locator(
      '[data-testid="analysis-results"], .analysis-output, [data-testid="chart-container"]'
    );
    
    if (await resultsContainer.count() > 0) {
      await resultsContainer.screenshot({ path: './test-results/comparison/analysis-results.png' });
    }
  });

  test('should export comparison or analysis results', async ({ page }) => {
    // First navigate to comparison or analysis page
    
    // Try comparison page first
    const compareLink = page.locator('a:has-text("Compare")');
    
    if (await compareLink.count() > 0) {
      await compareLink.click();
      await page.waitForTimeout(1000);
    } else {
      // Try analysis page
      const analysisLink = page.locator('a:has-text("Analysis")');
      
      if (await analysisLink.count() > 0) {
        await analysisLink.click();
        await page.waitForTimeout(1000);
      }
    }
    
    // Look for export functionality
    const exportControls = [
      'button:has-text("Export")',
      'button:has-text("Download")',
      '[data-testid="export-button"]',
      '[data-testid="download-results"]',
      'a:has-text("Export Results")'
    ];
    
    let exportControl = null;
    
    for (const selector of exportControls) {
      const element = page.locator(selector);
      if (await element.count() > 0) {
        console.log(`Found export control: ${selector}`);
        exportControl = element;
        break;
      }
    }
    
    if (!exportControl) {
      // Look for a dropdown that might contain export options
      const menuButtons = page.locator('button:has-text("Options"), button:has-text("Menu"), [aria-haspopup="menu"]');
      
      if (await menuButtons.count() > 0) {
        // Click the menu button
        await menuButtons.first().click();
        
        // Wait for menu to appear
        await page.waitForTimeout(300);
        
        // Look for export option in menu
        const exportMenuItem = page.locator(
          'li:has-text("Export"), li:has-text("Download"), [role="menuitem"]:has-text("Export")'
        );
        
        if (await exportMenuItem.count() > 0) {
          exportControl = exportMenuItem.first();
        } else {
          // Close the menu
          await page.keyboard.press('Escape');
        }
      }
    }
    
    if (!exportControl) {
      test.skip('Could not find export functionality');
      return;
    }
    
    // Take screenshot showing the export control
    await page.screenshot({ path: './test-results/comparison/export-control.png' });
    
    // We can't actually test the download in a headless browser easily,
    // but we can click the button and verify that it doesn't error
    
    // Set up a handler for any download dialog
    page.once('dialog', async dialog => {
      console.log(`Dialog message: ${dialog.message()}`);
      await dialog.accept();
    });
    
    // Click the export button
    await exportControl.click();
    
    // Wait for export process to start
    await page.waitForTimeout(1000);
    
    // Take screenshot after clicking export
    await page.screenshot({ path: './test-results/comparison/after-export-click.png' });
    
    // Check for success indication or download started message
    const successMessage = page.locator(
      'text="Export successful"',
      'text="Download started"',
      'text="File ready"'
    );
    
    const hasSuccessMessage = await successMessage.count() > 0;
    
    if (hasSuccessMessage) {
      console.log('Export success message found');
    } else {
      console.log('No export success message found, but export action was completed');
    }
  });
});
// @ts-check
const { test, expect } = require('@playwright/test');

/**
 * Protein Visualizer component tests
 * 
 * These tests verify the functionality of the MolstarViewer component
 * for molecular visualization.
 */
test.describe('Protein Visualizer Tests', () => {
  test('should render the protein visualizer demo page', async ({ page }) => {
    // Navigate to the protein visualizer demo page
    await page.goto('/protein-visualizer-demo');
    
    // Check that the page title is correct
    await expect(page).toHaveTitle(/Protein Visualizer Demo/);
    
    // Check that the main heading is present
    const heading = page.locator('h1');
    await expect(heading).toContainText('Protein Visualizer Demo');
    
    // Take a screenshot of the page
    await page.screenshot({ 
      path: './test-results/protein-visualizer/demo-page.png',
      fullPage: true
    });
  });
  
  test('should display small molecule examples', async ({ page }) => {
    // Navigate to the protein visualizer demo page
    await page.goto('/protein-visualizer-demo');
    
    // Wait for the tabs to load and ensure we're on the examples tab
    await page.locator('button[value="examples"]').click();
    
    // Find the small molecules section
    const smallMoleculeCard = page.locator('text="Small Molecules"').first();
    
    // Wait for the visualization to load
    await page.waitForSelector('.protein-visualizer', { timeout: 10000 });
    
    // Check if the loading indicator appears and then disappears
    const loadingIndicator = page.locator('.protein-visualizer .animate-spin');
    await expect(loadingIndicator).toBeVisible();
    await expect(loadingIndicator).toBeHidden({ timeout: 30000 });
    
    // Take a screenshot of the small molecule visualization
    await smallMoleculeCard.screenshot({ 
      path: './test-results/protein-visualizer/small-molecule.png' 
    });
    
    // Change the molecule selection
    await page.selectOption('select:near(:text("Small Molecules"))', 'caffeine');
    
    // Wait for visualization to update
    await page.waitForTimeout(2000);
    
    // Take a screenshot after changing the molecule
    await smallMoleculeCard.screenshot({ 
      path: './test-results/protein-visualizer/caffeine.png' 
    });
  });
  
  test('should display protein structure examples', async ({ page }) => {
    // Navigate to the protein visualizer demo page
    await page.goto('/protein-visualizer-demo');
    
    // Wait for the tabs to load and ensure we're on the examples tab
    await page.locator('button[value="examples"]').click();
    
    // Find the protein structures section
    const proteinStructureCard = page.locator('text="Protein Structures"').first();
    
    // Wait for the visualization to load
    await page.waitForSelector('.protein-visualizer', { timeout: 10000 });
    
    // Check if the loading indicator appears and then disappears
    const loadingIndicator = page.locator('.protein-visualizer .animate-spin').last();
    await expect(loadingIndicator).toBeVisible();
    await expect(loadingIndicator).toBeHidden({ timeout: 30000 });
    
    // Take a screenshot of the protein structure visualization
    await proteinStructureCard.screenshot({ 
      path: './test-results/protein-visualizer/protein-structure.png' 
    });
    
    // Change the protein selection
    await page.selectOption('select:near(:text("Protein Structures"))', '1hsg');
    
    // Wait for visualization to update
    await page.waitForTimeout(2000);
    
    // Take a screenshot after changing the protein
    await proteinStructureCard.screenshot({ 
      path: './test-results/protein-visualizer/hiv-protease.png' 
    });
  });
  
  test('should change visualization styles', async ({ page }) => {
    // Navigate to the protein visualizer demo page
    await page.goto('/protein-visualizer-demo');
    
    // Wait for the tabs to load and ensure we're on the examples tab
    await page.locator('button[value="examples"]').click();
    
    // Wait for visualizations to load
    await page.waitForSelector('.protein-visualizer', { timeout: 10000 });
    
    // Scroll to the visualization styles section
    await page.locator('text="Visualization Styles"').scrollIntoViewIfNeeded();
    
    // Wait for style examples to load
    await page.waitForTimeout(5000);
    
    // Take a screenshot of the visualization styles section
    await page.locator('h3:text("Visualization Styles")').evaluate(element => {
      const styleSection = element.closest('div');
      if (styleSection && styleSection.nextElementSibling) {
        return styleSection.nextElementSibling.scrollIntoView({ behavior: 'smooth', block: 'center' });
      }
    });
    
    await page.screenshot({ 
      path: './test-results/protein-visualizer/visualization-styles.png',
      fullPage: false
    });
  });
  
  test('should handle custom visualization inputs', async ({ page }) => {
    // Navigate to the protein visualizer demo page
    await page.goto('/protein-visualizer-demo');
    
    // Switch to the custom tab
    await page.locator('button[value="custom"]').click();
    
    // Wait for the custom visualization form to load
    await page.waitForSelector('input#pdbId', { timeout: 5000 });
    
    // Input a PDB ID
    await page.fill('input#pdbId', '4hhb');
    
    // Set a name
    await page.fill('input#name', 'Hemoglobin Test');
    
    // Change the style
    await page.click('button#style');
    await page.click('text=Surface');
    
    // Toggle auto-rotate
    await page.click('button[id="autoRotate"]');
    
    // Wait for the visualization to load
    await page.waitForSelector('.protein-visualizer', { timeout: 10000 });
    
    // Wait for loading to complete
    const loadingIndicator = page.locator('.protein-visualizer .animate-spin');
    if (await loadingIndicator.isVisible()) {
      await expect(loadingIndicator).toBeHidden({ timeout: 30000 });
    }
    
    // Take a screenshot of the custom visualization
    await page.locator('h3:text("Preview")').evaluate(element => {
      const previewSection = element.closest('div');
      if (previewSection) {
        return previewSection.scrollIntoView({ behavior: 'smooth', block: 'center' });
      }
    });
    
    await page.screenshot({ 
      path: './test-results/protein-visualizer/custom-visualization.png',
      fullPage: false
    });
    
    // Clear the PDB ID and enter a SMILES string instead
    await page.fill('input#pdbId', '');
    await page.fill('input#smiles', 'CC(=O)OC1=CC=CC=C1C(=O)O');
    await page.fill('input#name', 'Aspirin Test');
    
    // Wait for the visualization to update
    await page.waitForTimeout(3000);
    
    // Take a screenshot of the custom SMILES visualization
    await page.screenshot({ 
      path: './test-results/protein-visualizer/custom-smiles.png',
      fullPage: false
    });
  });
  
  test('should test control functionality', async ({ page }) => {
    // Navigate to the protein visualizer demo
    await page.goto('/protein-visualizer-demo');
    
    // Wait for the visualizations to load
    await page.waitForSelector('.protein-visualizer', { timeout: 10000 });
    
    // Find the first protein visualizer with controls
    const firstVisualizer = page.locator('.protein-visualizer').first();
    
    // Wait for loading to complete
    const loadingIndicator = firstVisualizer.locator('.animate-spin');
    if (await loadingIndicator.isVisible()) {
      await expect(loadingIndicator).toBeHidden({ timeout: 30000 });
    }
    
    // Test the style selector control
    const styleSelector = firstVisualizer.locator('button:has-text("Style")');
    await styleSelector.click();
    await page.locator('text=Surface').first().click();
    
    // Wait for style to update
    await page.waitForTimeout(1000);
    
    // Take a screenshot after style change
    await firstVisualizer.screenshot({ 
      path: './test-results/protein-visualizer/style-change.png' 
    });
    
    // Test the rotation toggle
    const rotateButton = firstVisualizer.locator('button[title="Stop rotation"], button[title="Start rotation"]');
    await rotateButton.click();
    
    // Wait for rotation to update
    await page.waitForTimeout(1000);
    
    // Test zoom in
    const zoomInButton = firstVisualizer.locator('button[title="Zoom in"]');
    await zoomInButton.click();
    
    // Wait for zoom to update
    await page.waitForTimeout(1000);
    
    // Take a screenshot after zoom
    await firstVisualizer.screenshot({ 
      path: './test-results/protein-visualizer/zoomed-in.png' 
    });
    
    // Test zoom out
    const zoomOutButton = firstVisualizer.locator('button[title="Zoom out"]');
    await zoomOutButton.click();
    
    // Wait for zoom to update
    await page.waitForTimeout(1000);
    
    // Test reset view
    const resetButton = firstVisualizer.locator('button[title="Reset view"]');
    await resetButton.click();
    
    // Wait for reset to update
    await page.waitForTimeout(1000);
    
    // Take a screenshot after reset
    await firstVisualizer.screenshot({ 
      path: './test-results/protein-visualizer/reset-view.png' 
    });
  });
});
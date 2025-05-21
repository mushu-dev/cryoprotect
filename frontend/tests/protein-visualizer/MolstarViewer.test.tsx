import { test, expect } from '@playwright/test';

// Tests for the MolstarViewer component
test.describe('MolstarViewer Component', () => {
  test('renders with SMILES string input', async ({ page }) => {
    // Go to the page that uses MolstarViewer with SMILES input
    await page.goto('/test-protein-visualizer?type=smiles');
    
    // Check that the viewer container is visible
    await expect(page.locator('.protein-visualizer')).toBeVisible();
    
    // Check that molecule name is displayed
    await expect(page.locator('text=Aspirin')).toBeVisible();
    
    // Check that controls are visible
    await expect(page.locator('button[title="Reset view"]')).toBeVisible();
    await expect(page.locator('button[title="Zoom in"]')).toBeVisible();
    await expect(page.locator('button[title="Zoom out"]')).toBeVisible();
  });
  
  test('renders with PDB ID input', async ({ page }) => {
    // Go to the page that uses MolstarViewer with PDB ID input
    await page.goto('/test-protein-visualizer?type=pdb');
    
    // Check that the viewer container is visible
    await expect(page.locator('.protein-visualizer')).toBeVisible();
    
    // Check that molecule name is displayed
    await expect(page.locator('text=Protein Structure')).toBeVisible();
    
    // Check that PDB ID is shown somewhere in the viewer
    await expect(page.locator('text=1cbs')).toBeVisible();
  });
  
  test('changes visualization style', async ({ page }) => {
    // Go to the test page
    await page.goto('/test-protein-visualizer?type=smiles');
    
    // Open the style dropdown
    await page.click('button:has-text("Cartoon")');
    
    // Select the Ball & Stick style
    await page.click('text=Ball & Stick');
    
    // Verify that the style was changed
    await expect(page.locator('text=Style: ball-and-stick')).toBeVisible();
  });
  
  test('toggles rotation', async ({ page }) => {
    // Go to the test page
    await page.goto('/test-protein-visualizer?type=smiles');
    
    // Get the rotation button
    const rotateButton = page.locator('button[title="Start rotation"]');
    
    // Click to start rotation
    await rotateButton.click();
    
    // Check that rotation button is now active (has primary color)
    await expect(rotateButton).toHaveClass(/text-primary/);
    
    // Click to stop rotation
    await rotateButton.click();
    
    // Check that rotation button is now inactive
    await expect(rotateButton).not.toHaveClass(/text-primary/);
  });
  
  test('resets view when reset button is clicked', async ({ page }) => {
    // Go to the test page
    await page.goto('/test-protein-visualizer?type=smiles');
    
    // Find and click the reset view button
    await page.click('button[title="Reset view"]');
    
    // This is a visual test, so we'll take a screenshot to verify later
    await page.screenshot({ path: 'test-results/protein-visualizer-reset.png' });
  });
  
  test('handles errors gracefully', async ({ page }) => {
    // Go to the test page with an invalid molecule
    await page.goto('/test-protein-visualizer?type=invalid');
    
    // Check that error message is shown
    await expect(page.locator('text=Error')).toBeVisible();
  });
  
  test('shows loading state', async ({ page }) => {
    // Go to the test page with a slow-loading molecule
    await page.goto('/test-protein-visualizer?type=slow');
    
    // Check that loading indicator is shown
    await expect(page.locator('text=Loading molecule')).toBeVisible();
    
    // Wait for loading to complete
    await page.waitForSelector('text=Loading molecule', { state: 'hidden' });
    
    // Check that the viewer is visible
    await expect(page.locator('.protein-visualizer')).toBeVisible();
  });
  
  test('supports highlighting specific residues', async ({ page }) => {
    // Go to the test page with residue highlighting
    await page.goto('/test-protein-visualizer?type=highlight');
    
    // Take a screenshot to verify highlighting
    await page.screenshot({ path: 'test-results/protein-visualizer-highlighting.png' });
  });
  
  test('handles download button click', async ({ page }) => {
    // Mock the download functionality since we can't easily test file downloads
    await page.goto('/test-protein-visualizer?type=smiles');
    
    // Set up download handler
    const downloadPromise = page.waitForEvent('download');
    
    // Click download button
    await page.click('button[title="Download image"]');
    
    // Wait for download to start
    const download = await downloadPromise;
    
    // Verify filename contains molecule name
    expect(download.suggestedFilename()).toContain('Aspirin');
  });
  
  test('supports custom height and width', async ({ page }) => {
    // Go to the test page with custom dimensions
    await page.goto('/test-protein-visualizer?type=custom-size');
    
    // Get the container element
    const container = page.locator('.protein-visualizer');
    
    // Check dimensions
    await expect(container).toHaveAttribute('style', /height: 600px/);
  });
  
  test('supports custom background color', async ({ page }) => {
    // Go to the test page with custom background
    await page.goto('/test-protein-visualizer?type=custom-bg');
    
    // Take a screenshot to verify background color
    await page.screenshot({ path: 'test-results/protein-visualizer-custom-bg.png' });
  });
  
  test('supports hiding controls', async ({ page }) => {
    // Go to the test page with hidden controls
    await page.goto('/test-protein-visualizer?type=no-controls');
    
    // Verify that control buttons are not visible
    await expect(page.locator('button[title="Reset view"]')).not.toBeVisible();
  });
});

// This is a helper test to set up a test page for the MolstarViewer tests
test.beforeAll(async () => {
  console.log('Setting up test pages for MolstarViewer tests...');
  
  // Note: In a real setup, you would create actual test pages
  // This is a placeholder for documentation purposes
});

// Test page setup (create pages at /pages/test-protein-visualizer.tsx):
/*
import { useRouter } from 'next/router';
import { MolstarViewer } from '@/components/protein-visualizer/MolstarViewer';

export default function TestProteinViewer() {
  const router = useRouter();
  const { type } = router.query;

  // Different test cases based on query parameter
  switch (type) {
    case 'pdb':
      return (
        <div className="container mx-auto p-4">
          <h1 className="text-2xl font-bold mb-4">Protein Structure Test</h1>
          <MolstarViewer 
            pdbId="1cbs" 
            name="Protein Structure" 
            height={400} 
            style="cartoon" 
          />
        </div>
      );
    
    case 'invalid':
      return (
        <div className="container mx-auto p-4">
          <h1 className="text-2xl font-bold mb-4">Invalid Structure Test</h1>
          <MolstarViewer 
            smiles="!!!INVALID!!!" 
            name="Invalid Molecule" 
            height={400} 
          />
        </div>
      );
    
    case 'slow':
      return (
        <div className="container mx-auto p-4">
          <h1 className="text-2xl font-bold mb-4">Slow Loading Test</h1>
          <MolstarViewer 
            url="/api/test/slow-molecule?delay=2000" 
            name="Slow Molecule" 
            height={400} 
          />
        </div>
      );
    
    case 'highlight':
      return (
        <div className="container mx-auto p-4">
          <h1 className="text-2xl font-bold mb-4">Residue Highlighting Test</h1>
          <MolstarViewer 
            pdbId="1hsg" 
            name="HIV-1 Protease" 
            style="cartoon" 
            highlightResidues={["A:25", "A:26", "A:27"]} 
            height={400} 
          />
        </div>
      );
    
    case 'custom-size':
      return (
        <div className="container mx-auto p-4">
          <h1 className="text-2xl font-bold mb-4">Custom Size Test</h1>
          <MolstarViewer 
            smiles="CC(=O)OC1=CC=CC=C1C(=O)O" 
            name="Aspirin" 
            height={600} 
            width="80%" 
          />
        </div>
      );
    
    case 'custom-bg':
      return (
        <div className="container mx-auto p-4">
          <h1 className="text-2xl font-bold mb-4">Custom Background Test</h1>
          <MolstarViewer 
            smiles="CC(=O)OC1=CC=CC=C1C(=O)O" 
            name="Aspirin" 
            height={400} 
            backgroundColor="#f0f9ff" 
          />
        </div>
      );
    
    case 'no-controls':
      return (
        <div className="container mx-auto p-4">
          <h1 className="text-2xl font-bold mb-4">No Controls Test</h1>
          <MolstarViewer 
            smiles="CC(=O)OC1=CC=CC=C1C(=O)O" 
            name="Aspirin" 
            height={400} 
            showControls={false} 
          />
        </div>
      );
    
    // Default: SMILES string example
    default:
      return (
        <div className="container mx-auto p-4">
          <h1 className="text-2xl font-bold mb-4">Molecule Viewer Test</h1>
          <MolstarViewer 
            smiles="CC(=O)OC1=CC=CC=C1C(=O)O" 
            name="Aspirin" 
            height={400} 
            style="ball-and-stick" 
          />
        </div>
      );
  }
}
*/
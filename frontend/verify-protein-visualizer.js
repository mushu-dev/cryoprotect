// Simple Node.js script to validate the Protein Visualizer component
// This doesn't require Jest or any complex testing frameworks

const fs = require('fs');
const path = require('path');

// Define color codes for output
const GREEN = '\x1b[32m';
const RED = '\x1b[31m';
const YELLOW = '\x1b[33m';
const RESET = '\x1b[0m';

console.log(`${YELLOW}Validating Protein Visualizer component...${RESET}`);

// Path to the component
const componentPath = path.join(
  __dirname, 
  'src/components/protein-visualizer/MolstarViewer.tsx'
);

try {
  // Read the component file
  const content = fs.readFileSync(componentPath, 'utf8');
  
  // Define validation checks
  const checks = [
    {
      name: 'Component exists',
      test: () => content.length > 0
    },
    {
      name: 'Exports MolstarViewer component',
      test: () => content.includes('export function MolstarViewer')
    },
    {
      name: 'Has proper interface definition',
      test: () => content.includes('interface MolstarViewerProps')
    },
    {
      name: 'Supports PDB ID input',
      test: () => content.includes('pdbId?:')
    },
    {
      name: 'Supports SMILES input',
      test: () => content.includes('smiles?:')
    },
    {
      name: 'Has style property',
      test: () => content.includes('style?:')
    },
    {
      name: 'Has visualization container',
      test: () => content.includes('protein-visualizer')
    },
    {
      name: 'Handles loading state',
      test: () => content.includes('isLoading') && content.includes('animate-spin')
    },
    {
      name: 'Handles error state',
      test: () => content.includes('error') && content.includes('bg-destructive')
    },
    {
      name: 'Supports styling controls',
      test: () => content.includes('showControls') && content.includes('SelectTrigger')
    },
    {
      name: 'Properly cleans up resources',
      test: () => content.includes('return () =>') && content.includes('dispose()')
    }
  ];
  
  // Run all checks
  let passedChecks = 0;
  let failedChecks = 0;
  
  console.log(`\n${YELLOW}Running validation checks:${RESET}\n`);
  
  checks.forEach((check, index) => {
    const result = check.test();
    
    if (result) {
      console.log(`${GREEN}✓ ${check.name}${RESET}`);
      passedChecks++;
    } else {
      console.log(`${RED}✗ ${check.name}${RESET}`);
      failedChecks++;
    }
  });
  
  // Print summary
  console.log(`\n${YELLOW}Validation summary:${RESET}`);
  console.log(`${GREEN}Passed: ${passedChecks}${RESET}`);
  console.log(`${RED}Failed: ${failedChecks}${RESET}`);
  
  if (failedChecks === 0) {
    console.log(`\n${GREEN}✓ All validation checks passed!${RESET}`);
    console.log(`${YELLOW}The component meets all basic requirements.${RESET}`);
    console.log(`${YELLOW}For a detailed manual test, run:${RESET}`);
    console.log(`npm run dev`);
    console.log(`Then navigate to: http://localhost:3000/protein-visualizer-demo`);
    process.exit(0);
  } else {
    console.log(`\n${RED}✗ Some validation checks failed.${RESET}`);
    console.log(`${YELLOW}Please review the component:${RESET}`);
    console.log(`${componentPath}`);
    process.exit(1);
  }
} catch (err) {
  console.error(`${RED}Error reading component file:${RESET}`, err);
  console.log(`${YELLOW}Please check if the file exists at:${RESET}`);
  console.log(componentPath);
  process.exit(1);
}
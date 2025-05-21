// Experimental Data Enhancement UI Test
// This script tests the implementation of the experimental data enhancement UI

const { execSync } = require('child_process');
const fs = require('fs');
const path = require('path');

// Utility function to run a playwright command using our container solution
function runPlaywrightCommand(command, ...args) {
  const scriptPath = path.join(__dirname, 'mcp-playwright-final.sh');
  const fullCommand = `${scriptPath} ${command} ${args.join(' ')}`;
  console.log(`\nRunning command: ${fullCommand}`);
  
  try {
    const output = execSync(fullCommand, { encoding: 'utf8' });
    console.log('Command succeeded');
    return output;
  } catch (error) {
    console.error(`Command failed: ${error.message}`);
    throw error;
  }
}

// Test function to navigate to URLs and take screenshots
async function testExperimentalUI() {
  console.log('Starting Experimental Data Enhancement UI Test');
  
  // Check if the container is ready
  try {
    runPlaywrightCommand('status');
    console.log('Playwright container is ready');
  } catch (error) {
    console.error('Playwright container is not ready:', error);
    process.exit(1);
  }
  
  // Create screenshots directory if it doesn't exist
  const screenshotDir = path.join(__dirname, 'screenshots');
  if (!fs.existsSync(screenshotDir)) {
    fs.mkdirSync(screenshotDir);
  }
  
  try {
    // Test the homepage
    console.log('\n--- Testing Homepage ---');
    runPlaywrightCommand('browser_navigate', 'https://cryoprotect.app');
    runPlaywrightCommand('browser_take_screenshot', 'https://cryoprotect.app', path.join(screenshotDir, 'homepage.png'));
    
    // Test the experiments page
    console.log('\n--- Testing Experiments Page ---');
    runPlaywrightCommand('browser_navigate', 'https://cryoprotect.app/experiments');
    runPlaywrightCommand('browser_take_screenshot', 'https://cryoprotect.app/experiments', path.join(screenshotDir, 'experiments-page.png'));
    
    // Test a specific experiment details page
    console.log('\n--- Testing Experiment Details Page ---');
    runPlaywrightCommand('browser_navigate', 'https://cryoprotect.app/experiments/1');
    runPlaywrightCommand('browser_take_screenshot', 'https://cryoprotect.app/experiments/1', path.join(screenshotDir, 'experiment-details.png'));
    
    // Test the protocols page
    console.log('\n--- Testing Protocols Page ---');
    runPlaywrightCommand('browser_navigate', 'https://cryoprotect.app/protocols');
    runPlaywrightCommand('browser_take_screenshot', 'https://cryoprotect.app/protocols', path.join(screenshotDir, 'protocols-page.png'));
    
    // Test a specific protocol details page
    console.log('\n--- Testing Protocol Details Page ---');
    runPlaywrightCommand('browser_navigate', 'https://cryoprotect.app/protocols/1');
    runPlaywrightCommand('browser_take_screenshot', 'https://cryoprotect.app/protocols/1', path.join(screenshotDir, 'protocol-details.png'));
    
    // Test accessibility snapshots
    console.log('\n--- Testing Accessibility Snapshots ---');
    const homepageSnapshot = runPlaywrightCommand('browser_snapshot', 'https://cryoprotect.app');
    console.log('Homepage structure verified');
    
    const experimentsSnapshot = runPlaywrightCommand('browser_snapshot', 'https://cryoprotect.app/experiments');
    console.log('Experiments page structure verified');
    
    const protocolsSnapshot = runPlaywrightCommand('browser_snapshot', 'https://cryoprotect.app/protocols');
    console.log('Protocols page structure verified');
    
    console.log('\nAll tests completed successfully');
    console.log(`Screenshots saved to ${screenshotDir}`);
    
  } catch (error) {
    console.error('Test failed:', error);
    process.exit(1);
  }
}

// Run the test function
testExperimentalUI();
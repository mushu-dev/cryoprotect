#!/bin/bash
# Quick verification script for the experimental data enhancement feature after deployment

# Configuration
NETLIFY_SITE_NAME=${NETLIFY_SITE_NAME:-cryoprotect}
NETLIFY_URL=${NETLIFY_URL:-"https://${NETLIFY_SITE_NAME}.netlify.app"}

echo "🔍 Verifying Experimental Data Enhancement Features"
echo "================================================="
echo "URL: $NETLIFY_URL"

# Create a temporary test file
TEST_FILE=$(mktemp)
cat > $TEST_FILE << EOF
const { chromium } = require('playwright');

async function verifyExperimentalFeatures() {
  console.log('Starting verification...');
  const browser = await chromium.launch();
  const page = await browser.newPage();
  
  try {
    // Go to the main page
    console.log('Navigating to: ${NETLIFY_URL}');
    await page.goto('${NETLIFY_URL}');
    console.log('✅ Main page loaded successfully');
    
    // Check if we can navigate to experiments page
    console.log('Navigating to experiments page...');
    await page.goto('${NETLIFY_URL}/experiments');
    
    // Wait for experiment list to appear
    console.log('Checking for experiment list...');
    const experimentListVisible = await page.waitForSelector('#experiment-list', { timeout: 10000 })
      .then(() => true)
      .catch(() => false);
    
    if (experimentListVisible) {
      console.log('✅ Experiment list loaded successfully');
    } else {
      console.log('❌ Experiment list not found');
    }
    
    // Check if experiment cards are visible
    const experimentCardsCount = await page.locator('.experiment-card').count();
    if (experimentCardsCount > 0) {
      console.log(\`✅ Found \${experimentCardsCount} experiment cards\`);
    } else {
      console.log('❌ No experiment cards found');
    }
    
    // Try to navigate to the first experiment
    if (experimentCardsCount > 0) {
      console.log('Attempting to navigate to experiment details...');
      await page.locator('.experiment-card').first().click();
      
      // Check if we're on the experiment detail page
      const detailPageLoaded = await page.waitForSelector('#experiment-details', { timeout: 10000 })
        .then(() => true)
        .catch(() => false);
        
      if (detailPageLoaded) {
        console.log('✅ Experiment detail page loaded successfully');
      } else {
        console.log('❌ Experiment detail page not found');
      }
    }
    
    // Check for experimental chart components
    const chartComponentsVisible = await page.locator('.experiment-chart').count();
    if (chartComponentsVisible > 0) {
      console.log(\`✅ Found \${chartComponentsVisible} chart components\`);
    } else {
      console.log('⚠️ No chart components found (may still be in progress)');
    }
    
    console.log('Verification complete!');
  } catch (error) {
    console.error('Error during verification:', error);
  } finally {
    await browser.close();
  }
}

verifyExperimentalFeatures();
EOF

# Run the verification
echo "Running verification checks..."
node $TEST_FILE

# Clean up
rm $TEST_FILE

echo ""
echo "For complete validation, run:"
echo "./validate-experimental-data-enhancement.sh"
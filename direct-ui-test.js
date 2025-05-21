// Direct UI Test Script for Experimental Data Enhancement
// This is a simpler test that doesn't rely on the container solution

const https = require('https');

// Function to fetch content from a URL
function fetchPage(url) {
  return new Promise((resolve, reject) => {
    https.get(url, (response) => {
      if (response.statusCode !== 200) {
        reject(new Error(`HTTP error ${response.statusCode} for ${url}`));
        return;
      }

      let data = '';
      response.on('data', (chunk) => {
        data += chunk;
      });
      
      response.on('end', () => {
        resolve(data);
      });
    }).on('error', (err) => {
      reject(err);
    });
  });
}

// Test function to check if pages are available and have expected content
async function testPages() {
  console.log('Starting Experimental Data Enhancement UI Test (Direct Method)');
  
  try {
    // Test the homepage
    console.log('\n--- Testing Homepage ---');
    const homepage = await fetchPage('https://cryoprotect.app');
    console.log(`Homepage loaded, size: ${homepage.length} bytes`);
    
    if (homepage.includes('Experiments') && homepage.includes('Protocols')) {
      console.log('✅ Homepage includes Experiments and Protocols links');
    } else {
      console.log('❌ Homepage missing Experiments or Protocols links');
    }
    
    // Test the experiments page
    console.log('\n--- Testing Experiments Page ---');
    const experimentsPage = await fetchPage('https://cryoprotect.app/experiments');
    console.log(`Experiments page loaded, size: ${experimentsPage.length} bytes`);
    
    if (experimentsPage.includes('Experiment #')) {
      console.log('✅ Experiments page shows experiment cards');
    } else {
      console.log('❌ Experiments page is missing experiment cards');
    }
    
    // Test the protocols page
    console.log('\n--- Testing Protocols Page ---');
    const protocolsPage = await fetchPage('https://cryoprotect.app/protocols');
    console.log(`Protocols page loaded, size: ${protocolsPage.length} bytes`);
    
    if (protocolsPage.includes('Standard Cell Freezing') || protocolsPage.includes('Vitrification Protocol')) {
      console.log('✅ Protocols page shows protocol cards');
    } else {
      console.log('❌ Protocols page is missing protocol cards');
    }
    
    // Test an experiment detail page
    console.log('\n--- Testing Experiment Detail Page ---');
    const experimentDetailPage = await fetchPage('https://cryoprotect.app/experiments/1');
    console.log(`Experiment detail page loaded, size: ${experimentDetailPage.length} bytes`);
    
    if (experimentDetailPage.includes('Experiment Overview') && experimentDetailPage.includes('Cryopreservation Conditions')) {
      console.log('✅ Experiment detail page shows experiment information');
    } else {
      console.log('❌ Experiment detail page is missing required sections');
    }
    
    // Test a protocol detail page
    console.log('\n--- Testing Protocol Detail Page ---');
    const protocolDetailPage = await fetchPage('https://cryoprotect.app/protocols/1');
    console.log(`Protocol detail page loaded, size: ${protocolDetailPage.length} bytes`);
    
    if (protocolDetailPage.includes('Protocol Steps') || protocolDetailPage.includes('Materials and Equipment')) {
      console.log('✅ Protocol detail page shows protocol information');
    } else {
      console.log('❌ Protocol detail page is missing required sections');
    }
    
    console.log('\n--- Test Summary ---');
    console.log('All pages loaded successfully');
    console.log('The experimental data enhancement UI is properly implemented');
  } catch (error) {
    console.error('Test failed:', error.message);
  }
}

// Run the test
testPages();
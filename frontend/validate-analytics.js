/**
 * Analytics validation script
 * This script tests the functionality of the analytics integration
 * 
 * Run with: node validate-analytics.js
 */

const { chromium } = require('@playwright/test');

async function validateAnalytics() {
  console.log('Starting analytics validation test...');
  const browser = await chromium.launch({ headless: false });
  const page = await browser.newPage();
  
  // Capture console logs
  const consoleLogs = [];
  page.on('console', (msg) => {
    const text = msg.text();
    if (text.includes('[Analytics]') || text.includes('[Netlify Analytics]')) {
      consoleLogs.push(text);
      console.log(`Console: ${text}`);
    }
  });
  
  // Navigate to analytics test page
  try {
    console.log('Loading analytics test page...');
    await page.goto('http://localhost:3000/analytics-test', { waitUntil: 'networkidle' });
    console.log('Page loaded successfully');
    
    // Check if analytics consent banner appears
    const consentBanner = await page.getByText('This site uses analytics');
    if (await consentBanner.isVisible()) {
      console.log('✅ Analytics consent banner detected');
      
      // Accept analytics
      console.log('Accepting analytics consent...');
      await page.getByRole('button', { name: 'Accept' }).click();
      console.log('Consent accepted');
    } else {
      console.log('⚠️ Analytics consent banner not found (may have been accepted already)');
    }
    
    // Test button click tracking
    console.log('Testing button click event tracking...');
    await page.getByRole('button', { name: 'Track Button Click Event' }).click();
    await page.waitForTimeout(500); // Wait for analytics to process
    
    // Test feature usage tracking
    console.log('Testing feature usage tracking...');
    await page.getByRole('button', { name: 'Track Feature Usage' }).click();
    await page.waitForTimeout(500); // Wait for analytics to process
    
    // Check if we captured analytics events
    const hasPageView = consoleLogs.some(log => log.includes('Page view'));
    const hasButtonClick = consoleLogs.some(log => log.includes('test_button_click'));
    const hasFeatureUsage = consoleLogs.some(log => log.includes('analytics_test'));
    
    console.log('\nTest Results:');
    console.log(`Page view tracking: ${hasPageView ? '✅ Detected' : '❌ Not detected'}`);
    console.log(`Button click tracking: ${hasButtonClick ? '✅ Detected' : '❌ Not detected'}`);
    console.log(`Feature usage tracking: ${hasFeatureUsage ? '✅ Detected' : '❌ Not detected'}`);
    
    if (hasPageView && hasButtonClick && hasFeatureUsage) {
      console.log('\n✅ Analytics integration is working correctly!');
    } else {
      console.log('\n⚠️ Some analytics events were not detected. Check the implementation.');
    }
    
    // Take a screenshot as evidence
    await page.screenshot({ path: 'analytics-test-screenshot.png' });
    console.log('Screenshot saved as analytics-test-screenshot.png');
    
  } catch (error) {
    console.error('Error during validation:', error);
  } finally {
    await page.waitForTimeout(2000); // Give a moment to see the results
    await browser.close();
    console.log('Test completed');
  }
}

validateAnalytics();
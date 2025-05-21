/**
 * Client-side analytics testing utilities
 * This script provides utility functions for testing analytics in the browser
 * 
 * Usage:
 * 1. Import in a component: import { testAnalytics } from '../utils/analytics-test-client';
 * 2. Call testAnalytics() to run diagnostics
 */

/**
 * Test analytics implementation in the browser
 * @returns {Object} Test results
 */
export function testAnalytics() {
  console.log("📊 Testing analytics integration...");
  
  const results = {
    consent: false,
    provider: false,
    identity: false,
    localStorage: false,
    events: false,
  };
  
  // Check if analytics provider is loaded
  if (window.__ANALYTICS_PROVIDER_LOADED) {
    console.log("✅ Analytics provider loaded");
    results.provider = true;
  } else {
    console.log("❌ Analytics provider not loaded");
  }
  
  // Check for Netlify Identity
  if (window.netlifyIdentity) {
    console.log("✅ Netlify Identity available");
    results.identity = true;
  } else {
    console.log("❌ Netlify Identity not available");
  }
  
  // Check localStorage for consent
  const consentStatus = localStorage.getItem('analytics-consent');
  if (consentStatus !== null) {
    console.log(`✅ Consent status in localStorage: ${consentStatus}`);
    results.localStorage = true;
    results.consent = consentStatus === 'true';
  } else {
    console.log("❌ Consent status not found in localStorage");
  }
  
  // Test event dispatching
  try {
    // Create and dispatch a test event
    const testEvent = new CustomEvent('test-analytics', { 
      detail: { timestamp: new Date().toISOString() } 
    });
    document.dispatchEvent(testEvent);
    console.log("✅ Test event dispatched");
    results.events = true;
  } catch (error) {
    console.log("❌ Error dispatching test event:", error);
  }
  
  console.log("📊 Analytics test complete.");
  return results;
}

/**
 * Manually track an analytics event (for testing)
 * @param {string} eventName - Name of the event to track
 * @param {Object} properties - Event properties
 */
export function trackTestEvent(eventName, properties = {}) {
  try {
    // Access the analytics context if available
    const analyticsContext = window.__ANALYTICS_CONTEXT;
    if (analyticsContext && analyticsContext.trackEvent) {
      analyticsContext.trackEvent(eventName, properties);
      console.log(`✅ Event tracked via context: ${eventName}`, properties);
      return true;
    }
    
    // Try fallback to global trackEvent method
    if (window.trackAnalyticsEvent) {
      window.trackAnalyticsEvent(eventName, properties);
      console.log(`✅ Event tracked via global method: ${eventName}`, properties);
      return true;
    }
    
    // Direct call to Netlify analytics if available
    if (window._nenl && window._nenl.push) {
      window._nenl.push(['track', eventName, properties]);
      console.log(`✅ Event tracked directly with Netlify: ${eventName}`, properties);
      return true;
    }
    
    console.log("❌ Unable to track event: no tracking method available");
    return false;
  } catch (error) {
    console.error("❌ Error tracking test event:", error);
    return false;
  }
}

/**
 * Toggle analytics consent (for testing)
 * @param {boolean} enabled - Whether to enable or disable analytics
 */
export function toggleAnalyticsConsent(enabled) {
  localStorage.setItem('analytics-consent', String(enabled));
  console.log(`✅ Analytics consent set to: ${enabled}`);
  
  // Try to notify the analytics system
  try {
    const analyticsContext = window.__ANALYTICS_CONTEXT;
    if (analyticsContext && analyticsContext.setEnabled) {
      analyticsContext.setEnabled(enabled);
      console.log("✅ Analytics context notified of consent change");
    } else {
      console.log("ℹ️ Analytics context not available, reload page to apply changes");
    }
  } catch (error) {
    console.error("❌ Error updating analytics consent:", error);
  }
  
  return true;
}

/**
 * Clear analytics consent (for testing)
 */
export function clearAnalyticsConsent() {
  localStorage.removeItem('analytics-consent');
  console.log("✅ Analytics consent cleared from localStorage");
  console.log("ℹ️ Reload the page to see the consent banner");
  return true;
}
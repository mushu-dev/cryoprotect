/**
 * Netlify analytics integration script
 * 
 * This script provides a React component to integrate analytics on Netlify sites.
 * It works with both Netlify Analytics and optionally other analytics providers.
 */

'use client';

import Script from 'next/script';
import { useEffect } from 'react';

/**
 * Custom Netlify analytics integration supporting page view tracking
 */
function trackPageView() {
  if (typeof window !== 'undefined') {
    // Only run on Netlify environment
    const isNetlify = process.env.NEXT_PUBLIC_NETLIFY === 'true' || 
                    window.location.hostname.includes('netlify.app');
    
    if (!isNetlify) return;
    
    // Netlify Analytics is automatically enabled if you've purchased it
    // but we can add additional tracking functionality here
    
    console.log('[Netlify Analytics] Page view:', window.location.pathname);
    
    // Custom event tracking if needed
    if (window.plausible) {
      window.plausible('pageview');
    }
  }
}

/**
 * NetlifyAnalytics component to add analytics tracking to your app
 * This component should be used in the root layout
 */
export function NetlifyAnalytics() {
  useEffect(() => {
    // Track initial page view
    trackPageView();
    
    // Track page changes
    const handleRouteChange = () => {
      trackPageView();
    };
    
    // Add event listeners for client-side navigation
    window.addEventListener('popstate', handleRouteChange);
    
    // Clean up event listeners
    return () => {
      window.removeEventListener('popstate', handleRouteChange);
    };
  }, []);
  
  return (
    <>
      {/* Optional: Add Plausible analytics for enhanced tracking */}
      {process.env.NEXT_PUBLIC_USE_PLAUSIBLE === 'true' && (
        <Script
          src="https://plausible.io/js/script.js"
          data-domain={process.env.NEXT_PUBLIC_PLAUSIBLE_DOMAIN || window.location.hostname}
          strategy="afterInteractive"
        />
      )}
    </>
  );
}

/**
 * Helper function to track custom events
 * @param {string} eventName - Name of the event to track
 * @param {object} properties - Additional event properties
 */
export function trackEvent(eventName, properties = {}) {
  if (typeof window === 'undefined') return;
  
  // Log event to console in development
  if (process.env.NODE_ENV === 'development') {
    console.log(`[Analytics Event] ${eventName}`, properties);
  }
  
  // Track with Plausible if available
  if (window.plausible) {
    window.plausible(eventName, { props: properties });
  }
}

export default NetlifyAnalytics;
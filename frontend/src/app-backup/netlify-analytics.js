/**
 * Simple Netlify analytics tracking module
 * Provides basic event tracking for Netlify analytics
 */

/**
 * Track an event with Netlify analytics
 * @param {string} eventName - Name of the event to track
 * @param {object} properties - Additional properties for the event
 */
export function trackEvent(eventName, properties = {}) {
  // Skip tracking in development environment
  if (process.env.NODE_ENV === 'development') {
    console.log(`[Analytics] Event tracked: ${eventName}`, properties);
    return;
  }

  // Track event with Netlify analytics if available
  if (typeof window !== 'undefined' && window.netlifyIdentity) {
    try {
      // Note: This is a simplified version, actual implementation may vary
      if (eventName === 'pageview') {
        // Track pageview
        window.netlifyIdentity.on('pageview', properties.path);
      } else {
        // Track custom event
        window.netlifyIdentity.on('event', {
          type: eventName,
          ...properties
        });
      }
    } catch (error) {
      console.error('[Analytics] Error tracking event:', error);
    }
  }
}

/**
 * Initialize Netlify analytics
 */
export function initNetlifyAnalytics() {
  if (typeof window !== 'undefined' && !window.netlifyIdentityInitialized) {
    const script = document.createElement('script');
    script.src = 'https://identity.netlify.com/v1/netlify-identity-widget.js';
    script.async = true;
    script.onload = () => {
      window.netlifyIdentityInitialized = true;
    };
    document.head.appendChild(script);
  }
}
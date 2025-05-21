/**
 * Netlify Analytics integration
 * This file provides utility functions for sending events to Netlify Analytics
 */

/**
 * Track an event in Netlify Analytics
 * @param eventName - Name of the event
 * @param properties - Additional properties to track
 */
export function trackEvent(eventName: string, properties: Record<string, any> = {}) {
  // Log in development
  if (process.env.NODE_ENV === 'development') {
    console.log(`[Netlify Analytics] Event: ${eventName}`, properties);
    return;
  }

  // Make sure the window object exists
  if (typeof window === 'undefined') return;

  // Check if Netlify Identity is available
  if (window.netlifyIdentity) {
    try {
      // Track the event
      window.netlifyIdentity.on('init', (user) => {
        const userId = user ? user.id : 'anonymous';
        const eventData = {
          ...properties,
          userId,
          timestamp: new Date().toISOString(),
        };
        
        // Send to Netlify Analytics
        if (window._nenl) {
          window._nenl.push(['track', eventName, eventData]);
        }
      });
    } catch (error) {
      console.error('[Netlify Analytics] Error tracking event:', error);
    }
  } else {
    // Netlify Identity not loaded yet - use queue
    window.netlifyAnalyticsQueue = window.netlifyAnalyticsQueue || [];
    window.netlifyAnalyticsQueue.push({
      eventName,
      properties,
    });
  }
}

// Add to window type
declare global {
  interface Window {
    netlifyIdentity: any;
    _nenl: any;
    netlifyAnalyticsQueue: Array<{
      eventName: string;
      properties: Record<string, any>;
    }>;
  }
}
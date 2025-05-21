import React, { createContext, useContext, useEffect, useState } from 'react';
import { useAnalytics } from '../../hooks/useAnalytics';

// Create the context with default values
const AnalyticsContext = createContext({
  isEnabled: false,
  trackEvent: () => {},
  trackPageView: () => {},
  trackFeatureUsage: () => {},
  setEnabled: () => {},
});

// Hook to use the analytics context
export const useAnalyticsContext = () => useContext(AnalyticsContext);

// Analytics Provider Component
export function AnalyticsProvider({ children }) {
  const analytics = useAnalytics();
  const [userConsent, setUserConsent] = useState(null);
  
  // Initialize user consent from localStorage if available
  useEffect(() => {
    if (typeof window !== 'undefined') {
      const storedConsent = localStorage.getItem('analytics-consent');
      setUserConsent(storedConsent === 'true');
    }
  }, []);
  
  // Determine if analytics is enabled based on environment and user consent
  const isEnabled = analytics.isAnalyticsEnabled && userConsent !== false;
  
  // Track page views on route changes
  useEffect(() => {
    if (isEnabled && typeof window !== 'undefined') {
      // Get current path
      const path = window.location.pathname + window.location.search;
      analytics.trackPageView(path);
    }
  }, [isEnabled, analytics]);
  
  // Update user consent
  const setEnabled = (enabled) => {
    setUserConsent(enabled);
    if (typeof window !== 'undefined') {
      localStorage.setItem('analytics-consent', String(enabled));
    }
  };
  
  // Value for the context provider
  const contextValue = {
    isEnabled,
    trackEvent: (eventName, properties) => {
      if (isEnabled) {
        analytics.trackAction(eventName, 'custom', properties);
      }
    },
    trackPageView: analytics.trackPageView,
    trackFeatureUsage: analytics.trackFeatureUsage,
    setEnabled,
  };
  
  return (
    <AnalyticsContext.Provider value={contextValue}>
      {children}
    </AnalyticsContext.Provider>
  );
}
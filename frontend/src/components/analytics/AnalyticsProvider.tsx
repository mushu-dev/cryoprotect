'use client';

import React, { createContext, useContext, useEffect, useState } from 'react';
import { usePathname, useSearchParams } from 'next/navigation';
import { useAnalytics } from '@/hooks/useAnalytics';

// Define the context type
type AnalyticsContextType = {
  isEnabled: boolean;
  trackEvent: (eventName: string, properties?: Record<string, any>) => void;
  trackPageView: (path?: string) => void;
  trackFeatureUsage: (feature: string, properties?: Record<string, any>) => void;
  // Enable/disable analytics (for user preference compliance)
  setEnabled: (enabled: boolean) => void;
};

// Create the context with default values
const AnalyticsContext = createContext<AnalyticsContextType>({
  isEnabled: false,
  trackEvent: () => {},
  trackPageView: () => {},
  trackFeatureUsage: () => {},
  setEnabled: () => {},
});

// Hook to use the analytics context
export const useAnalyticsContext = () => useContext(AnalyticsContext);

// Analytics Provider Component
export function AnalyticsProvider({ children }: { children: React.ReactNode }) {
  const pathname = usePathname();
  const searchParams = useSearchParams();
  const analytics = useAnalytics();
  const [userConsent, setUserConsent] = useState<boolean | null>(null);
  
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
    if (isEnabled && pathname) {
      // Include search params if present
      const path = searchParams.toString() 
        ? `${pathname}?${searchParams.toString()}`
        : pathname;
        
      analytics.trackPageView(path);
    }
  }, [pathname, searchParams, isEnabled, analytics]);
  
  // Update user consent
  const setEnabled = (enabled: boolean) => {
    setUserConsent(enabled);
    if (typeof window !== 'undefined') {
      localStorage.setItem('analytics-consent', String(enabled));
    }
  };
  
  // Value for the context provider
  const contextValue: AnalyticsContextType = {
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
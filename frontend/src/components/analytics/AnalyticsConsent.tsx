'use client';

import { useState, useEffect } from 'react';
import { useAnalyticsContext } from './AnalyticsProvider';

export function AnalyticsConsent() {
  const { isEnabled, setEnabled } = useAnalyticsContext();
  const [isVisible, setIsVisible] = useState(false);
  const [hasInteracted, setHasInteracted] = useState(false);
  
  // Show consent banner if user hasn't made a choice yet
  useEffect(() => {
    const hasConsented = localStorage.getItem('analytics-consent');
    if (hasConsented === null) {
      // Show banner after a short delay
      const timer = setTimeout(() => {
        setIsVisible(true);
      }, 1500);
      
      return () => clearTimeout(timer);
    } else {
      setHasInteracted(true);
    }
  }, []);
  
  const handleAccept = () => {
    setEnabled(true);
    setIsVisible(false);
    setHasInteracted(true);
  };
  
  const handleDecline = () => {
    setEnabled(false);
    setIsVisible(false);
    setHasInteracted(true);
  };
  
  // Don't render anything if the banner shouldn't be visible
  if (!isVisible) {
    return null;
  }
  
  return (
    <div className="fixed bottom-4 left-4 right-4 bg-white dark:bg-gray-800 p-4 rounded-lg shadow-lg z-50 max-w-2xl mx-auto border border-gray-200 dark:border-gray-700">
      <div className="flex flex-col sm:flex-row gap-4 items-start sm:items-center justify-between">
        <div className="flex-1">
          <h3 className="font-medium text-sm">This site uses analytics</h3>
          <p className="text-sm text-gray-600 dark:text-gray-300 mt-1">
            We use analytics to improve your experience. This helps us understand how you use our platform.
          </p>
        </div>
        <div className="flex gap-2 mt-2 sm:mt-0">
          <button
            onClick={handleDecline}
            className="text-sm px-3 py-1.5 border border-gray-300 dark:border-gray-600 rounded-md hover:bg-gray-100 dark:hover:bg-gray-700 transition-colors"
          >
            Decline
          </button>
          <button
            onClick={handleAccept}
            className="text-sm px-3 py-1.5 bg-blue-600 text-white rounded-md hover:bg-blue-700 transition-colors"
          >
            Accept
          </button>
        </div>
      </div>
    </div>
  );
}

// Settings toggle component for analytics preferences
export function AnalyticsToggle() {
  const { isEnabled, setEnabled } = useAnalyticsContext();
  
  return (
    <div className="flex items-center justify-between py-3">
      <div>
        <h3 className="text-sm font-medium">Analytics Tracking</h3>
        <p className="text-sm text-gray-500 dark:text-gray-400">
          Allow anonymous usage data collection to help improve our service
        </p>
      </div>
      <button
        type="button"
        role="switch"
        aria-checked={isEnabled}
        className={`relative inline-flex h-6 w-11 items-center rounded-full transition-colors focus:outline-none focus:ring-2 focus:ring-blue-500 focus:ring-offset-2 ${
          isEnabled ? 'bg-blue-600' : 'bg-gray-200 dark:bg-gray-700'
        }`}
        onClick={() => setEnabled(!isEnabled)}
      >
        <span className="sr-only">Enable analytics</span>
        <span
          className={`inline-block h-4 w-4 transform rounded-full bg-white transition-transform ${
            isEnabled ? 'translate-x-6' : 'translate-x-1'
          }`}
        />
      </button>
    </div>
  );
}
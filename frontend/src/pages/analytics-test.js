import { useEffect, useState } from 'react';
import { useAnalyticsContext } from '../components/analytics/AnalyticsProvider';
import { AnalyticsToggle } from '../components/analytics/AnalyticsConsent';

export default function AnalyticsTestPage() {
  const { isEnabled, trackEvent, trackFeatureUsage } = useAnalyticsContext();
  const [eventCount, setEventCount] = useState(0);
  
  // Track page view on mount
  useEffect(() => {
    if (isEnabled) {
      console.log('Analytics page view tracked automatically');
    }
  }, [isEnabled]);
  
  const handleTrackButtonClick = () => {
    const newCount = eventCount + 1;
    setEventCount(newCount);
    trackEvent('test_button_click', { count: newCount });
    console.log('Custom event tracked:', 'test_button_click', { count: newCount });
  };
  
  const handleTrackFeatureUsage = () => {
    trackFeatureUsage('analytics_test', { timestamp: new Date().toISOString() });
    console.log('Feature usage tracked:', 'analytics_test');
  };

  return (
    <div className="container mx-auto px-4 py-8">
      <h1 className="text-2xl font-bold mb-4">Analytics Test Page</h1>
      
      <div className="bg-white dark:bg-gray-800 p-6 rounded-lg shadow-md mb-6">
        <h2 className="text-xl font-semibold mb-4">Analytics Status</h2>
        <p className="mb-4">
          Analytics tracking is currently: <span className="font-medium">{isEnabled ? 'Enabled' : 'Disabled'}</span>
        </p>
        
        <div className="mb-6">
          <AnalyticsToggle />
        </div>
        
        <div className="space-y-4">
          <h3 className="text-lg font-medium">Test Analytics Events</h3>
          
          <div>
            <button
              onClick={handleTrackButtonClick}
              disabled={!isEnabled}
              className={`px-4 py-2 rounded-md transition-colors ${
                isEnabled 
                  ? 'bg-blue-600 text-white hover:bg-blue-700' 
                  : 'bg-gray-300 text-gray-500 cursor-not-allowed'
              }`}
            >
              Track Button Click Event
            </button>
            {eventCount > 0 && (
              <p className="mt-2 text-sm text-gray-600 dark:text-gray-400">
                Button clicked {eventCount} time{eventCount !== 1 ? 's' : ''}
              </p>
            )}
          </div>
          
          <div>
            <button
              onClick={handleTrackFeatureUsage}
              disabled={!isEnabled}
              className={`px-4 py-2 rounded-md transition-colors ${
                isEnabled 
                  ? 'bg-green-600 text-white hover:bg-green-700' 
                  : 'bg-gray-300 text-gray-500 cursor-not-allowed'
              }`}
            >
              Track Feature Usage
            </button>
          </div>
        </div>
      </div>
      
      <div className="bg-white dark:bg-gray-800 p-6 rounded-lg shadow-md">
        <h2 className="text-xl font-semibold mb-4">Analytics Events Monitor</h2>
        <p className="text-sm text-gray-600 dark:text-gray-400">
          Open your browser console to see the analytics events being logged.
        </p>
        
        <div className="mt-4 p-4 bg-gray-100 dark:bg-gray-700 rounded-md">
          <pre className="text-xs overflow-auto max-h-48">
            {`
// Example analytics events:
[Analytics] Page view: /analytics-test
[Analytics] Event: test_button_click { count: 1 }
[Analytics] Event: feature_used { feature: "analytics_test", timestamp: "2025-05-20T12:34:56.789Z" }
            `}
          </pre>
        </div>
      </div>
    </div>
  );
}
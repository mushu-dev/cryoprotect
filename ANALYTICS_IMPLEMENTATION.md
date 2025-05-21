# Analytics Implementation Documentation

This document outlines the implementation of analytics tracking in the CryoProtect application.

## Overview

The CryoProtect application uses a provider-agnostic analytics system that currently integrates with Netlify Analytics. The implementation includes user consent management, event tracking, and page view tracking across the application.

## Key Components

### 1. Analytics Provider

The `AnalyticsProvider` component (`src/components/analytics/AnalyticsProvider.tsx`) serves as the central hub for all analytics functionality:

- Exposes a React context for analytics functionality
- Manages user consent (stored in localStorage)
- Handles automatic page view tracking
- Provides methods for custom event tracking
- Integrates with the analytics implementation

### 2. Analytics Hook

The `useAnalytics` hook (`src/hooks/useAnalytics.ts`) provides the core tracking functionality:

- Determines if analytics is enabled based on environment variables
- Provides tracking methods for different types of events:
  - Page views
  - Feature usage
  - Errors
  - Searches
  - Custom actions
- Logs events to console in development mode for easier debugging

### 3. User Consent Management

The `AnalyticsConsent` component (`src/components/analytics/AnalyticsConsent.tsx`) handles user consent:

- Shows a consent banner when users first visit the site
- Allows users to accept or decline analytics tracking
- Stores user preferences in localStorage
- Provides a settings toggle component for user preference management

### 4. Netlify Analytics Integration

The Netlify Analytics integration (`src/app/netlify-analytics.ts`) handles sending events to Netlify:

- Automatically tracks page views
- Supports custom event tracking
- Handles development vs. production environments
- Works with Netlify Identity for user tracking
- Implements queue system for events when Netlify Identity is loading

## Implementation Details

### Integration in the Application

The analytics system is integrated at the application root level in `_app.js`:

1. The `AnalyticsProvider` wraps the entire application
2. The Netlify Identity script is loaded dynamically when enabled
3. The `AnalyticsConsent` component is rendered to manage user consent
4. Page views are automatically tracked on route changes
5. CSP headers in `netlify.toml` are configured to allow analytics scripts

### Environment Configuration

Analytics can be enabled/disabled via environment variables:

- `NEXT_PUBLIC_NETLIFY=true` - Enables Netlify Analytics
- `NEXT_PUBLIC_VERCEL=true` - Enables Vercel Analytics (future implementation)

These variables can be set in `.env.production` or in the Netlify/Vercel deployment settings.

## Testing Analytics

A test page and validation script are provided to verify the analytics implementation:

1. **Test Page**: `/analytics-test` provides a UI to test analytics features
   - Shows current analytics status
   - Allows toggling analytics consent
   - Provides buttons to test various tracking events
   - Displays analytics events in the console

2. **Validation Script**: `validate-analytics.js` automatically tests the implementation
   - Navigates to the test page
   - Interacts with analytics consent features
   - Triggers various tracking events
   - Captures and validates analytics events in the console

To test analytics:

```bash
# Start the application
cd frontend && npm run dev

# In another terminal, run the validation script
cd frontend && node validate-analytics.js
```

## Adding Custom Event Tracking

To track custom events in your components:

```typescript
import { useAnalyticsContext } from '../components/analytics/AnalyticsProvider';

function MyComponent() {
  const { trackEvent, trackFeatureUsage } = useAnalyticsContext();

  // Track a custom event
  const handleAction = () => {
    trackEvent('button_click', { component: 'MyComponent', action: 'save' });
  };

  // Track feature usage
  const handleFeatureUse = () => {
    trackFeatureUsage('advanced_search', { query_type: 'molecular' });
  };

  return (
    <div>
      <button onClick={handleAction}>Save</button>
      <button onClick={handleFeatureUse}>Advanced Search</button>
    </div>
  );
}
```

## Best Practices

1. **User Privacy**: Always respect user consent settings
2. **Performance**: Keep tracking code lightweight to avoid impacting user experience
3. **Data Collection**: Only collect necessary data that provides actionable insights
4. **Debugging**: Use console logging in development to verify events are firing correctly
5. **Documentation**: Document what events are being tracked and why
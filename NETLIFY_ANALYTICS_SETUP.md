# Netlify Analytics Setup Guide

This guide explains how to setup and use analytics with your Netlify deployment of the CryoProtect frontend application.

## Overview

We've implemented a dual-track analytics approach:

1. **Netlify Analytics**: Server-side analytics provided as a Netlify add-on service
2. **Custom Client Analytics**: Custom event tracking in the application

## Netlify Analytics Setup

Netlify Analytics is a paid add-on that provides privacy-focused, server-side analytics without affecting site performance.

### Enabling Netlify Analytics

1. Go to your Netlify site dashboard
2. Navigate to the "Analytics" tab
3. Click "Enable Netlify Analytics"
4. Confirm the purchase ($9/month per site)

Once enabled, Netlify Analytics works automatically without any code changes. It provides:

- Page views
- Top pages
- Visitor locations
- Device and browser data
- Bandwidth usage
- 404 tracking
- No cookies or client-side JavaScript required

## Custom Analytics Integration

We've also implemented a custom analytics system that provides more detailed tracking while respecting user privacy preferences.

### Features

- **User Consent**: Displays a consent banner and respects user choices
- **Privacy Controls**: Users can toggle analytics in settings
- **Page View Tracking**: Automatically tracks page navigation
- **Event Tracking**: API for tracking custom events
- **Feature Usage**: Track specific feature utilization
- **Error Logging**: Monitor application errors
- **Search Analytics**: Track search patterns and effectiveness

### Usage in Code

#### Tracking Page Views

Page views are automatically tracked by the `AnalyticsProvider`. For custom page tracking:

```typescript
import { useAnalyticsContext } from '@/components/analytics/AnalyticsProvider';

function MyComponent() {
  const { trackPageView } = useAnalyticsContext();
  
  // Track a custom page view
  trackPageView('/custom-path');
  
  return <div>...</div>;
}
```

#### Tracking Events

```typescript
import { useAnalyticsContext } from '@/components/analytics/AnalyticsProvider';

function MyComponent() {
  const { trackEvent } = useAnalyticsContext();
  
  const handleButtonClick = () => {
    // Track a custom event
    trackEvent('button_clicked', { buttonId: 'save-profile' });
    
    // Rest of your handler
  };
  
  return <button onClick={handleButtonClick}>Save Profile</button>;
}
```

#### Tracking Feature Usage

```typescript
import { useAnalyticsContext } from '@/components/analytics/AnalyticsProvider';

function FeatureComponent() {
  const { trackFeatureUsage } = useAnalyticsContext();
  
  const handleFeatureActivation = () => {
    trackFeatureUsage('molecular_viewer', { 
      mode: '3d', 
      molecule_id: 'ABC123' 
    });
    
    // Feature logic
  };
  
  return <div>...</div>;
}
```

## Optional Enhancements

### Enabling Plausible Analytics

We've included optional support for Plausible Analytics, a privacy-focused alternative:

1. Sign up for a Plausible account at [plausible.io](https://plausible.io)
2. Add your site domain in the Plausible dashboard
3. Set environment variables in Netlify:
   - `NEXT_PUBLIC_USE_PLAUSIBLE`: Set to `"true"`
   - `NEXT_PUBLIC_PLAUSIBLE_DOMAIN`: Your domain (e.g., `"cryoprotect.netlify.app"`)
4. Deploy your site

Alternatively, build with Plausible enabled:

```bash
npm run analytics:enable-plausible
```

## Static Export Compatibility

For static exports with Netlify (using `output: 'export'` in Next.js), we've ensured analytics remains compatible:

1. All analytics calls are wrapped in `useEffect` to ensure they only run client-side
2. We use a global event listener for navigation changes with `popstate`
3. The CSP in netlify.toml allows connections to Netlify analytics domains
4. Environment variables are passed through the build process for static generation
5. The `unoptimized: true` setting for images ensures compatibility with static export

For more details on handling client-side analytics with dynamic routes and static export, see [NETLIFY_DYNAMIC_ROUTES.md](./NETLIFY_DYNAMIC_ROUTES.md).

## File Structure

The analytics implementation consists of:

- `/src/app/netlify-analytics.js`: Base analytics integration
- `/src/hooks/useAnalytics.ts`: Hook for analytics functionality
- `/src/components/analytics/AnalyticsProvider.tsx`: Context provider
- `/src/components/analytics/AnalyticsConsent.tsx`: Consent UI components
- `/src/components/analytics/index.ts`: Convenient exports

## Environment Variables

The following environment variables control analytics behavior:

- `NEXT_PUBLIC_NETLIFY`: Set to `"true"` on Netlify deployments
- `NEXT_PUBLIC_VERCEL`: Set to `"true"` on Vercel deployments
- `NEXT_PUBLIC_USE_PLAUSIBLE`: Enable Plausible integration
- `NEXT_PUBLIC_PLAUSIBLE_DOMAIN`: Domain for Plausible

## Testing Analytics

1. Deploy to Netlify
2. Open the browser console
3. Navigate through the application
4. Look for "[Analytics] Page view:" logs in development mode
5. Check Netlify Analytics dashboard after 24 hours for data

## Privacy Considerations

Our analytics implementation:

- Respects user preferences with opt-in consent
- Stores consent in localStorage
- Provides settings controls for users
- Does not track personally identifiable information
- Complies with GDPR, CCPA, and other privacy regulations
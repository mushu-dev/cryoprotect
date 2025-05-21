# Vercel Analytics & Speed Insights Setup Guide

This guide explains how Vercel Analytics and Speed Insights have been implemented in the CryoProtect application.

## Integration Steps Completed

### Vercel Analytics

The following steps have been completed to add Vercel Analytics to the project:

1. Installed the `@vercel/analytics` package:
   ```bash
   npm install @vercel/analytics
   ```

2. Added the Analytics component to the root layout file at `frontend/src/app/layout.tsx`:
   ```jsx
   import { Analytics } from '@vercel/analytics/react'
   
   // ...inside the RootLayout component
   return (
     <html lang="en">
       <body>
         {/* ... other components */}
         <Analytics />
       </body>
     </html>
   )
   ```

3. Enabled analytics in the Vercel deployment script by adding the environment variable:
   ```bash
   -e VERCEL_ANALYTICS_ID=true
   ```

### Vercel Speed Insights

The following steps have been completed to add Vercel Speed Insights to the project:

1. Installed the `@vercel/speed-insights` package:
   ```bash
   npm install @vercel/speed-insights
   ```

2. Added the SpeedInsights component to the root layout file at `frontend/src/app/layout.tsx`:
   ```jsx
   import { SpeedInsights } from '@vercel/speed-insights/next'
   
   // ...inside the RootLayout component
   return (
     <html lang="en">
       <body>
         {/* ... other components */}
         <SpeedInsights />
       </body>
     </html>
   )
   ```

3. Enabled Speed Insights in the Vercel deployment script by adding the environment variable:
   ```bash
   -e VERCEL_SPEED_INSIGHTS=true
   ```

## Verifying the Setup

### Verifying Analytics

To verify that Vercel Analytics is working properly:

1. Deploy your application using the updated deployment script:
   ```bash
   cd frontend && ./deploy-to-vercel.sh
   ```

2. After deployment, visit your Vercel dashboard:
   - Go to https://vercel.com/dashboard
   - Select your CryoProtect project
   - Click on "Analytics" in the navigation

3. Verify that you can see the "Web Analytics" panel with data about page views, visitors, etc.

4. If you don't see data immediately, don't worry. It might take some time for analytics data to populate.

### Verifying Speed Insights

To verify that Speed Insights is working properly:

1. After deploying, visit your Vercel dashboard:
   - Go to https://vercel.com/dashboard
   - Select your CryoProtect project
   - Click on "Speed Insights" in the navigation

2. You should see performance metrics such as:
   - Core Web Vitals (LCP, FID, CLS)
   - Page load performance stats
   - Breakdown by page and device

3. As with Analytics, data may take some time to appear after users start visiting your site.

## Debugging Issues

### Analytics Issues

If analytics data isn't appearing:

1. Check browser console for any errors related to Vercel Analytics
2. Verify that the `VERCEL_ANALYTICS_ID` environment variable is set to `true` in your deployment
3. Make sure the `<Analytics />` component is correctly placed in your layout file
4. Use your browser's network tab to confirm that analytics requests are being sent to Vercel

### Speed Insights Issues

If Speed Insights data isn't appearing:

1. Check browser console for any errors related to Speed Insights
2. Verify that the `VERCEL_SPEED_INSIGHTS` environment variable is set to `true` in your deployment
3. Make sure the `<SpeedInsights />` component is correctly placed in your layout file
4. Check if you need to enable Speed Insights in your Vercel project settings

## Custom Events and Measurements

### Analytics Custom Events

You can track custom events with Vercel Analytics:

```javascript
import { track } from '@vercel/analytics';

// Track a custom event
track('user-login');

// Track with custom properties
track('search-executed', { query: 'cryoprotectant', results: 42 });
```

### Speed Insights Custom Measurements

You can record custom performance measurements with Speed Insights:

```javascript
import { SpeedInsights } from '@vercel/speed-insights/next';
import { useSpeedInsights } from '@vercel/speed-insights/react';

// Within a component using the hook:
function SearchComponent() {
  const { trackCustomMetric } = useSpeedInsights();
  
  const handleSearch = async (query) => {
    const startTime = performance.now();
    const results = await searchApi(query);
    const endTime = performance.now();
    
    // Track the search time
    trackCustomMetric('search-time', endTime - startTime);
    
    return results;
  };
  
  // ...component implementation
}
```

## Learn More

For more information on Vercel tools, visit:

### Analytics
- [Vercel Analytics Documentation](https://vercel.com/docs/analytics)
- [Web Analytics API Reference](https://vercel.com/docs/concepts/analytics/api)
- [Custom Events](https://vercel.com/docs/concepts/analytics/custom-events)

### Speed Insights
- [Speed Insights Documentation](https://vercel.com/docs/speed-insights)
- [Core Web Vitals](https://vercel.com/docs/speed-insights/web-vitals)
- [Custom Metrics](https://vercel.com/docs/speed-insights/custom-metrics)
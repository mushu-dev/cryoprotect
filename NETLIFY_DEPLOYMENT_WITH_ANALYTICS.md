# Netlify Deployment Guide with Analytics

This guide walks you through deploying the CryoProtect application to Netlify with analytics integration.

## Prerequisites

- Node.js v18+ installed
- Netlify CLI installed (`npm install -g netlify-cli`)
- Netlify account with access to your site

## Step 1: Validate Configuration

Run the validation script to ensure everything is set up correctly:

```bash
cd /home/mushu/Projects/cryoprotect
./validate-analytics.js
```

If all checks pass, proceed to deployment. If not, fix any reported issues.

## Step 2: Deploy to Netlify

Deploy the application to Netlify:

```bash
# Make sure you're logged in
netlify login

# If this is your first deployment, initialize the site
netlify init

# Deploy a preview first to test
netlify deploy

# If preview looks good, deploy to production
netlify deploy --prod
```

Alternatively, use our deployment script:

```bash
# From the frontend directory
npm run deploy:netlify
```

## Step 3: Enable Netlify Analytics

1. Go to your Netlify dashboard: https://app.netlify.com/teams/mushu-dev/sites/cryoprotect
2. Navigate to the "Analytics" tab
3. Click "Enable Netlify Analytics"
4. Confirm the purchase ($9/month per site)

Once enabled, Netlify Analytics will start collecting data immediately. It provides:

- Page views
- Top pages
- Visitor locations
- Device and browser data
- Bandwidth usage
- 404 tracking
- No cookies or client-side JavaScript required

## Step 4: Verify Analytics Integration

1. Visit your deployed site
2. Open browser developer tools (F12)
3. Check the console for analytics messages
   - In development, you should see "[Analytics] Page view: /path" logs
   - In production, these logs are suppressed
4. Navigate between pages to verify page view tracking
5. Go to Settings > Privacy to test the analytics toggle

## Step 5: Testing Custom Events

To test custom event tracking:

1. Add the following code to any component where you want to track events:

```tsx
import { useAnalyticsContext } from '@/components/analytics/AnalyticsProvider';

function MyComponent() {
  const { trackEvent } = useAnalyticsContext();
  
  const handleAction = () => {
    // Track a custom event
    trackEvent('button_clicked', { buttonId: 'action-button' });
    
    // Rest of your handler
  };
  
  return <button onClick={handleAction}>Action</button>;
}
```

2. Test by clicking the button and checking browser console

## Step 6: Monitor Analytics Data

It takes about 24 hours for initial analytics data to appear in the Netlify dashboard. Check back regularly to see:

- Overall traffic patterns
- Popular pages
- User demographics
- 404 errors
- Bandwidth usage

## Troubleshooting

### Analytics Not Appearing

1. **Check Environment Variables**: Make sure `NEXT_PUBLIC_NETLIFY` is set to "true" in netlify.toml
2. **Verify Analytics Purchase**: Netlify Analytics is a paid add-on ($9/month)
3. **Wait 24 Hours**: Initial data may take up to 24 hours to appear
4. **Check Browser Console**: Look for errors related to analytics

### Deployment Issues

1. **Build Failures**: 
   - Check Node.js version (should be v18+)
   - Ensure all dependencies are installed (`npm install`)
   - Review build logs for specific errors

2. **API Connection Issues**:
   - Verify the Heroku backend is running
   - Check that API redirects are properly configured in netlify.toml
   - Test API connectivity with `curl https://cryoprotect.netlify.app/api/v1/health/connectivity`

## Additional Resources

- [Netlify Analytics Documentation](https://docs.netlify.com/monitor-sites/analytics/)
- [Netlify CLI Documentation](https://docs.netlify.com/cli/get-started/)
- [CryoProtect Analytics Setup Guide](./NETLIFY_ANALYTICS_SETUP.md)
- [Original Netlify Deployment Guide](./NETLIFY_DEPLOYMENT_GUIDE.md)
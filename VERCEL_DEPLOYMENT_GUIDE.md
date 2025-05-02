# Vercel Deployment Guide for CryoProtect v2

This guide explains how to deploy CryoProtect v2 on Vercel's Hobby plan while staying within the limit of 12 Serverless Functions.

## Optimization Strategy

To stay within Vercel's Hobby plan limits, we've implemented the following optimizations:

1. **Consolidated API Routes**: All API endpoints are now routed through a single serverless function, significantly reducing the function count.

2. **Static Generation**: Where possible, static content is pre-generated during build time rather than generated dynamically.

3. **Optimized Build Process**: The build process has been streamlined to focus on necessary assets.

4. **Region Optimization**: Deployment is set to a single region for better performance.

## Files Modified/Added

- `vercel.json`: Configuration file for Vercel deployment
- `api/index.py`: Consolidated API router that handles all API requests
- `build_vercel.py`: Build script for Vercel deployment
- `package.json`: Updated build scripts

## Deploying to Vercel

Follow these steps to deploy the application to Vercel:

1. **Install Vercel CLI**:
   ```bash
   npm install -g vercel
   ```

2. **Login to Vercel**:
   ```bash
   vercel login
   ```

3. **Deploy to Vercel**:
   ```bash
   vercel
   ```

4. **Production Deployment**:
   ```bash
   vercel --prod
   ```

## Environment Variables

Make sure to set the following environment variables in your Vercel project settings:

- `SUPABASE_URL`: Your Supabase instance URL
- `SUPABASE_KEY`: Your Supabase API key
- `SUPABASE_SERVICE_KEY`: Your Supabase service role key
- `SECRET_KEY`: Secret key for Flask sessions
- `API_VERSION`: API version (default: "1.0.0")
- `FLASK_ENV`: Environment (development, production)

## Monitoring Serverless Functions

You can monitor your serverless function usage in the Vercel dashboard under the "Functions" tab. This will show you how many functions are being deployed and their usage statistics.

## Troubleshooting

### Function Count Still Exceeds Limit

If you're still exceeding the function limit:

1. Check `vercel.json` to ensure all routes are properly configured
2. Review the deployment logs for any missed routes
3. Consider further consolidating API endpoints

### Other Common Issues

- **Build Failures**: Check the build logs for specific errors
- **API 500 Errors**: Check your environment variables
- **Missing Static Assets**: Ensure build scripts are correctly copying assets

## Further Optimizations

For even more efficiency, consider:

1. Implementing API caching with Vercel Edge Caching
2. Using Incremental Static Regeneration for semi-dynamic content
3. Implementing stale-while-revalidate patterns for data fetching
4. Optimizing image assets with Vercel Image Optimization

## Need Help?

Contact the CryoProtect support team or open an issue on the GitHub repository.
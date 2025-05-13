# Vercel Deployment Guide for CryoProtect

This guide explains how to deploy CryoProtect on Vercel's platform, including the challenges faced and solutions implemented during the deployment process. 

## May 2025 Deployment Status

The most recent deployment was completed on May 13, 2025. The application is available at:
https://cryoprotect-uis-pop2.vercel.app

## Optimization Strategy

To stay within Vercel's Hobby plan limits, we've implemented the following optimizations:

1. **Consolidated API Routes**: All API endpoints are now routed through a single serverless function, significantly reducing the function count.

2. **Static Generation**: Where possible, static content is pre-generated during build time rather than generated dynamically.

3. **Optimized Build Process**: The build process has been streamlined to focus on necessary assets.

4. **Region Optimization**: Deployment is set to a single region for better performance.

## Files Modified/Added

- `vercel.json`: Configuration file for Vercel deployment
- `api/index.py`: Consolidated API router that handles all API requests
- `api/app.py`: Main application landing page handler
- `api/requirements.txt`: Python dependencies for serverless functions
- `build_vercel.py`: Build script for Vercel deployment
- `package.json`: Updated build scripts with Vercel-specific commands
- `supabase/config.toml`: Supabase configuration for Vercel deployment

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

3. **Prepare the deployment**:
   - Run the build script to prepare static assets:
     ```bash
     python build_vercel.py
     ```
   - Ensure the `.env.production` file contains the correct environment variables

4. **Deploy to Vercel**:
   ```bash
   vercel deploy --prod --yes --archive=tgz
   ```

   Note: The `--archive=tgz` flag is necessary because the repository has more than 15,000 files, which exceeds Vercel's default limit.

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

### Deployment Challenges

During the deployment process, we encountered several challenges:

1. **File Count Limit**: Vercel has a 15,000 file limit for direct deployment. We resolved this by using the `--archive=tgz` flag to upload a compressed archive instead.

2. **Function Format**: Vercel expects a specific format for Python serverless functions. We updated our handler functions to return objects with the following structure:
   ```python
   {
       'statusCode': 200,
       'headers': {'Content-Type': 'application/json'},
       'body': json.dumps(data)
   }
   ```

3. **Configuration Conflicts**: Vercel doesn't allow mixing certain configuration options in `vercel.json`. We had to consolidate headers into the routes configuration.

4. **Cron Job Limitations**: The Hobby plan only allows daily cron jobs, so we had to adjust the schedule for our health check.

### Common Issues and Solutions

- **404 Errors**: Make sure your routes in `vercel.json` are correctly configured and point to the right handler functions.
  
- **Build Failures**: Check that your Python dependencies are correctly specified in `api/requirements.txt`.
  
- **Environment Variables**: Ensure all required environment variables are set in the Vercel dashboard.
  
- **API Errors**: If API calls fail, check the Supabase connection details and ensure the service role key has the necessary permissions.

## Further Optimizations

For even more efficiency, consider:

1. Implementing API caching with Vercel Edge Caching
2. Using Incremental Static Regeneration for semi-dynamic content
3. Implementing stale-while-revalidate patterns for data fetching
4. Optimizing image assets with Vercel Image Optimization

## Need Help?

Contact the CryoProtect support team or open an issue on the GitHub repository.
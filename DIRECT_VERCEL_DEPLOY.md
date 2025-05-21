# Direct Vercel Deployment Guide

This guide provides instructions for deploying the CryoProtect project directly using the Vercel CLI.

## Preparation

I've made the following changes to prevent file conflicts:

1. Renamed conflicting files:
   - Renamed `/api/index.py` to `/api/index.python.bak`
   - This prevents the conflict with `/api/index.js`

2. Created a comprehensive `.vercelignore` file that:
   - Excludes Python-related files and environments
   - Ignores large data directories
   - Reduces the total file count for deployment

3. Fixed incorrect serverless function paths in vercel.json

## Deployment Options

### Option 1: Clean Next.js Only Deployment (Recommended)

This approach creates a clean deployment with only the essential Next.js files:

```bash
cd /home/mushu/Projects/cryoprotect
./deploy-nextjs-only.sh
```

This script:
- Creates a temporary directory with only the necessary Next.js files
- Configures a clean vercel.json
- Deploys just the frontend application
- Avoids all file conflicts and limits

### Option 2: Deploy Directly from Root Directory

```bash
cd /home/mushu/Projects/cryoprotect
vercel deploy --archive=tgz
```

The `--archive=tgz` flag reduces the file size by compressing the upload.

### Option 3: Deploy to Production

When you're ready to go live:

```bash
cd /home/mushu/Projects/cryoprotect
vercel deploy --archive=tgz --prod
```

## Environment Variables

You'll need to set these environment variables during deployment:

- `NEXT_PUBLIC_API_URL` - Backend API URL
- `NEXT_PUBLIC_USE_MOCK_DATA` - Set to "false" in production
- `NEXT_PUBLIC_ENABLE_API_LOGGING` - Set to "true" to enable logging
- `NEXT_PUBLIC_ENVIRONMENT` - Set to "production" for live deployments
- `NEXTAUTH_URL` - Your frontend URL
- `NEXTAUTH_SECRET` - A secure random string
- `PROTECTION_BYPASS` - Your protection bypass token
- `VERCEL_ANALYTICS_ID` - Set to "true" to enable analytics
- `VERCEL_SPEED_INSIGHTS` - Set to "true" to enable speed insights

## Troubleshooting

### If You Encounter File Conflicts Again

1. Check for any additional conflicting file names:
   ```bash
   find . -name "*.js" | grep -i api
   find . -name "*.py" | grep -i api
   ```

2. Rename or move any conflicting files.

3. Update the `.vercelignore` file to exclude problematic files.

### If You Exceed the File Limit

If you still hit the file limit even with compression:

1. Consider deploying just the frontend directory:
   ```bash
   cd frontend
   vercel deploy --archive=tgz
   ```

2. Further optimize the `.vercelignore` file to exclude more files.

## Restoring Renamed Files

After deployment, you may restore the original files:

```bash
mv /home/mushu/Projects/cryoprotect/api/index.python.bak /home/mushu/Projects/cryoprotect/api/index.py
```

This ensures your development environment remains intact.
# Frontend Deployment Guide

This guide explains how to deploy the CryoProtect frontend to Vercel, including the Analytics and Speed Insights configuration.

## Deployment Options

You have two options for deploying the frontend:

### Option 1: Deploy Frontend Without API Conflicts (Recommended)

This is the preferred approach as it avoids both file limit issues and API conflicts:

```bash
cd frontend
./deploy-frontend-without-conflicts.sh
```

This script:
- Creates a temporary clean deployment directory
- Deploys only the frontend files
- Excludes conflicting API files
- Uses TGZ compression to reduce upload size
- Sets all required environment variables
- Includes Analytics and Speed Insights configuration

### Option 2: Deploy Frontend Only

This approach works but might encounter API conflicts:

```bash
cd frontend
./deploy-frontend-only.sh
```

### Option 3: Deploy Root Project

This option is not recommended due to file count limits but can be used if needed:

```bash
cd /home/mushu/Projects/cryoprotect
vercel deploy --archive=tgz
```

## Troubleshooting

### File Count Limits

If you encounter this error:
```
Error: Invalid request: `files` should NOT have more than 15000 items
```

Use the frontend-only deployment approach or further reduce the files using `.vercelignore`.

### API Conflicts

If you encounter this error:
```
Error: Two or more files have conflicting paths or names. Please make sure path segments and filenames, without their extension, are unique. The path "api/index.js" has conflicts with "api/index.py".
```

Use the `deploy-frontend-without-conflicts.sh` script which creates a clean temporary directory excluding the conflicting files.

### Vercel Configuration

The deployment is configured by:
- `frontend/vercel.json` - Configuration for the frontend deployment
- `frontend/.vercelignore` - Files to exclude from deployment
- `VERCEL_ANALYTICS_ID=true` - Enables Analytics
- `VERCEL_SPEED_INSIGHTS=true` - Enables Speed Insights

## Verifying Analytics & Speed Insights

After deployment:

1. Visit your Vercel dashboard
2. Select the CryoProtect project
3. Check the "Analytics" and "Speed Insights" sections

For more details, see the [VERCEL_ANALYTICS_GUIDE.md](./VERCEL_ANALYTICS_GUIDE.md) file.

## Environment Variables

Important environment variables set during deployment:

| Variable | Purpose |
|----------|---------|
| `NEXT_PUBLIC_API_URL` | Backend API URL |
| `NEXTAUTH_URL` | Frontend URL for authentication |
| `PROTECTION_BYPASS` | Security token for API access |
| `VERCEL_ANALYTICS_ID` | Enables Vercel Analytics |
| `VERCEL_SPEED_INSIGHTS` | Enables Speed Insights |

## Post-Deployment

After successful deployment, the script will output the frontend URL. Visit this URL to verify that your application is running correctly with analytics enabled.
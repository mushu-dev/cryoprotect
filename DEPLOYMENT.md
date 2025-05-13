# CryoProtect Deployment Guide

This guide provides instructions for deploying the CryoProtect application to Vercel, with specific focus on resolving common issues related to the Next.js frontend and Python Flask backend combination.

## Project Structure

CryoProtect consists of two main components:

1. **Backend API** - A Flask application serving as the API
2. **Frontend** - A Next.js application providing the user interface

These components can be deployed separately or together, depending on your needs.

## Deployment Options for Frontend

The frontend can be deployed in two ways:

1. **Standard Deployment** - Uses Next.js server capabilities for API routes and server-side rendering.
   - Supports full authentication and dynamic routes
   - Requires Vercel to support Next.js server functions

2. **Static Export** - A fully static build that works on any hosting platform.
   - Better compatibility with simple static hosts
   - Limited server-side functionality (auth pages work client-side only)
   - No dynamic API routes (AUTH API must be hosted elsewhere)

## Deployment Scripts

### Option 1: Separate Frontend and Backend Deployments (Recommended)

This approach deploys the frontend and backend as separate applications, which provides better scalability and isolation.

#### Backend Deployment

1. Navigate to the project root directory
2. Deploy to Vercel:
   ```bash
   vercel deploy --prod
   ```

#### Frontend Deployment

1. Navigate to the project root
2. Run our new deployment script:
   ```bash
   ./deploy-vercel.sh
   ```
   This script will offer you a choice between standard and static deployment

### Option 2: Combined Deployment Script

We've created a script that can deploy both components with a single command:

```bash
./deploy.sh
```

This script will prompt you to choose whether to deploy the frontend, backend, or both.

## Environment Variables

### Backend Environment Variables

- `FLASK_ENV` - Set to `production` for production deployment
- `FLASK_APP` - Set to `app.py`
- `API_VERSION` - The API version (e.g., `v1`)
- `SUPABASE_URL` - Supabase project URL
- `SUPABASE_KEY` - Supabase API key

### Frontend Environment Variables

- `NEXT_PUBLIC_API_URL` - URL of the backend API (e.g., `https://api.cryoprotect.app/v1`)
- `NEXTAUTH_URL` - URL of the frontend application (e.g., `https://www.cryoprotect.app`)
- `NEXTAUTH_SECRET` - Secret key for NextAuth.js session encryption

You can generate a secure NEXTAUTH_SECRET with:
```bash
openssl rand -base64 32
```

## Domains and URLs

The application is configured to work with the following domains:

- Frontend: `https://www.cryoprotect.app`
- Backend API: `https://api.cryoprotect.app`

If you're using different domains, make sure to update the appropriate environment variables and configuration files.

## Common Issues

### 1. The "app.py download" Issue

**Problem**: When accessing the website, the server sends `app.py` as a downloadable file instead of routing to the Next.js application.

**Solution**: This is fixed by adding this route in `vercel.json`:
```json
{
  "src": "^/app.py$",
  "status": 404,
  "dest": "/404.html"
}
```

Our deployment scripts automatically add this route to the configuration.

### 2. Next.js build errors with authentication

**Problem**: Server-side rendering (SSR) issues with authentication can cause build errors like "Cannot read properties of undefined (reading 'user')".

**Solutions**:
- Use proper optional chaining for session properties: `session?.user?.name`
- Implement safety checks for browser-specific code during SSR
- Use dynamic routes with client-side rendering for auth pages
- For static export, use the provided `build-static.sh` script which skips dynamic routes

### 3. API routes not working

**Problem**: Next.js API routes won't work with static exports.

**Solution**: For static export, use external APIs or serverless functions instead of relying on Next.js API routes.

## Local Development

For local development, you can use the provided script to run both frontend and backend servers:

```bash
./run-dev.sh
```

This script will:
1. Start the backend server on port 5000
2. Start the frontend development server on port 3000
3. Configure proper routing between them

## Build Scripts

We've created several scripts to help with building and deploying:

### 1. `frontend/build-prod.sh`

This script builds the frontend for production with modified settings to avoid SSR issues:
- Forces dynamic rendering for authenticated pages
- Sets proper configuration for NextAuth
- Preserves API routes functionality

### 2. `frontend/build-static.sh`

This script creates a fully static export:
- Generates HTML files for all static pages
- Skips dynamic routes for compatibility
- No server-side rendering or API routes

### 3. `deploy-vercel.sh`

A comprehensive deployment script that:
- Updates Vercel configuration to fix routing issues
- Offers a choice between standard and static deployment
- Handles proper file paths and environment setup

## Additional Resources

- [Vercel Documentation](https://vercel.com/docs)
- [Next.js Deployment Guide](https://nextjs.org/docs/deployment)
- [Next.js Static Export Guide](https://nextjs.org/docs/advanced-features/static-html-export)
- [Flask Deployment Guide](https://flask.palletsprojects.com/en/2.0.x/deploying/)
- [Frontend Deployment Improvements](./frontend/DEPLOYMENT_IMPROVEMENTS.md)
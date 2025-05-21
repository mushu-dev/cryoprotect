# Vercel Protection Bypass

This document explains how to use the protection bypass for Vercel deployments of the CryoProtect application.

## Overview

To protect development and staging environments from public access, we've implemented a protection bypass system. This ensures that only authorized users with the correct token can access protected deployments.

## How It Works

The protection system works in two parts:

1. **Server-side validation**: All API requests are checked for a valid protection bypass token
2. **Frontend integration**: The frontend automatically includes the token in all API requests

## Using the Protection Bypass

### For API Requests

All API requests to protected endpoints must include the protection bypass token. This can be done in one of two ways:

1. **HTTP Header**: Include the `x-protection-bypass` header with the token value
   ```
   curl -H "x-protection-bypass: YOUR_TOKEN" https://your-vercel-deployment.vercel.app/api/v1/endpoint
   ```

2. **Query Parameter**: Append the `bypass` query parameter to the URL
   ```
   curl https://your-vercel-deployment.vercel.app/api/v1/endpoint?bypass=YOUR_TOKEN
   ```

### For Frontend Users

The frontend application automatically handles the token for you. The token is:

- Set during deployment via environment variables
- Included in all API requests made from the frontend
- Not visible to end users in normal operation

## Configuration

### Environment Variables

- `PROTECTION_BYPASS`: The bypass token stored as a server-side environment variable
- `NEXT_PUBLIC_FRONTEND_PROTECTION_BYPASS`: The same token, available to the frontend client

### Deployment Configuration

The token is set in the Vercel deployment using our deployment script:

```bash
./frontend/deploy-to-vercel.sh
```

This script sets the bypass token in both environment variables.

## Development Workflow

1. For development and testing, use the standard token that's stored in the deployment script
2. For production, you may disable the protection bypass or use a more secure token

## Current Bypass Token

The current bypass token is: `TAt23KbtFE8dkZobJU3hpgTP4L5ja07V`

Note: This token should be kept confidential and shared only with authorized team members.

## Health Check Endpoint

The `/api/v1/health` endpoint is excluded from protection to allow monitoring systems to check the application status without needing the bypass token.

## Troubleshooting

If you encounter "Access Denied" errors when accessing the application:

1. Check that you're using the correct bypass token
2. Verify the token is being sent in the request header or URL
3. Check the browser console for any CORS or API errors

## Security Considerations

- The bypass token should be treated as a shared secret
- Rotate the token periodically for better security
- In production environments, consider using more robust authentication
# Convex Production Database Guide

This guide provides instructions for using Convex as the production database for CryoProtect. Convex offers real-time synchronization, automatic caching, TypeScript integration, and scalable database infrastructure.

## Table of Contents

1. [Overview](#overview)
2. [Configuration](#configuration)
3. [Development Workflow](#development-workflow)
4. [Production Deployment](#production-deployment)
5. [Data Migration](#data-migration)
6. [Testing](#testing)
7. [Monitoring](#monitoring)
8. [Troubleshooting](#troubleshooting)

## Overview

Convex is now configured as the primary database for CryoProtect in production. This implementation provides:

- Real-time data synchronization across clients
- Automatic caching and optimization
- TypeScript integration for type safety
- Simplified authentication flow
- Improved performance for complex queries

The application is now configured to switch between Supabase and Convex based on environment variables, allowing for a gradual transition.

## Configuration

### Environment Variables

The following environment variables control Convex usage:

- `NEXT_PUBLIC_CONVEX_URL` - URL of the Convex deployment
- `NEXT_PUBLIC_USE_CONVEX` - Set to `true` to use Convex, `false` to use Supabase

### File Structure

- `/frontend/convex.json` - Convex configuration file
- `/frontend/src/convex/client.ts` - Convex client initialization
- `/frontend/src/convex/ConvexClientProvider.tsx` - React provider for Convex
- `/frontend/src/convex/hooks.ts` - Custom React hooks for Convex
- `/frontend/src/convex/generated/` - Generated TypeScript types

## Development Workflow

### Starting Development Server with Convex

```bash
# Start Next.js with Convex enabled
npm run dev:with-convex

# Or manually set the environment variable
NEXT_PUBLIC_USE_CONVEX=true npm run dev
```

### Generating TypeScript Types

```bash
# Generate TypeScript types for Convex API
npm run convex:codegen
```

### Deploying Convex Functions

```bash
# Deploy Convex functions during development
npm run convex:deploy
```

### Checking Convex Status

Visit `/convex-status` in your browser to check the Convex connection status.

## Production Deployment

### Building for Production

```bash
# Build with Convex enabled
npm run build:with-convex

# Or use the deploy script for full deployment
./deploy-with-convex.sh
```

### Deployment Steps

1. Generate Convex TypeScript types
2. Deploy Convex functions
3. Build Next.js with Convex enabled
4. Generate static export
5. Deploy to Netlify

## Data Migration

### Migrating from Supabase to Convex

The application supports running both databases in parallel, allowing for gradual data migration:

1. Enable both databases: `NEXT_PUBLIC_USE_CONVEX=true`
2. Use the migration utilities to copy data from Supabase to Convex
3. Verify data integrity in the new database
4. Switch to Convex-only mode when ready

### Migration Scripts

Currently, migration must be performed manually. Automated migration scripts are being developed.

## Testing

### Running Tests

```bash
# Run Convex integration tests
npm run test:convex

# Run full test suite (includes Convex tests)
npm run test:all
```

### Test Files

- `/frontend/tests/e2e/convex-integration.spec.js` - Convex integration tests

## Monitoring

### Production Monitoring

1. Access Convex dashboard at [https://dashboard.convex.dev](https://dashboard.convex.dev)
2. Monitor:
   - Function executions
   - Database queries
   - Authentication attempts
   - Error rates

### Status Checks

The application includes a built-in status page at `/convex-status` to verify Convex connectivity.

## Troubleshooting

### Common Issues

1. **Connection Failures**
   - Verify Convex URL is correct in environment variables
   - Check network connectivity to Convex servers
   - Ensure authentication is properly configured

2. **Type Errors**
   - Run `npm run convex:codegen` to regenerate TypeScript types
   - Ensure you're importing from the correct paths

3. **Authentication Problems**
   - Verify that the ConvexClientProvider is correctly set up
   - Check auth token format in the provider

### Support

For additional help, refer to:
- Convex documentation: [https://docs.convex.dev](https://docs.convex.dev)
- CryoProtect internal documentation
- File an issue on the project repository
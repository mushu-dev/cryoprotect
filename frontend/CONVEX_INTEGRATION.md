# Convex Integration Guide

This guide explains how to use the Convex database integration with the CryoProtect frontend.

## Overview

The CryoProtect application now supports [Convex](https://convex.dev) as a backend database solution. Convex provides:

- Real-time data synchronization
- Document-based database structure
- Built-in authentication with Clerk
- TypeScript support with strong typing
- Serverless functions
- Auto-generated React hooks

## Setup Instructions

### 1. Install Dependencies

First, make sure all dependencies are installed:

```bash
npm run install-deps
```

### 2. Configure Environment Variables

Run the setup script to configure the necessary environment variables:

```bash
npm run convex:setup
```

This will prompt you for:
- Convex deployment URL
- Whether to enable Convex integration
- Clerk publishable key (if using Clerk for auth)

### 3. Generate TypeScript Types

Generate the TypeScript types for the Convex API:

```bash
npm run convex:codegen
```

Alternatively, you can run the initialization script which combines setup and codegen:

```bash
npm run convex:init
```

### 4. Development

To start both the Next.js development server and the Convex development environment:

```bash
npm run dev:with-convex
```

This will start the Next.js app and connect it to your Convex backend.

### 5. Building for Production

To build the application with Convex integration:

```bash
npm run build:with-convex
```

## Using Convex in Components

The integration provides React hooks for accessing Convex data and mutations.

### Example Usage

```tsx
import { useMolecules, useCreateMolecule } from '@/convex/hooks';

function MyComponent() {
  // Query data from Convex
  const molecules = useMolecules();
  
  // Get a mutation function
  const createMolecule = useCreateMolecule();
  
  // Create a new molecule
  const handleCreate = async () => {
    await createMolecule({
      name: "New Molecule",
      status: "active"
    });
  };
  
  return (
    <div>
      <button onClick={handleCreate}>Create Molecule</button>
      <ul>
        {molecules?.map(molecule => (
          <li key={molecule._id}>{molecule.name}</li>
        ))}
      </ul>
    </div>
  );
}
```

## Testing the Integration

Visit the `/convex-test` page in your application to test the Convex integration. This page allows you to:

- View molecules from the Convex database
- Create new molecules
- Verify real-time updates

## Structure

The Convex integration consists of these key files:

- `/src/convex/client.ts` - The Convex client configuration
- `/src/convex/ConvexClientProvider.tsx` - React provider for Convex
- `/src/convex/hooks.ts` - Custom React hooks for Convex data
- `/src/convex/generated/` - Auto-generated TypeScript types
- `/frontend/convex.json` - Convex configuration for the frontend

## Switching Between Backends

The integration is designed to allow switching between the original backend and Convex. Set the `NEXT_PUBLIC_USE_CONVEX` environment variable to `true` to use Convex.

## Advanced Usage

For more advanced usage, refer to the [Convex documentation](https://docs.convex.dev) for details on:

- Access control
- Data validation
- Pagination
- Indexes
- Scheduled functions
- File storage
- Deployment options
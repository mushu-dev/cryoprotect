# Component Migration Plan

This document outlines the plan for migrating key components from the main frontend to the minimal frontend, ensuring Convex integration works consistently across both.

## Components to Migrate

Based on analysis of both frontends, we've identified the following key components to migrate from the main frontend to the minimal frontend:

### 1. Core Data Access Components

- **Molecule Service**: The service that handles molecule data access from both API and Convex
  - Source: `/frontend/src/features/molecules/services/molecule-service.ts`
  - Target: `/minimal-frontend/src/services/molecule-service.ts`

- **Mixture Service**: The service that handles mixture data access from both API and Convex
  - Source: `/frontend/src/features/mixtures/services/mixture-service.ts`
  - Target: `/minimal-frontend/src/services/mixture-service.ts`

### 2. Resilience Components

- **Circuit Breaker**: Implements the circuit breaker pattern for API calls
  - Source: `/frontend/src/components/circuit-breaker/`
  - Target: `/minimal-frontend/src/components/circuit-breaker/`

- **Timeout Service**: Provides configurable timeouts for API calls
  - Source: `/frontend/src/services/timeout/`
  - Target: `/minimal-frontend/src/services/timeout/`

### 3. UI Components

- **Molecule List**: A component for displaying a list of molecules
  - Source: `/frontend/src/features/molecules/components/molecules-list.tsx.bak`
  - Target: `/minimal-frontend/src/components/molecules/molecules-list.tsx`

- **Molecule Card**: A card component for displaying molecule details
  - Source: `/frontend/src/features/molecules/components/molecule-card.tsx.bak`
  - Target: `/minimal-frontend/src/components/molecules/molecule-card.tsx`

- **Molecule Viewer**: Component for viewing molecule structures
  - Source: `/frontend/src/features/molecules/components/molecule-viewer-3d.tsx.bak`
  - Target: `/minimal-frontend/src/components/molecules/molecule-viewer-3d.tsx`

### 4. React Hooks

- **use-molecules**: Hook for fetching molecule data
  - Source: `/frontend/src/features/molecules/hooks/use-molecules.ts`
  - Target: `/minimal-frontend/src/hooks/use-molecules.ts`

- **use-mixtures**: Hook for fetching mixture data
  - Source: `/frontend/src/features/mixtures/hooks/use-mixtures.ts`
  - Target: `/minimal-frontend/src/hooks/use-mixtures.ts`

## Migration Strategy

To ensure smooth migration of these components, we'll follow these steps:

1. **Create Required Directory Structure**: Set up the necessary directories in the minimal frontend

2. **Create Shared TypeScript Interfaces**: Ensure both frontends use the same data types

3. **Migrate Services First**: Start with the data access services, as other components depend on them

4. **Migrate Resilience Components**: Add resilience patterns to the minimal frontend

5. **Migrate UI Components**: Adapt the UI components to work in the minimal frontend

6. **Migrate Hooks**: Update the hooks to work with the migrated services

7. **Update Pages**: Modify the existing pages in the minimal frontend to use the new components

## Implementation Notes

- **TypeScript Support**: Some components may need to be converted from TypeScript to JavaScript if the minimal frontend doesn't use TypeScript
- **Dependencies**: Check for additional dependencies that may need to be installed in the minimal frontend
- **Configuration**: Ensure the migrated components use the shared Convex configuration

## Post-Migration Testing

After migrating the components, we'll test:

1. Data fetching from both API and Convex
2. Resilience patterns (circuit breaker, timeouts)
3. UI component rendering
4. End-to-end functionality in the minimal frontend

This migration plan will ensure that the minimal frontend has the core functionality needed to work with both the API backend and Convex database, while maintaining consistency with the main frontend.
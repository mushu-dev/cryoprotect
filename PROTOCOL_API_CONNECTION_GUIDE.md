# Protocol API Connection Implementation Guide

This guide documents the changes made to connect the Protocol UI components to real backend API endpoints instead of mock data.

## Changes Made

1. **Updated Service Factory**
   - Added protocol service factory functions in `service-factory.ts`:
     - `getProtocolService()`: Returns appropriate protocol service (mock or API) based on environment
     - `getContextProtocolService(apiClient)`: Returns context-aware protocol service using the resilient API client

2. **Refactored Protocol Hooks**
   - Updated all protocol hooks in `use-protocols.ts` to use the service factory
   - Removed direct instantiation of service implementations
   - Used proper TypeScript casting for API-specific methods in specialized hooks

3. **Enhanced Context Protocol Hooks**
   - Updated `use-context-protocols.ts` to use the service factory
   - Improved resilience by leveraging the API context's error handling and connection tracking

## Architecture

The protocol service architecture follows these patterns:

### Service Interface Layer
- `ProtocolService` interface defines the contract for all protocol operations
- Service implementations include `MockProtocolService` for testing and `ApiProtocolService` for production

### Factory Pattern
- `getProtocolService()` provides the appropriate implementation
- Uses environment checks to determine whether to use mock or real services
- Follows consistent factory pattern used by other services

### Context-Enhanced Services
- `ContextProtocolService` wraps protocol operations with resilient API client
- Provides automatic retry, circuit breaking, offline caching, and connection status

### Hook Integration Layer
- Regular hooks (`useProtocols`, etc.) use basic service implementation
- Context-enabled hooks (`useContextProtocols`) leverage the API context for additional resilience
- Component hooks provide UI-specific data handling and state management

## Usage Guide

### Basic Protocol Operations

```tsx
// In a component:
import { useProtocols } from '@/features/protocols/hooks';

function ProtocolsList() {
  const { 
    protocols, 
    loading, 
    error, 
    totalCount,
    fetchProtocols
  } = useProtocols();
  
  // Use the protocols data in your UI
}
```

### With Resilient Connection Handling

```tsx
// In a component:
import { useContextProtocols } from '@/features/protocols/hooks';

function ResilientProtocolsList() {
  const { 
    protocols, 
    loading, 
    error, 
    connectionStatus,
    refreshConnection
  } = useContextProtocols();
  
  // Handle connection status changes
  // Show offline/degraded connection UI
  // Provide refresh option
}
```

### For Specialized Protocol Features

```tsx
// For mixture-specific protocol operations:
import { useProtocolsForMixture } from '@/features/protocols/hooks';

function MixtureProtocolsList({ mixtureId }) {
  const {
    protocols,
    mixtureName,
    loading,
    error
  } = useProtocolsForMixture(mixtureId);
  
  // Show mixture-specific protocols
}
```

## Error Handling

All hooks include standardized error handling:

1. Errors are captured and exposed through the hook's `error` state
2. Network and API errors are properly formatted with clear messages
3. Context-enabled hooks provide connection status awareness
4. Automatic retry with exponential backoff for transient issues
5. Circuit breaking for persistent API failures
6. Offline mode fallback with cached responses when available

## Testing

To test components with mock data:

1. Set environment variable: `NEXT_PUBLIC_USE_MOCK_DATA=true`
2. Or use URL parameter: `?mock=true`
3. Or manually mark API connection as failed: `markApiConnectionFailed()`

## Next Steps

1. Update all components to use context-enabled hooks where appropriate
2. Add loading and error states to all protocol UI components
3. Implement offline indicators when connection is lost
4. Add proper validation for protocol forms before submission
5. Enhance protocol search with real-time filtering
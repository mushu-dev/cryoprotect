# Application Resilience System

This module provides a comprehensive resilience system for handling API failures gracefully, ensuring the application remains functional even during network issues or service outages.

## Core Components

### 1. Circuit Breaker

The circuit breaker pattern prevents cascading failures by stopping requests to failing services:

- **Closed State**: Normal operation, all requests proceed
- **Open State**: Service is failing, all requests are rejected to prevent overload
- **Half-Open State**: Testing if service has recovered by allowing limited requests

```ts
import { getCircuitBreaker } from '@/services/circuit-breaker';

// Get or create a circuit breaker for a specific service
const apiCircuit = getCircuitBreaker('api');

// Execute a function with circuit breaker protection
try {
  const result = await apiCircuit.execute('get-user', async () => {
    return await fetchUserData(userId);
  });
  // Success
} catch (error) {
  // Handle failure
}
```

### 2. Retry with Exponential Backoff

Automatically retry failed requests with increasing delays to handle transient failures:

```ts
import { retryWithBackoff } from '@/services/retry-with-backoff';

try {
  const result = await retryWithBackoff(
    async () => fetchData(),
    {
      maxRetries: 3,
      initialDelay: 300,
      backoffFactor: 2,
      jitterFactor: 0.1,
      onRetry: (error, attempt, delay) => {
        console.log(`Retrying (${attempt}): ${error.message}, delay: ${delay}ms`);
      }
    }
  );
  // Success
} catch (error) {
  // All retries failed
}
```

### 3. Resilient API Client

A comprehensive API client that combines circuit breakers, retries, and fallbacks:

```ts
import { ResilientApiClient } from '@/services/resilient-api-client';

const api = new ResilientApiClient({
  baseURL: '/api',
  enableCache: true,
  logRequests: process.env.NODE_ENV === 'development'
});

// Use the client
try {
  const userData = await api.get('/users/123');
  // Success
} catch (error) {
  // Handle failure after all resilience mechanisms have been attempted
}
```

### 4. UI Components

React components for visualizing circuit breaker status:

```tsx
import { 
  CircuitBreakerStatus, 
  CircuitBreakerDashboard, 
  CircuitBreakerIndicator 
} from '@/components/circuit-breaker';

// Show status of a specific circuit breaker
<CircuitBreakerStatus circuitName="api" />

// Show all circuit breakers in a dashboard
<CircuitBreakerDashboard />

// Compact indicator for navigation headers
<CircuitBreakerIndicator circuitName="api" showLabel={false} />
```

### 5. Context Provider

A context provider for tracking circuit breaker status across the application:

```tsx
import { CircuitBreakerProvider } from '@/components/circuit-breaker';

// Wrap your application with the provider
<CircuitBreakerProvider initialCircuits={['api', 'database']}>
  <App />
</CircuitBreakerProvider>
```

### 6. Hooks

React hooks for using the resilient API client and circuit breakers:

```tsx
import { useResilientApi } from '@/hooks/use-resilient-api';

function YourComponent() {
  const { 
    api, 
    connectionStatus,
    isOnline,
    circuitState,
    resetCircuit
  } = useResilientApi();
  
  // Use the API client with resilience
  async function fetchData() {
    try {
      const data = await api.get('/data');
      return data;
    } catch (error) {
      // Handle error
    }
  }
  
  return (
    <div>
      <div>Connection: {connectionStatus}</div>
      <div>Circuit: {circuitState}</div>
      <button onClick={resetCircuit}>Reset Circuit</button>
    </div>
  );
}
```

## Common Usage Patterns

### Basic API Requests with Resilience

```tsx
import { useResilientApi } from '@/hooks/use-resilient-api';

function DataComponent() {
  const { api, isOnline, isCircuitOpen } = useResilientApi();
  const [data, setData] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  
  const fetchData = async () => {
    if (isCircuitOpen) {
      setError('Service is currently unavailable. Please try again later.');
      return;
    }
    
    setLoading(true);
    try {
      const result = await api.get('/data');
      setData(result);
      setError(null);
    } catch (error) {
      setError('Failed to fetch data: ' + error.message);
    } finally {
      setLoading(false);
    }
  };
  
  return (
    <div>
      {isOnline ? (
        <button onClick={fetchData} disabled={loading}>
          {loading ? 'Loading...' : 'Fetch Data'}
        </button>
      ) : (
        <p>You are currently offline.</p>
      )}
      
      {error && <div className="error">{error}</div>}
      {data && <div className="data">{/* Render data */}</div>}
    </div>
  );
}
```

### Showing Connection Status in Header

```tsx
import { ContextAwareIndicator } from '@/components/circuit-breaker';

function Header() {
  return (
    <header className="flex justify-between items-center p-4 border-b">
      <h1>Your App</h1>
      
      <div className="flex items-center space-x-2">
        <ContextAwareIndicator />
        <UserMenu />
      </div>
    </header>
  );
}
```

### Managing Resilience at App Level

```tsx
// _app.js or equivalent
import { CircuitBreakerProvider } from '@/components/circuit-breaker';

function MyApp({ Component, pageProps }) {
  return (
    <CircuitBreakerProvider initialCircuits={['api', 'database']}>
      <Layout>
        <Component {...pageProps} />
      </Layout>
    </CircuitBreakerProvider>
  );
}
```

## Configuration

### Circuit Breaker Configuration

```ts
const options = {
  failureThreshold: 5,              // Failures before opening circuit
  failureThresholdTimeWindow: 60000, // Time window for failures (ms)
  resetTimeout: 30000,              // Time before half-open (ms)
  successThreshold: 2,              // Successes needed to close circuit
  halfOpenMaxConcurrent: 1,         // Max concurrent requests in half-open
  maxQueueSize: 10,                 // Queue size for delayed execution
  enableFallbacks: true             // Enable fallbacks
};

const circuit = getCircuitBreaker('api', options);
```

### Retry Configuration

```ts
const options = {
  maxRetries: 3,                    // Maximum retry attempts
  initialDelay: 300,                // Initial delay (ms)
  maxDelay: 10000,                  // Maximum delay (ms)
  backoffFactor: 2,                 // Backoff multiplier
  jitterFactor: 0.1,                // Random jitter (0-1)
  isRetryable: (error) => true,     // Function to determine if retryable
  onRetry: (error, attempt, delay) => {}, // Called before each retry
  respectNavigatorOnline: true      // Respect browser online status
};

retryWithBackoff(operation, options);
```

### Resilient API Client Configuration

```ts
const options = {
  baseURL: '/api',                  // Base URL for requests
  timeout: 10000,                   // Request timeout (ms)
  circuitBreakerName: 'api',        // Circuit breaker to use
  circuitBreakerOptions: {},        // Circuit breaker options
  retryOptions: {},                 // Retry options
  enableFallbackCache: true,        // Enable response caching
  fallbackCacheMaxAge: 300000,      // Cache TTL (ms)
  enableOfflineDetection: true      // Enable offline detection
};

const api = createResilientApiClient(options);
```

## Demo

We provide a comprehensive demo page at `/resilience-demo` that showcases all resilience components and allows testing different failure scenarios:

- Successful/failed requests
- Circuit breaker tripping
- Online/offline simulation
- Request retry behavior
- Fallback responses

This demo is useful for understanding how the resilience system works and testing its configuration.

## Best Practices

1. **Add Circuit Breakers for Each Service**: Create separate circuit breakers for different services or API endpoints
2. **Use the Context Provider**: Wrap your app with `CircuitBreakerProvider` to track all circuit breakers
3. **Show Status in UI**: Use the indicators to show users when services are degraded
4. **Configure Retry Policies**: Adjust retry settings based on the operation (short for UI requests, longer for background tasks)
5. **Implement Fallbacks**: Provide fallback responses for critical operations
6. **Use the Resilient API Client**: Prefer this over direct fetch/axios calls for external API requests
7. **Cache Important Data**: Enable caching for frequently accessed data that doesn't change often

By following these practices, your application will be more resilient to network issues and service outages, providing a better user experience even when backend services are experiencing problems.
# Circuit Breaker UI Components

This module provides UI components for displaying and interacting with circuit breakers in the application.

## Overview

Circuit breakers are a design pattern used to prevent cascading failures and allow graceful degradation when services fail. 
This module provides React components for visualizing and managing circuit breakers in your application.

## Components

### CircuitBreakerProvider

A context provider that manages and tracks all circuit breakers in the application. Wrap your application with this provider to use the context-aware components.

```jsx
import { CircuitBreakerProvider } from '@/components/circuit-breaker';

function App() {
  return (
    <CircuitBreakerProvider 
      initialCircuits={['api', 'database']} 
      refreshInterval={1000}
    >
      {/* Your app components */}
    </CircuitBreakerProvider>
  );
}
```

### ContextAwareIndicator

A context-aware indicator that displays the aggregate status of all circuit breakers. Ideal for headers and status bars.

```jsx
import { ContextAwareIndicator } from '@/components/circuit-breaker';

function Header() {
  return (
    <header>
      <div className="flex items-center justify-between">
        <h1>My App</h1>
        <ContextAwareIndicator showLabels={false} />
      </div>
    </header>
  );
}
```

### CircuitBreakerIndicator

A simple indicator for a specific circuit breaker. Useful when you need to show the status of a single circuit breaker.

```jsx
import { CircuitBreakerIndicator } from '@/components/circuit-breaker';

function ApiStatus() {
  return (
    <div>
      <h2>API Status</h2>
      <CircuitBreakerIndicator circuitName="api" />
    </div>
  );
}
```

### CircuitBreakerStatus

A detailed card that shows the status and metrics of a specific circuit breaker.

```jsx
import { CircuitBreakerStatus } from '@/components/circuit-breaker';

function ApiMonitor() {
  return (
    <div>
      <h2>API Monitor</h2>
      <CircuitBreakerStatus circuitName="api" showDetails={true} />
    </div>
  );
}
```

### CircuitBreakerDashboard

A dashboard that displays all circuit breakers in the application.

```jsx
import { CircuitBreakerDashboard } from '@/components/circuit-breaker';

function MonitoringPage() {
  return (
    <div>
      <h1>System Monitoring</h1>
      <CircuitBreakerDashboard />
    </div>
  );
}
```

## Hooks

### useCircuitBreakers

A hook that provides access to the circuit breaker context.

```jsx
import { useCircuitBreakers } from '@/components/circuit-breaker';

function ServiceStatus() {
  const { 
    circuitBreakers, 
    hasOpenCircuits, 
    hasHalfOpenCircuits,
    resetAllCircuits 
  } = useCircuitBreakers();
  
  return (
    <div>
      <h2>Service Status</h2>
      {hasOpenCircuits ? (
        <p>Some services are unavailable</p>
      ) : hasHalfOpenCircuits ? (
        <p>Some services are degraded</p>
      ) : (
        <p>All services are operational</p>
      )}
      <button onClick={resetAllCircuits}>Reset All Circuits</button>
    </div>
  );
}
```

## Integration with Navigation

For the best user experience, integrate the `ContextAwareIndicator` into your application's navigation header:

```jsx
import { ContextAwareIndicator } from '@/components/circuit-breaker';

function NavigationHeader() {
  return (
    <header className="sticky top-0 z-40 border-b bg-background">
      <div className="container mx-auto flex h-16 items-center justify-between px-4">
        <div className="flex items-center gap-6">
          <Link href="/">
            <a className="font-bold text-xl">My App</a>
          </Link>
          <nav className="hidden md:flex items-center gap-6">
            {/* Navigation items */}
          </nav>
        </div>
        
        <div className="flex items-center gap-4">
          {/* Circuit Breaker Status */}
          <ContextAwareIndicator 
            showLabels={false} 
            onClick={() => router.push('/system-status')}
          />
          
          {/* User menu, etc. */}
        </div>
      </div>
    </header>
  );
}
```

## Circuit Breaker Demo Page

We provide a demo page that showcases all circuit breaker UI components and allows testing the circuit breaker functionality:

```jsx
// pages/circuit-breaker-demo.tsx
import { CircuitBreakerDashboard } from '@/components/circuit-breaker';

export default function CircuitBreakerDemo() {
  return (
    <div className="container mx-auto py-8 px-4">
      <h1 className="text-3xl font-bold mb-6">Circuit Breaker Status</h1>
      <CircuitBreakerDashboard />
    </div>
  );
}
```

## Implementation Steps

1. Wrap your application with `CircuitBreakerProvider` in `_app.js` or equivalent.
2. Add the `ContextAwareIndicator` to your navigation header.
3. Create a circuit breaker status page using `CircuitBreakerDashboard`.
4. Use the `useCircuitBreakers` hook in components that need to access circuit breaker status.

## Circuit Breaker States

Circuit breakers can be in one of three states:

- **CLOSED**: Normal operation, all requests are allowed through.
- **OPEN**: Failure state, all requests are rejected to prevent cascading failures.
- **HALF_OPEN**: Testing state, a limited number of requests are allowed through to check if the service has recovered.

The UI components use the following color coding:

- Green: Circuit is closed and operational
- Yellow: Circuit is half-open or has some failures
- Red: Circuit is open and service is unavailable

## Configuring Circuit Breakers

Circuit breakers can be configured with various options:

```js
import { getCircuitBreaker } from '@/services/circuit-breaker';

// Configure the API circuit breaker
const apiCircuit = getCircuitBreaker('api', {
  failureThreshold: 5,            // Number of failures before opening circuit
  failureThresholdTimeWindow: 60000, // Time window for failures (ms)
  resetTimeout: 30000,            // Time before attempting to half-open (ms)
  successThreshold: 2,            // Successful requests needed to close circuit
  halfOpenMaxConcurrent: 1,       // Max concurrent requests in half-open state
  maxQueueSize: 10,               // Max queue size for delayed execution
  enableFallbacks: true           // Whether to enable fallbacks
});
```

This is typically done during application initialization.
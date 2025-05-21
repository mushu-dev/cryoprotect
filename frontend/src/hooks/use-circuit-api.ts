/**
 * Circuit API Hook
 * 
 * This hook provides easy access to the circuit-breaker-enabled API client,
 * with state management for circuit status and fallbacks.
 */

import { useState, useEffect, useCallback } from 'react';
import { CircuitResilientApiClient, createCircuitResilientApiClient, CircuitClientOptions } from '@/services/circuit-resilient-api-client';
import { CircuitState } from '@/services/circuit-breaker';

// Circuit API hook options
interface UseCircuitApiOptions extends Partial<CircuitClientOptions> {
  /** Whether to use a shared client instance (default: true) */
  useSharedClient?: boolean;
  
  /** Circuit breaker name for shared instance */
  circuitName?: string;
}

// Shared client instances by name
const sharedClients = new Map<string, CircuitResilientApiClient>();

/**
 * Hook for working with circuit-breaker-enabled API client
 */
export function useCircuitApi(options: UseCircuitApiOptions = {}) {
  // Get or create API client
  const [apiClient] = useState<CircuitResilientApiClient>(() => {
    // Check if we should use a shared client
    if (options.useSharedClient !== false) {
      const clientName = options.circuitName || 'default';
      
      // Return existing shared client if available
      if (sharedClients.has(clientName)) {
        return sharedClients.get(clientName)!;
      }
      
      // Create a new shared client
      const client = createCircuitResilientApiClient({
        baseURL: options.baseURL || '/api',
        circuitBreakerName: clientName,
        ...options,
      });
      
      // Save as shared client
      sharedClients.set(clientName, client);
      return client;
    }
    
    // Create a new client instance
    return createCircuitResilientApiClient({
      baseURL: options.baseURL || '/api',
      ...options,
    });
  });
  
  // Track circuit state
  const [circuitState, setCircuitState] = useState<CircuitState>(
    apiClient.getCircuitBreaker().getSnapshot().state
  );
  
  // Track circuit metrics
  const [circuitMetrics, setCircuitMetrics] = useState(
    apiClient.getCircuitBreaker().getSnapshot()
  );
  
  // Listen for circuit state changes
  useEffect(() => {
    // Subscribe to circuit state changes
    const unsubscribe = apiClient.onCircuitStateChange(state => {
      setCircuitState(state);
      setCircuitMetrics(apiClient.getCircuitBreaker().getSnapshot());
    });
    
    // Update metrics on interval
    const intervalId = setInterval(() => {
      setCircuitMetrics(apiClient.getCircuitBreaker().getSnapshot());
    }, 1000);
    
    return () => {
      unsubscribe();
      clearInterval(intervalId);
    };
  }, [apiClient]);
  
  // Reset the circuit breaker
  const resetCircuit = useCallback(() => {
    apiClient.resetCircuit();
  }, [apiClient]);
  
  // Register a fallback for a specific operation
  const registerFallback = useCallback(
    <T>(operationKey: string, fallback: (error: Error) => Promise<T>) => {
      apiClient.registerFallback(operationKey, fallback);
    },
    [apiClient]
  );
  
  // Clear fallback cache
  const clearFallbackCache = useCallback(() => {
    apiClient.clearFallbackCache();
  }, [apiClient]);
  
  // Cache a response for fallback
  const cacheFallbackResponse = useCallback(
    <T>(method: string, url: string, data: T) => {
      apiClient.cacheFallbackResponse(method, url, data);
    },
    [apiClient]
  );
  
  return {
    // API client
    api: apiClient,
    
    // Circuit state
    circuitState,
    isCircuitOpen: circuitState === CircuitState.OPEN,
    isCircuitHalfOpen: circuitState === CircuitState.HALF_OPEN,
    isCircuitClosed: circuitState === CircuitState.CLOSED,
    
    // Circuit metrics
    circuitMetrics,
    
    // Actions
    resetCircuit,
    registerFallback,
    clearFallbackCache,
    cacheFallbackResponse,
    
    // Constants
    CircuitState,
  };
}

export default useCircuitApi;
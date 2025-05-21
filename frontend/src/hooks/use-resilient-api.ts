/**
 * Resilient API Hook
 * 
 * This hook provides a single interface for accessing a resilient API client
 * with circuit breaker integration and UI status components.
 */

import { useState, useEffect, useCallback } from 'react';
import { ResilientApiClient, ConnectionStatus } from '@/services/resilient-api-client';
import { getCircuitBreaker, CircuitState } from '@/services/circuit-breaker';

// Types for the hook result
interface UseResilientApiResult {
  // API client
  api: ResilientApiClient;
  
  // Connection status
  connectionStatus: ConnectionStatus;
  isOnline: boolean;
  isDegraded: boolean;
  isOffline: boolean;
  
  // Circuit breaker status
  circuitState: CircuitState;
  isCircuitOpen: boolean;
  isCircuitHalfOpen: boolean;
  isCircuitClosed: boolean;
  
  // Actions
  resetCircuit: () => void;
  clearCache: () => void;
  retry: () => void;
}

// Interface for the hook options
interface UseResilientApiOptions {
  /** Base URL for API requests */
  baseURL?: string;
  
  /** Circuit breaker name to use */
  circuitName?: string;
  
  /** Whether to automatically retry when connection is restored */
  autoRetryOnReconnect?: boolean;
  
  /** Mock data for offline fallbacks */
  mockData?: Record<string, any>;
  
  /** Whether to enable response caching */
  enableCache?: boolean;
}

/**
 * Hook for using a resilient API client with circuit breaker integration
 */
export function useResilientApi(options: UseResilientApiOptions = {}): UseResilientApiResult {
  // Default options
  const {
    baseURL = '/api',
    circuitName = 'api',
    autoRetryOnReconnect = true,
    mockData = {},
    enableCache = true,
  } = options;
  
  // Create or get the API client instance
  const [apiClient] = useState<ResilientApiClient>(() => {
    // Check if we already have a client in the window object
    if (typeof window !== 'undefined' && (window as any).__resilientApiClient) {
      return (window as any).__resilientApiClient;
    }
    
    // Create a new client
    const client = new ResilientApiClient({
      baseURL,
      enableCache,
      mockData,
      logRequests: process.env.NODE_ENV === 'development',
    });
    
    // Store for reuse
    if (typeof window !== 'undefined') {
      (window as any).__resilientApiClient = client;
    }
    
    return client;
  });
  
  // Get the circuit breaker
  const circuit = getCircuitBreaker(circuitName);
  
  // Track state
  const [connectionStatus, setConnectionStatus] = useState<ConnectionStatus>(
    apiClient.getConnectionStatus()
  );
  
  const [circuitState, setCircuitState] = useState<CircuitState>(
    circuit.getSnapshot().state
  );
  
  // Set up API client listeners
  useEffect(() => {
    // Listen for connection status changes
    const unsubscribeConnection = apiClient.onConnectionStatusChange(status => {
      setConnectionStatus(status);
      
      // Auto-retry on reconnect if enabled
      if (
        autoRetryOnReconnect && 
        status === 'connected' && 
        circuit.getSnapshot().state === CircuitState.OPEN
      ) {
        setTimeout(() => circuit.reset(), 1000);
      }
    });
    
    // Listen for circuit breaker changes
    const unsubscribeCircuit = circuit.onStateChange(state => {
      setCircuitState(state);
    });
    
    // Clean up listeners on unmount
    return () => {
      unsubscribeConnection();
      unsubscribeCircuit();
    };
  }, [apiClient, circuit, autoRetryOnReconnect]);
  
  // Online status listeners
  useEffect(() => {
    const handleOnline = () => {
      if (connectionStatus !== 'connected') {
        setConnectionStatus('connected');
        
        // Reset circuit on reconnect if auto-retry is enabled
        if (autoRetryOnReconnect && circuit.getSnapshot().state === CircuitState.OPEN) {
          setTimeout(() => circuit.reset(), 1000);
        }
      }
    };
    
    const handleOffline = () => {
      if (connectionStatus !== 'offline') {
        setConnectionStatus('offline');
        
        // Trip circuit when offline
        if (circuit.getSnapshot().state !== CircuitState.OPEN) {
          circuit.trip();
        }
      }
    };
    
    // Add window event listeners for online/offline events
    if (typeof window !== 'undefined') {
      window.addEventListener('online', handleOnline);
      window.addEventListener('offline', handleOffline);
    }
    
    return () => {
      if (typeof window !== 'undefined') {
        window.removeEventListener('online', handleOnline);
        window.removeEventListener('offline', handleOffline);
      }
    };
  }, [connectionStatus, circuit, autoRetryOnReconnect]);
  
  // Reset the circuit breaker
  const resetCircuit = useCallback(() => {
    circuit.reset();
  }, [circuit]);
  
  // Clear the API cache
  const clearCache = useCallback(() => {
    apiClient.clearCache();
  }, [apiClient]);
  
  // Retry connection
  const retry = useCallback(() => {
    if (typeof navigator !== 'undefined' && navigator.onLine) {
      resetCircuit();
    }
  }, [resetCircuit]);
  
  // Compute derived values
  const isOnline = connectionStatus === 'connected';
  const isDegraded = connectionStatus === 'degraded';
  const isOffline = connectionStatus === 'offline';
  
  const isCircuitOpen = circuitState === CircuitState.OPEN;
  const isCircuitHalfOpen = circuitState === CircuitState.HALF_OPEN;
  const isCircuitClosed = circuitState === CircuitState.CLOSED;
  
  return {
    api: apiClient,
    connectionStatus,
    isOnline,
    isDegraded,
    isOffline,
    circuitState,
    isCircuitOpen,
    isCircuitHalfOpen,
    isCircuitClosed,
    resetCircuit,
    clearCache,
    retry,
  };
}

export default useResilientApi;
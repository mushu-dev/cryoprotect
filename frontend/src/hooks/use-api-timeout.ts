/**
 * API Timeout Hook
 * 
 * This hook provides easy access to API timeout functionality, including:
 * - Tracking active requests
 * - Cancelling requests
 * - Setting timeouts for specific endpoints
 * - Monitoring request status
 */

import { useState, useEffect, useCallback } from 'react';
import { getTimeoutController, TimeoutConfig } from '@/services/timeout-controller';
import { createEnhancedApiClient, EnhancedResilientApiClient } from '@/services/enhanced-resilient-api-client';

interface RequestStatus {
  /** Total number of active requests */
  activeCount: number;
  
  /** Number of requests in progress (not completed) */
  inProgressCount: number;
  
  /** Number of requests that are nearly timing out (>80% of timeout elapsed) */
  nearTimeoutCount: number;
  
  /** Number of requests being cancelled */
  cancellingCount: number;
  
  /** Detailed information about active requests */
  requests: Array<{
    id: string;
    url: string;
    method: string;
    elapsedTime: number;
    timeout: number;
    completed: boolean;
    pendingCancellation: boolean;
    progress: number;
  }>;
  
  /** Whether any critical request is nearly timing out */
  hasCriticalRequests: boolean;
}

interface UseApiTimeoutOptions {
  /** Base URL for API requests */
  baseURL?: string;
  
  /** Initial timeout configuration */
  timeoutConfig?: Partial<TimeoutConfig>;
  
  /** Whether to create a new API client instance */
  createClient?: boolean;
  
  /** URL patterns that are considered critical (timeout warnings will be shown) */
  criticalUrls?: RegExp[];
}

/**
 * Hook for working with API timeouts
 */
export function useApiTimeout(options: UseApiTimeoutOptions = {}) {
  // Get the timeout controller
  const timeoutController = getTimeoutController(options.timeoutConfig);
  
  // State for request status
  const [requestStatus, setRequestStatus] = useState<RequestStatus>({
    activeCount: 0,
    inProgressCount: 0,
    nearTimeoutCount: 0,
    cancellingCount: 0,
    requests: [],
    hasCriticalRequests: false,
  });
  
  // Create API client if needed
  const [apiClient] = useState<EnhancedResilientApiClient>(() => {
    if (options.createClient) {
      return createEnhancedApiClient({
        baseURL: options.baseURL || '/api',
        timeoutConfig: options.timeoutConfig,
      });
    }
    
    return null as any; // Will use the singleton TimeoutController without an API client
  });
  
  // Update request status
  const updateRequestStatus = useCallback(() => {
    const requests = timeoutController.getActiveRequests();
    
    // Count various request types
    const activeCount = requests.length;
    const inProgressCount = requests.filter(r => !r.completed).length;
    const nearTimeoutCount = requests.filter(r => !r.completed && r.progress > 80).length;
    const cancellingCount = requests.filter(r => r.pendingCancellation).length;
    
    // Check for critical requests
    const criticalPatterns = options.criticalUrls || [];
    const hasCriticalRequests = requests.some(r => 
      !r.completed && 
      r.progress > 80 && 
      criticalPatterns.some(pattern => pattern.test(r.url))
    );
    
    setRequestStatus({
      activeCount,
      inProgressCount,
      nearTimeoutCount,
      cancellingCount,
      requests,
      hasCriticalRequests,
    });
  }, [options.criticalUrls, timeoutController]);
  
  // Listen for changes to active requests
  useEffect(() => {
    const unsubscribe = timeoutController.onRequestStatusChange(() => {
      updateRequestStatus();
    });
    
    // Initial update
    updateRequestStatus();
    
    return unsubscribe;
  }, [timeoutController, updateRequestStatus]);
  
  // Cancel a specific request
  const cancelRequest = useCallback((requestId: string) => {
    timeoutController.cancelRequest(requestId);
  }, [timeoutController]);
  
  // Cancel all active requests
  const cancelAllRequests = useCallback(() => {
    timeoutController.cancelAllRequests();
  }, [timeoutController]);
  
  // Set timeout for a specific endpoint
  const setEndpointTimeout = useCallback((endpoint: string, timeout: number) => {
    timeoutController.setEndpointTimeout(endpoint, timeout);
  }, [timeoutController]);
  
  // Update timeout configuration
  const updateTimeoutConfig = useCallback((config: Partial<TimeoutConfig>) => {
    timeoutController.updateConfig(config);
  }, [timeoutController]);
  
  return {
    // Request status
    requestStatus,
    
    // API client (if created)
    apiClient,
    
    // Timeout controller
    timeoutController,
    
    // Actions
    cancelRequest,
    cancelAllRequests,
    setEndpointTimeout,
    updateTimeoutConfig,
    
    // Helpers
    hasActiveRequests: requestStatus.inProgressCount > 0,
    hasNearTimeoutRequests: requestStatus.nearTimeoutCount > 0,
    hasCriticalRequests: requestStatus.hasCriticalRequests,
  };
}

export default useApiTimeout;
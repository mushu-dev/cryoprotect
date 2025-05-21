/**
 * Timeout Controller Service
 * 
 * This service provides enhanced timeout handling capabilities:
 * - Per-endpoint timeout configuration
 * - Adaptive timeouts based on endpoint performance
 * - Manual cancellation of long-running requests
 * - Request timeout tracking
 */

import { AxiosRequestConfig, CancelTokenSource } from 'axios';

// Types for timeout configuration
export interface TimeoutConfig {
  /** Base timeout in milliseconds */
  timeout: number;
  
  /** Whether to use adaptive timeouts */
  adaptive?: boolean;
  
  /** Maximum timeout in milliseconds (for adaptive timeouts) */
  maxTimeout?: number;
  
  /** Minimum timeout in milliseconds (for adaptive timeouts) */
  minTimeout?: number;
  
  /** Percentile to use for adaptive timeout calculation (default: 95) */
  adaptivePercentile?: number;
  
  /** Whether to allow user cancellation of requests */
  cancelable?: boolean;
  
  /** Custom timeout for specific endpoints */
  endpointTimeouts?: Record<string, number>;
}

// Interface for request tracking
interface RequestTracker {
  /** URL being requested */
  url: string;
  
  /** HTTP method */
  method: string;
  
  /** Start time of the request */
  startTime: number;
  
  /** Timeout for this request */
  timeout: number;
  
  /** Axios cancel token source */
  cancelToken: CancelTokenSource;
  
  /** Whether this request has completed */
  completed: boolean;
  
  /** Whether this request is pending user cancellation */
  pendingCancellation: boolean;
}

/**
 * Service for managing request timeouts and cancellations
 */
export class TimeoutController {
  // Default configuration
  private config: TimeoutConfig = {
    timeout: 15000, // 15 seconds default
    adaptive: true,
    maxTimeout: 30000, // 30 seconds max
    minTimeout: 5000, // 5 seconds min
    adaptivePercentile: 95,
    cancelable: true,
    endpointTimeouts: {},
  };
  
  // Request tracking
  private activeRequests: Map<string, RequestTracker> = new Map();
  
  // Performance history for adaptive timeouts
  private endpointPerformance: Record<string, number[]> = {};
  
  // Maximum history length for performance tracking
  private readonly MAX_HISTORY_LENGTH = 50;
  
  // Status change listeners
  private requestStatusListeners: Array<(activeRequests: Map<string, RequestTracker>) => void> = [];
  
  /**
   * Create a new timeout controller with the given config
   */
  constructor(config?: Partial<TimeoutConfig>) {
    this.config = { ...this.config, ...config };
  }
  
  /**
   * Get a timeout for a specific request
   * @param url The URL to request
   * @param method The HTTP method
   * @returns Timeout configuration with cancel token
   */
  getRequestConfig(url: string, method: string): {
    timeout: number;
    cancelToken: CancelTokenSource;
    requestId: string;
  } {
    const axios = require('axios');
    const cancelToken = axios.CancelToken.source();
    const requestId = this.generateRequestId(url, method);
    
    // Calculate the appropriate timeout for this request
    const timeout = this.calculateTimeout(url, method);
    
    // Track the request
    this.activeRequests.set(requestId, {
      url,
      method,
      startTime: Date.now(),
      timeout,
      cancelToken,
      completed: false,
      pendingCancellation: false,
    });
    
    // Notify listeners
    this.notifyStatusChange();
    
    // Set up automatic cleanup after timeout
    setTimeout(() => {
      this.checkRequestTimeout(requestId);
    }, timeout + 1000); // Add 1s buffer
    
    return {
      timeout,
      cancelToken,
      requestId,
    };
  }
  
  /**
   * Mark a request as completed
   * @param requestId The request ID
   * @param success Whether the request was successful
   */
  completeRequest(requestId: string, success: boolean = true): void {
    const request = this.activeRequests.get(requestId);
    
    if (request) {
      // Calculate request duration
      const duration = Date.now() - request.startTime;
      
      // Update performance history for successful requests
      if (success) {
        this.updatePerformanceHistory(request.url, request.method, duration);
      }
      
      // Mark as completed
      request.completed = true;
      
      // Remove from active requests after a short delay
      // (keeping it briefly allows UI to show completion status)
      setTimeout(() => {
        this.activeRequests.delete(requestId);
        this.notifyStatusChange();
      }, 2000);
      
      // Immediate status update to show completion
      this.notifyStatusChange();
    }
  }
  
  /**
   * Cancel a specific request
   * @param requestId The request ID
   */
  cancelRequest(requestId: string): void {
    const request = this.activeRequests.get(requestId);
    
    if (request && !request.completed) {
      request.pendingCancellation = true;
      this.notifyStatusChange();
      
      // Cancel the request
      request.cancelToken.cancel('Request cancelled by user');
      
      // Mark as completed after a short delay
      setTimeout(() => {
        if (request) {
          request.completed = true;
          this.activeRequests.delete(requestId);
          this.notifyStatusChange();
        }
      }, 1000);
    }
  }
  
  /**
   * Cancel all active requests
   */
  cancelAllRequests(): void {
    for (const [requestId, request] of this.activeRequests.entries()) {
      if (!request.completed) {
        this.cancelRequest(requestId);
      }
    }
  }
  
  /**
   * Get the status of all active requests
   */
  getActiveRequests(): Array<{
    id: string;
    url: string;
    method: string;
    elapsedTime: number;
    timeout: number;
    completed: boolean;
    pendingCancellation: boolean;
    progress: number; // 0-100 percentage of timeout elapsed
  }> {
    const now = Date.now();
    const result = [];
    
    for (const [id, request] of this.activeRequests.entries()) {
      const elapsedTime = now - request.startTime;
      const progress = Math.min(100, Math.round((elapsedTime / request.timeout) * 100));
      
      result.push({
        id,
        url: request.url,
        method: request.method,
        elapsedTime,
        timeout: request.timeout,
        completed: request.completed,
        pendingCancellation: request.pendingCancellation,
        progress,
      });
    }
    
    return result;
  }
  
  /**
   * Add a listener for request status changes
   */
  onRequestStatusChange(listener: (activeRequests: Map<string, RequestTracker>) => void): () => void {
    this.requestStatusListeners.push(listener);
    
    // Return a function to remove the listener
    return () => {
      const index = this.requestStatusListeners.indexOf(listener);
      if (index !== -1) {
        this.requestStatusListeners.splice(index, 1);
      }
    };
  }
  
  /**
   * Update the configuration
   */
  updateConfig(newConfig: Partial<TimeoutConfig>): void {
    this.config = { ...this.config, ...newConfig };
  }
  
  /**
   * Configure a specific endpoint timeout
   */
  setEndpointTimeout(endpoint: string, timeout: number): void {
    if (!this.config.endpointTimeouts) {
      this.config.endpointTimeouts = {};
    }
    
    this.config.endpointTimeouts[endpoint] = timeout;
  }
  
  /**
   * Get a timeout configuration object for axios
   */
  getAxiosConfig(url: string, method: string, existingConfig?: AxiosRequestConfig): AxiosRequestConfig {
    const { timeout, cancelToken, requestId } = this.getRequestConfig(url, method);
    
    return {
      ...existingConfig,
      timeout,
      cancelToken: cancelToken.token,
      metadata: {
        ...(existingConfig?.metadata || {}),
        timeoutRequestId: requestId,
      },
    };
  }
  
  /**
   * Create an axios request interceptor that adds timeout config
   */
  createRequestInterceptor() {
    return (config: AxiosRequestConfig) => {
      const url = config.url || '';
      const method = config.method || 'get';
      
      // Skip if already has a cancel token and timeout
      if (config.cancelToken && config.timeout) {
        return config;
      }
      
      const { timeout, cancelToken, requestId } = this.getRequestConfig(url, method);
      
      // Add timeout and cancel token
      config.timeout = timeout;
      config.cancelToken = cancelToken.token;
      
      // Add metadata
      config.metadata = {
        ...(config.metadata || {}),
        timeoutRequestId: requestId,
      };
      
      return config;
    };
  }
  
  /**
   * Create an axios response interceptor that tracks request completion
   */
  createResponseInterceptor() {
    return {
      onSuccess: (response: any) => {
        // Mark request as completed
        if (response.config?.metadata?.timeoutRequestId) {
          this.completeRequest(response.config.metadata.timeoutRequestId, true);
        }
        
        return response;
      },
      onError: (error: any) => {
        // Mark request as completed (with failure)
        if (error.config?.metadata?.timeoutRequestId) {
          this.completeRequest(error.config.metadata.timeoutRequestId, false);
        }
        
        throw error;
      },
    };
  }
  
  /**
   * Get the endpoint key from URL and method
   * @private
   */
  private getEndpointKey(url: string, method: string): string {
    // Extract the base endpoint (e.g., /api/users/123 -> /api/users)
    const urlPath = url.split('?')[0]; // Remove query params
    const pathParts = urlPath.split('/').filter(Boolean);
    
    // For REST-like endpoints, we don't want to include IDs in the endpoint key
    // Keep only the first 2-3 path parts, or less if there aren't that many
    const basePath = '/' + pathParts.slice(0, Math.min(3, pathParts.length)).join('/');
    
    return `${method.toLowerCase()}:${basePath}`;
  }
  
  /**
   * Calculate an appropriate timeout for a request
   * @private
   */
  private calculateTimeout(url: string, method: string): number {
    const endpointKey = this.getEndpointKey(url, method);
    
    // Check for endpoint-specific timeout
    if (this.config.endpointTimeouts?.[endpointKey]) {
      return this.config.endpointTimeouts[endpointKey];
    }
    
    // Use adaptive timeout if enabled and we have performance history
    if (this.config.adaptive && this.endpointPerformance[endpointKey]?.length >= 5) {
      return this.calculateAdaptiveTimeout(endpointKey);
    }
    
    // Default to base timeout
    return this.config.timeout;
  }
  
  /**
   * Calculate an adaptive timeout based on endpoint performance history
   * @private
   */
  private calculateAdaptiveTimeout(endpointKey: string): number {
    const history = this.endpointPerformance[endpointKey];
    
    // Sort the history to find the percentile
    const sortedTimes = [...history].sort((a, b) => a - b);
    
    // Calculate the index for the percentile
    const percentile = this.config.adaptivePercentile || 95;
    const index = Math.ceil((sortedTimes.length * percentile) / 100) - 1;
    
    // Get the time at the percentile
    const percentileTime = sortedTimes[index];
    
    // Add a buffer (50% extra time)
    const adaptiveTimeout = percentileTime * 1.5;
    
    // Ensure it's within min/max limits
    return Math.max(
      this.config.minTimeout || 5000,
      Math.min(
        this.config.maxTimeout || 30000,
        adaptiveTimeout
      )
    );
  }
  
  /**
   * Update the performance history for an endpoint
   * @private
   */
  private updatePerformanceHistory(url: string, method: string, duration: number): void {
    const endpointKey = this.getEndpointKey(url, method);
    
    if (!this.endpointPerformance[endpointKey]) {
      this.endpointPerformance[endpointKey] = [];
    }
    
    // Add the new duration
    this.endpointPerformance[endpointKey].push(duration);
    
    // Limit the history length
    if (this.endpointPerformance[endpointKey].length > this.MAX_HISTORY_LENGTH) {
      this.endpointPerformance[endpointKey].shift();
    }
  }
  
  /**
   * Generate a unique ID for a request
   * @private
   */
  private generateRequestId(url: string, method: string): string {
    return `${method.toLowerCase()}:${url}:${Date.now()}:${Math.random().toString(36).substr(2, 9)}`;
  }
  
  /**
   * Check if a request has timed out and clean up if needed
   * @private
   */
  private checkRequestTimeout(requestId: string): void {
    const request = this.activeRequests.get(requestId);
    
    if (request && !request.completed) {
      const elapsed = Date.now() - request.startTime;
      
      // If elapsed time exceeds timeout, clean up
      if (elapsed >= request.timeout) {
        // Request already timed out at the axios level, just clean up our tracking
        this.activeRequests.delete(requestId);
        this.notifyStatusChange();
      }
    }
  }
  
  /**
   * Notify all listeners of status changes
   * @private
   */
  private notifyStatusChange(): void {
    for (const listener of this.requestStatusListeners) {
      try {
        listener(this.activeRequests);
      } catch (error) {
        console.error('Error in request status change listener:', error);
      }
    }
  }
}

// Create a singleton instance
let timeoutControllerInstance: TimeoutController | null = null;

/**
 * Get the timeout controller instance
 */
export function getTimeoutController(config?: Partial<TimeoutConfig>): TimeoutController {
  if (!timeoutControllerInstance) {
    timeoutControllerInstance = new TimeoutController(config);
  } else if (config) {
    timeoutControllerInstance.updateConfig(config);
  }
  
  return timeoutControllerInstance;
}
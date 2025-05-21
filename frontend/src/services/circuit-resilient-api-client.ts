/**
 * Circuit-Enabled Resilient API Client
 * 
 * This client combines the EnhancedResilientApiClient with Circuit Breaker
 * pattern support for maximum resilience.
 */

import axios, { AxiosError, AxiosRequestConfig, AxiosResponse } from 'axios';
import { EnhancedResilientApiClient, EnhancedClientOptions } from './enhanced-resilient-api-client';
import { CircuitBreaker, CircuitBreakerOptions, CircuitState, getCircuitBreaker } from './circuit-breaker';

// Extended options interface
export interface CircuitClientOptions extends EnhancedClientOptions {
  circuitBreakerOptions?: CircuitBreakerOptions;
  enableFallbackCache?: boolean;
  fallbackCacheExpiry?: number;
  circuitBreakerName?: string;
}

/**
 * API Client with Circuit Breaker pattern support
 */
export class CircuitResilientApiClient extends EnhancedResilientApiClient {
  private circuitBreaker: CircuitBreaker;
  private fallbackCache: Map<string, { data: any; timestamp: number; expiry: number; }> = new Map();
  private enableFallbackCache: boolean;
  private fallbackCacheExpiry: number;
  private circuitBreakerListeners: Array<(state: CircuitState) => void> = [];
  
  /**
   * Create a new circuit-enabled API client
   */
  constructor(options: CircuitClientOptions) {
    super(options);
    
    // Get circuit breaker
    const circuitBreakerName = options.circuitBreakerName || 'default';
    this.circuitBreaker = getCircuitBreaker(circuitBreakerName, options.circuitBreakerOptions);
    
    // Set fallback cache options
    this.enableFallbackCache = options.enableFallbackCache !== false; // Default to true
    this.fallbackCacheExpiry = options.fallbackCacheExpiry || 5 * 60 * 1000; // 5 minutes default
    
    // Listen for circuit state changes
    this.circuitBreaker.onStateChange((state) => {
      // Notify circuit state change listeners
      this.notifyCircuitStateChange(state);
    });
  }
  
  /**
   * Register for circuit state changes
   */
  onCircuitStateChange(listener: (state: CircuitState) => void): () => void {
    this.circuitBreakerListeners.push(listener);
    
    // Immediately notify with current state
    listener(this.circuitBreaker.getSnapshot().state);
    
    // Return unsubscribe function
    return () => {
      const index = this.circuitBreakerListeners.indexOf(listener);
      if (index !== -1) {
        this.circuitBreakerListeners.splice(index, 1);
      }
    };
  }
  
  /**
   * Get the circuit breaker
   */
  getCircuitBreaker(): CircuitBreaker {
    return this.circuitBreaker;
  }
  
  /**
   * Override GET request to use circuit breaker
   */
  async get<T = any>(url: string, config?: AxiosRequestConfig): Promise<T> {
    const cacheKey = this.getCacheKey('get', url, config);
    
    return this.circuitBreaker.execute<T>(
      `get:${url}`,
      () => super.get<T>(url, config),
      () => this.getFallback<T>(cacheKey)
    );
  }
  
  /**
   * Override POST request to use circuit breaker
   */
  async post<T = any>(url: string, data?: any, config?: AxiosRequestConfig): Promise<T> {
    return this.circuitBreaker.execute<T>(
      `post:${url}`,
      () => super.post<T>(url, data, config)
    );
  }
  
  /**
   * Override PUT request to use circuit breaker
   */
  async put<T = any>(url: string, data?: any, config?: AxiosRequestConfig): Promise<T> {
    return this.circuitBreaker.execute<T>(
      `put:${url}`,
      () => super.put<T>(url, data, config)
    );
  }
  
  /**
   * Override PATCH request to use circuit breaker
   */
  async patch<T = any>(url: string, data?: any, config?: AxiosRequestConfig): Promise<T> {
    return this.circuitBreaker.execute<T>(
      `patch:${url}`,
      () => super.patch<T>(url, data, config)
    );
  }
  
  /**
   * Override DELETE request to use circuit breaker
   */
  async delete<T = any>(url: string, config?: AxiosRequestConfig): Promise<T> {
    return this.circuitBreaker.execute<T>(
      `delete:${url}`,
      () => super.delete<T>(url, config)
    );
  }
  
  /**
   * Register a fallback function for a specific operation
   */
  registerFallback<T>(operationKey: string, fallback: (error: Error) => Promise<T>): void {
    this.circuitBreaker.registerFallback(operationKey, fallback);
  }
  
  /**
   * Reset the circuit breaker
   */
  resetCircuit(): void {
    this.circuitBreaker.reset();
  }
  
  /**
   * Clear the fallback cache
   */
  clearFallbackCache(): void {
    this.fallbackCache.clear();
  }
  
  /**
   * Create a cache key for a request
   * @private
   */
  private getCacheKey(method: string, url: string, config?: AxiosRequestConfig): string {
    const params = config?.params ? JSON.stringify(config.params) : '';
    return `${method}:${url}:${params}`;
  }
  
  /**
   * Get a fallback response from cache
   * @private
   */
  private async getFallback<T>(cacheKey: string): Promise<T> {
    // Check if we have a cached response
    if (this.enableFallbackCache) {
      const cached = this.fallbackCache.get(cacheKey);
      
      // If cache exists and is not expired
      if (cached && Date.now() - cached.timestamp < cached.expiry) {
        return cached.data as T;
      }
    }
    
    // No cached fallback, throw error
    throw new Error(`Service unavailable (circuit open) and no fallback data available for ${cacheKey}`);
  }
  
  /**
   * Cache a successful response for fallback
   * @private
   */
  cacheFallbackResponse<T>(method: string, url: string, data: T, config?: AxiosRequestConfig): void {
    if (this.enableFallbackCache) {
      const cacheKey = this.getCacheKey(method, url, config);
      
      this.fallbackCache.set(cacheKey, {
        data,
        timestamp: Date.now(),
        expiry: this.fallbackCacheExpiry,
      });
    }
  }
  
  /**
   * Notify circuit state change listeners
   * @private
   */
  private notifyCircuitStateChange(state: CircuitState): void {
    for (const listener of this.circuitBreakerListeners) {
      try {
        listener(state);
      } catch (error) {
        console.error('Error in circuit state change listener:', error);
      }
    }
  }
}

/**
 * Create a circuit-enabled API client
 */
export function createCircuitResilientApiClient(options: CircuitClientOptions): CircuitResilientApiClient {
  return new CircuitResilientApiClient(options);
}
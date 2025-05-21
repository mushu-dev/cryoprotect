/**
 * Resilient API Client
 * 
 * A robust API client with built-in resilience mechanisms including:
 * - Retry with exponential backoff
 * - Circuit breaker pattern
 * - Timeout management
 * - Service health tracking
 * - Offline mode with local storage caching
 * - Comprehensive error handling
 */

import axios, { AxiosError, AxiosInstance, AxiosRequestConfig, AxiosResponse } from 'axios';

// Constants for resilience configurations
const DEFAULT_TIMEOUT = 15000; // 15 seconds default timeout
const MAX_RETRIES = 3; // Maximum number of retries
const INITIAL_BACKOFF = 300; // Initial backoff in milliseconds
const MAX_BACKOFF = 10000; // Maximum backoff in milliseconds
const CIRCUIT_BREAKER_THRESHOLD = 5; // Number of failures before opening circuit
const CIRCUIT_BREAKER_RESET_TIMEOUT = 30000; // Time before attempting to close circuit (30 seconds)
const CACHE_EXPIRY = 24 * 60 * 60 * 1000; // 24 hours cache expiry

// Error types that are considered retryable
const RETRYABLE_ERRORS = [
  'ECONNABORTED', // Timeout
  'ETIMEDOUT', // Connection timeout
  'ECONNRESET', // Connection reset
  'ECONNREFUSED', // Connection refused
  'ENETUNREACH', // Network unreachable
  'EHOSTUNREACH', // Host unreachable
  'NETWORK_ERROR', // General network error
];

// Status codes that are considered retryable
const RETRYABLE_STATUS_CODES = [
  408, // Request Timeout
  429, // Too Many Requests
  500, // Internal Server Error
  502, // Bad Gateway
  503, // Service Unavailable
  504, // Gateway Timeout
];

// Interface for service health tracking
interface ServiceHealth {
  failures: number;
  successes: number;
  lastFailure: number | null;
  circuitOpen: boolean;
  circuitOpenTime: number | null;
}

// Interface for cache entry
interface CacheEntry<T> {
  data: T;
  timestamp: number;
  expiry: number;
}

// Options for the resilient client
export interface ResilientClientOptions {
  baseURL: string;
  timeout?: number;
  authToken?: string;
  retries?: number;
  enableCache?: boolean;
  cachePrefix?: string;
  circuitBreakerThreshold?: number;
  circuitBreakerResetTimeout?: number;
  defaultHeaders?: Record<string, string>;
  logRequests?: boolean;
  mockMode?: boolean;
  mockData?: Record<string, any>;
}

// Connection status types
export type ConnectionStatus = 'connected' | 'degraded' | 'offline';

/**
 * Resilient API Client class that handles network errors, retries, and circuit breaking
 */
export class ResilientApiClient {
  private apiClient: AxiosInstance;
  private serviceHealth: Record<string, ServiceHealth>;
  private cacheEnabled: boolean;
  private cachePrefix: string;
  private mockMode: boolean;
  private mockData: Record<string, any>;
  private maxRetries: number;
  private circuitBreakerThreshold: number;
  private circuitBreakerResetTimeout: number;
  private connectionListeners: Array<(status: ConnectionStatus) => void>;
  private currentConnectionStatus: ConnectionStatus;
  
  /**
   * Creates a new instance of the Resilient API Client
   * @param options Client configuration options
   */
  constructor(options: ResilientClientOptions) {
    // Initialize service health tracker
    this.serviceHealth = {};
    
    // Set configuration options
    this.cacheEnabled = options.enableCache ?? true;
    this.cachePrefix = options.cachePrefix ?? 'resilient_api_cache_';
    this.mockMode = options.mockMode ?? false;
    this.mockData = options.mockData ?? {};
    this.maxRetries = options.retries ?? MAX_RETRIES;
    this.circuitBreakerThreshold = options.circuitBreakerThreshold ?? CIRCUIT_BREAKER_THRESHOLD;
    this.circuitBreakerResetTimeout = options.circuitBreakerResetTimeout ?? CIRCUIT_BREAKER_RESET_TIMEOUT;
    this.connectionListeners = [];
    this.currentConnectionStatus = 'connected';
    
    // Create axios instance with default config
    this.apiClient = axios.create({
      baseURL: options.baseURL,
      timeout: options.timeout ?? DEFAULT_TIMEOUT,
      headers: {
        'Content-Type': 'application/json',
        ...(options.authToken ? { 'Authorization': `Bearer ${options.authToken}` } : {}),
        ...(options.defaultHeaders ?? {})
      },
    });
    
    // Add request interceptor
    this.apiClient.interceptors.request.use(
      (config) => {
        // Generate unique service key for this request
        const serviceKey = this.getServiceKey(config);
        
        // Create service health entry if it doesn't exist
        if (!this.serviceHealth[serviceKey]) {
          this.serviceHealth[serviceKey] = {
            failures: 0,
            successes: 0,
            lastFailure: null,
            circuitOpen: false,
            circuitOpenTime: null,
          };
        }
        
        // Check circuit breaker status
        if (this.isCircuitOpen(serviceKey)) {
          const shouldAttemptReset = this.shouldAttemptCircuitReset(serviceKey);
          
          if (!shouldAttemptReset) {
            // Circuit is open and shouldn't reset yet
            return Promise.reject(new Error(`Circuit breaker is open for ${serviceKey}`));
          }
          
          // This request will be a test request to see if the service is healthy again
          if (options.logRequests) {
            console.log(`Attempting to reset circuit breaker for ${serviceKey} with test request`);
          }
        }
        
        // Check if we should use a cached response
        if (this.cacheEnabled && config.method?.toLowerCase() === 'get') {
          config.metadata = {
            ...(config.metadata || {}),
            cacheKey: this.getCacheKey(config),
          };
        }
        
        // Add request logging if needed
        if (options.logRequests) {
          console.log(`API Request: ${config.method?.toUpperCase()} ${config.url}`);
        }
        
        return config;
      },
      (error) => Promise.reject(error)
    );
    
    // Add response interceptor
    this.apiClient.interceptors.response.use(
      (response) => {
        // Track successful request
        const serviceKey = this.getServiceKey(response.config);
        this.trackSuccess(serviceKey);
        
        // Cache the response if applicable
        if (
          this.cacheEnabled && 
          response.config.method?.toLowerCase() === 'get' && 
          response.config.metadata?.cacheKey
        ) {
          this.setCacheResponse(
            response.config.metadata.cacheKey,
            response.data
          );
        }
        
        // Update connection status
        if (this.currentConnectionStatus !== 'connected') {
          this.updateConnectionStatus('connected');
        }
        
        // Add response logging if needed
        if (options.logRequests) {
          console.log(`API Response: ${response.status} ${response.config.url}`);
        }
        
        return response;
      },
      async (error: AxiosError) => {
        // Extract the original request config
        const config = error.config as AxiosRequestConfig & { 
          _retryCount?: number;
          metadata?: {
            cacheKey?: string;
            serviceKey?: string;
          };
        };
        
        // Track failure if we have a service key
        if (config && config.url) {
          const serviceKey = this.getServiceKey(config);
          this.trackFailure(serviceKey);
        }
        
        // Try to get cached response if it's a GET request and cache is enabled
        if (
          this.cacheEnabled && 
          config && 
          config.method?.toLowerCase() === 'get' && 
          config.metadata?.cacheKey
        ) {
          const cachedResponse = this.getCachedResponse(config.metadata.cacheKey);
          if (cachedResponse) {
            if (options.logRequests) {
              console.log(`Using cached response for ${config.url}`);
            }
            
            // Create a mock successful response with the cached data
            return Promise.resolve({
              data: cachedResponse,
              status: 200,
              statusText: 'OK (Cached)',
              headers: {},
              config,
              request: null,
              cached: true,
            });
          }
        }
        
        // Check if we should retry this request
        if (config && this.isRetryableError(error) && this.canRetry(config)) {
          // Increment retry count
          const retryCount = config._retryCount || 0;
          config._retryCount = retryCount + 1;
          
          // Calculate backoff time with jitter
          const backoffTime = Math.min(
            INITIAL_BACKOFF * Math.pow(2, retryCount) + Math.random() * 100,
            MAX_BACKOFF
          );
          
          if (options.logRequests) {
            console.log(`Retrying request (${config._retryCount}/${this.maxRetries}) after ${backoffTime}ms`);
          }
          
          // Return a promise that resolves after backoff time
          return new Promise(resolve => {
            setTimeout(() => {
              resolve(this.apiClient(config));
            }, backoffTime);
          });
        }
        
        // Check if we should update connection status
        if (this.isNetworkError(error)) {
          this.updateConnectionStatus('offline');
        } else if (this.isServiceDegradation(error)) {
          this.updateConnectionStatus('degraded');
        }
        
        // Check if we should fall back to mock mode
        if (this.shouldEnableMockMode(error)) {
          if (options.logRequests) {
            console.log('Enabling mock mode due to persistent API failures');
          }
          this.enableMockMode();
          
          // If we're in mock mode, handle the request with mock data
          if (config && config.url) {
            const mockResponse = this.getMockResponse(config.url, config.method || 'get');
            if (mockResponse) {
              if (options.logRequests) {
                console.log('Returning mock data for failed request');
              }
              return Promise.resolve({
                data: mockResponse,
                status: 200,
                statusText: 'OK (Mock)',
                headers: {},
                config,
                request: null,
                mock: true,
              });
            }
          }
        }
        
        // If we can't retry, use cache, or mock, reject with error
        return Promise.reject(error);
      }
    );
  }

  /**
   * GET request with resilience
   * @param url URL to request
   * @param config Optional axios config
   */
  async get<T = any>(url: string, config?: AxiosRequestConfig): Promise<T> {
    try {
      const response = await this.apiClient.get<T>(url, config);
      return response.data;
    } catch (error) {
      if (axios.isAxiosError(error)) {
        // Format axios errors nicely
        const status = error.response?.status;
        const errorData = error.response?.data as any;
        const errorMessage = errorData?.message || errorData?.error || error.message;
        
        throw new Error(
          `API Error (GET ${url}): ${status ? `${status} ` : ''}${errorMessage}`
        );
      }
      
      // Re-throw other errors
      throw error;
    }
  }

  /**
   * POST request with resilience
   * @param url URL to request
   * @param data Data to send
   * @param config Optional axios config
   */
  async post<T = any>(url: string, data?: any, config?: AxiosRequestConfig): Promise<T> {
    try {
      const response = await this.apiClient.post<T>(url, data, config);
      return response.data;
    } catch (error) {
      if (axios.isAxiosError(error)) {
        // Format axios errors nicely
        const status = error.response?.status;
        const errorData = error.response?.data as any;
        const errorMessage = errorData?.message || errorData?.error || error.message;
        
        throw new Error(
          `API Error (POST ${url}): ${status ? `${status} ` : ''}${errorMessage}`
        );
      }
      
      // Re-throw other errors
      throw error;
    }
  }

  /**
   * PUT request with resilience
   * @param url URL to request
   * @param data Data to send
   * @param config Optional axios config
   */
  async put<T = any>(url: string, data?: any, config?: AxiosRequestConfig): Promise<T> {
    try {
      const response = await this.apiClient.put<T>(url, data, config);
      return response.data;
    } catch (error) {
      if (axios.isAxiosError(error)) {
        // Format axios errors nicely
        const status = error.response?.status;
        const errorData = error.response?.data as any;
        const errorMessage = errorData?.message || errorData?.error || error.message;
        
        throw new Error(
          `API Error (PUT ${url}): ${status ? `${status} ` : ''}${errorMessage}`
        );
      }
      
      // Re-throw other errors
      throw error;
    }
  }

  /**
   * PATCH request with resilience
   * @param url URL to request
   * @param data Data to send
   * @param config Optional axios config
   */
  async patch<T = any>(url: string, data?: any, config?: AxiosRequestConfig): Promise<T> {
    try {
      const response = await this.apiClient.patch<T>(url, data, config);
      return response.data;
    } catch (error) {
      if (axios.isAxiosError(error)) {
        // Format axios errors nicely
        const status = error.response?.status;
        const errorData = error.response?.data as any;
        const errorMessage = errorData?.message || errorData?.error || error.message;
        
        throw new Error(
          `API Error (PATCH ${url}): ${status ? `${status} ` : ''}${errorMessage}`
        );
      }
      
      // Re-throw other errors
      throw error;
    }
  }

  /**
   * DELETE request with resilience
   * @param url URL to request
   * @param config Optional axios config
   */
  async delete<T = any>(url: string, config?: AxiosRequestConfig): Promise<T> {
    try {
      const response = await this.apiClient.delete<T>(url, config);
      return response.data;
    } catch (error) {
      if (axios.isAxiosError(error)) {
        // Format axios errors nicely
        const status = error.response?.status;
        const errorData = error.response?.data as any;
        const errorMessage = errorData?.message || errorData?.error || error.message;
        
        throw new Error(
          `API Error (DELETE ${url}): ${status ? `${status} ` : ''}${errorMessage}`
        );
      }
      
      // Re-throw other errors
      throw error;
    }
  }

  /**
   * Register a connection status listener
   * @param listener Function to call when connection status changes
   * @returns Function to unregister the listener
   */
  onConnectionStatusChange(listener: (status: ConnectionStatus) => void): () => void {
    this.connectionListeners.push(listener);
    
    // Call the listener with the current status
    listener(this.currentConnectionStatus);
    
    // Return a function to unregister the listener
    return () => {
      const index = this.connectionListeners.indexOf(listener);
      if (index !== -1) {
        this.connectionListeners.splice(index, 1);
      }
    };
  }

  /**
   * Get the current connection status
   */
  getConnectionStatus(): ConnectionStatus {
    return this.currentConnectionStatus;
  }

  /**
   * Updates the connection status and notifies listeners
   * @param status New connection status
   */
  private updateConnectionStatus(status: ConnectionStatus): void {
    if (this.currentConnectionStatus !== status) {
      this.currentConnectionStatus = status;
      
      // Notify all listeners
      this.connectionListeners.forEach(listener => {
        try {
          listener(status);
        } catch (error) {
          console.error('Error in connection status listener:', error);
        }
      });
      
      // Store status in localStorage for persistence across page loads
      if (typeof window !== 'undefined') {
        localStorage.setItem('api_connection_status', status);
      }
    }
  }

  /**
   * Generates a unique key for a service based on the request
   * @param config Request configuration
   */
  private getServiceKey(config: AxiosRequestConfig): string {
    // Create a service key based on the URL
    const baseUrl = new URL(config.url || '/', config.baseURL).origin;
    const pathname = new URL(config.url || '/', config.baseURL).pathname.split('/')[1];
    
    return `${baseUrl}/${pathname || 'root'}`;
  }

  /**
   * Checks if the circuit breaker is currently open for a service
   * @param serviceKey Service identifier
   */
  private isCircuitOpen(serviceKey: string): boolean {
    return this.serviceHealth[serviceKey]?.circuitOpen || false;
  }
  
  /**
   * Determines if we should attempt to reset the circuit
   * @param serviceKey Service identifier
   */
  private shouldAttemptCircuitReset(serviceKey: string): boolean {
    const health = this.serviceHealth[serviceKey];
    if (!health || !health.circuitOpen || !health.circuitOpenTime) {
      return false;
    }
    
    const now = Date.now();
    return (now - health.circuitOpenTime) >= this.circuitBreakerResetTimeout;
  }
  
  /**
   * Tracks a successful API request
   * @param serviceKey Service identifier
   */
  private trackSuccess(serviceKey: string): void {
    if (!this.serviceHealth[serviceKey]) {
      this.serviceHealth[serviceKey] = {
        failures: 0,
        successes: 1,
        lastFailure: null,
        circuitOpen: false,
        circuitOpenTime: null,
      };
      return;
    }
    
    this.serviceHealth[serviceKey].successes += 1;
    
    // If the circuit is half-open and we get a success, close it
    if (
      this.serviceHealth[serviceKey].circuitOpen && 
      this.shouldAttemptCircuitReset(serviceKey)
    ) {
      this.serviceHealth[serviceKey].circuitOpen = false;
      this.serviceHealth[serviceKey].circuitOpenTime = null;
      this.serviceHealth[serviceKey].failures = 0;
    }
  }
  
  /**
   * Tracks a failed API request
   * @param serviceKey Service identifier
   */
  private trackFailure(serviceKey: string): void {
    if (!this.serviceHealth[serviceKey]) {
      this.serviceHealth[serviceKey] = {
        failures: 1,
        successes: 0,
        lastFailure: Date.now(),
        circuitOpen: false,
        circuitOpenTime: null,
      };
      return;
    }
    
    this.serviceHealth[serviceKey].failures += 1;
    this.serviceHealth[serviceKey].lastFailure = Date.now();
    
    // Check if we should open the circuit
    if (
      !this.serviceHealth[serviceKey].circuitOpen && 
      this.serviceHealth[serviceKey].failures >= this.circuitBreakerThreshold
    ) {
      this.serviceHealth[serviceKey].circuitOpen = true;
      this.serviceHealth[serviceKey].circuitOpenTime = Date.now();
    }
  }
  
  /**
   * Determines if an error is retryable
   * @param error Axios error
   */
  private isRetryableError(error: AxiosError): boolean {
    // Check if it's a network error or timeout
    if (error.code && RETRYABLE_ERRORS.includes(error.code)) {
      return true;
    }
    
    // Check if it's a retryable status code
    if (error.response?.status && RETRYABLE_STATUS_CODES.includes(error.response.status)) {
      return true;
    }
    
    // Check for specific error conditions in the response
    const errorMessage = error.response?.data?.error || '';
    return errorMessage.includes('timeout') || 
           errorMessage.includes('overloaded') || 
           errorMessage.includes('rate limit');
  }

  /**
   * Determines if an error is a network error
   * @param error Axios error
   */
  private isNetworkError(error: AxiosError): boolean {
    return !error.response && !!error.request;
  }

  /**
   * Determines if an error indicates service degradation
   * @param error Axios error
   */
  private isServiceDegradation(error: AxiosError): boolean {
    // Check if it's a server error (5xx)
    if (error.response?.status && error.response.status >= 500 && error.response.status < 600) {
      return true;
    }
    
    // Check if we have multiple services with open circuits
    const openCircuits = Object.values(this.serviceHealth).filter(health => health.circuitOpen);
    return openCircuits.length > 1;
  }
  
  /**
   * Determines if a request can be retried based on retry count
   * @param config Request configuration
   */
  private canRetry(config?: AxiosRequestConfig & { _retryCount?: number }): boolean {
    if (!config) return false;
    
    const retryCount = config._retryCount || 0;
    return retryCount < this.maxRetries;
  }
  
  /**
   * Determines if we should enable mock mode due to persistent API failures
   * @param error Axios error
   */
  private shouldEnableMockMode(error: AxiosError): boolean {
    // Enable mock mode if:
    // 1. Multiple services have open circuits OR
    // 2. The API seems completely inaccessible
    
    const openCircuits = Object.values(this.serviceHealth).filter(health => health.circuitOpen);
    const apiInaccessible = error.code === 'ECONNREFUSED' || 
                           error.code === 'ENOTFOUND' ||
                           (error.message && error.message.includes('Network Error'));
    
    return openCircuits.length >= 2 || apiInaccessible;
  }
  
  /**
   * Enables mock mode for offline operation
   */
  private enableMockMode(): void {
    if (!this.mockMode) {
      console.log('Enabling mock mode for offline operation');
      this.mockMode = true;
      
      // Notify offline mode
      this.updateConnectionStatus('offline');
      
      // Store offline flag in localStorage
      if (typeof window !== 'undefined') {
        localStorage.setItem('api_offline_mode', 'true');
      }
    }
  }
  
  /**
   * Gets a mock response for a given URL and method
   * @param url URL
   * @param method HTTP method
   */
  private getMockResponse(url: string, method: string): any {
    // First, try to get a mock by exact URL and method
    const exactKey = `${method.toLowerCase()}:${url}`;
    if (this.mockData[exactKey]) {
      return this.mockData[exactKey];
    }
    
    // If no exact match, try pattern matching
    for (const key in this.mockData) {
      if (key.startsWith(`${method.toLowerCase()}:`)) {
        const pattern = key.substring(method.length + 1);
        const regex = new RegExp(pattern);
        if (regex.test(url)) {
          return this.mockData[key];
        }
      }
    }
    
    // No mock available
    return null;
  }
  
  /**
   * Generates a cache key for a request
   * @param config Request configuration
   */
  private getCacheKey(config: AxiosRequestConfig): string {
    // Create a unique key based on URL and params
    const url = config.url || '';
    const params = config.params ? JSON.stringify(config.params) : '';
    
    return `${this.cachePrefix}${url}:${params}`;
  }
  
  /**
   * Gets a cached response if available and not expired
   * @param cacheKey Cache key
   */
  private getCachedResponse(cacheKey: string): any {
    if (!this.cacheEnabled || typeof window === 'undefined') {
      return null;
    }
    
    try {
      const cachedItem = localStorage.getItem(cacheKey);
      if (!cachedItem) {
        return null;
      }
      
      const cacheEntry: CacheEntry<any> = JSON.parse(cachedItem);
      
      // Check if cache is expired
      if (Date.now() > cacheEntry.timestamp + cacheEntry.expiry) {
        // Remove expired cache
        localStorage.removeItem(cacheKey);
        return null;
      }
      
      return cacheEntry.data;
    } catch (error) {
      console.error('Error retrieving from cache:', error);
      return null;
    }
  }
  
  /**
   * Sets a response in the cache
   * @param cacheKey Cache key
   * @param data Response data
   * @param expiry Cache expiry in milliseconds (defaults to 24 hours)
   */
  private setCacheResponse(cacheKey: string, data: any, expiry: number = CACHE_EXPIRY): void {
    if (!this.cacheEnabled || typeof window === 'undefined') {
      return;
    }
    
    try {
      const cacheEntry: CacheEntry<any> = {
        data,
        timestamp: Date.now(),
        expiry,
      };
      
      localStorage.setItem(cacheKey, JSON.stringify(cacheEntry));
    } catch (error) {
      console.error('Error saving to cache:', error);
    }
  }
  
  /**
   * Clears all cached responses
   */
  clearCache(): void {
    if (typeof window === 'undefined') {
      return;
    }
    
    // Remove all items with this cache prefix
    for (let i = 0; i < localStorage.length; i++) {
      const key = localStorage.key(i);
      if (key && key.startsWith(this.cachePrefix)) {
        localStorage.removeItem(key);
      }
    }
  }
  
  /**
   * Adds mock data for testing or offline fallbacks
   * @param mockData Record of URL patterns to mock responses
   */
  addMockData(mockData: Record<string, any>): void {
    this.mockData = {
      ...this.mockData,
      ...mockData,
    };
  }
  
  /**
   * Checks the overall health of the API
   */
  getApiHealth(): {
    healthy: boolean;
    degraded: boolean;
    services: Record<string, {
      healthy: boolean;
      failures: number;
      successes: number;
      circuitOpen: boolean;
    }>;
  } {
    const services: Record<string, {
      healthy: boolean;
      failures: number;
      successes: number;
      circuitOpen: boolean;
    }> = {};
    
    let openCircuitCount = 0;
    
    // Compile service health information
    for (const [key, health] of Object.entries(this.serviceHealth)) {
      services[key] = {
        healthy: !health.circuitOpen,
        failures: health.failures,
        successes: health.successes,
        circuitOpen: health.circuitOpen,
      };
      
      if (health.circuitOpen) {
        openCircuitCount++;
      }
    }
    
    // Overall health assessment
    const healthy = openCircuitCount === 0;
    const degraded = openCircuitCount > 0;
    
    return {
      healthy,
      degraded,
      services,
    };
  }
}
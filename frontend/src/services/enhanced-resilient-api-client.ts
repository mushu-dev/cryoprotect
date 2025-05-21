/**
 * Enhanced Resilient API Client
 * 
 * A robust API client that extends the ResilientApiClient with
 * additional resiliency features:
 * - Enhanced timeout handling with adaptive timeouts
 * - Request cancellation support
 * - Live request tracking
 * - Resource prioritization
 */

import axios, { AxiosError, AxiosInstance, AxiosRequestConfig, AxiosResponse } from 'axios';
import { ResilientApiClient, ResilientClientOptions, ConnectionStatus } from './resilient-api-client';
import { TimeoutController, TimeoutConfig, getTimeoutController } from './timeout-controller';

// Enhanced options interface
export interface EnhancedClientOptions extends ResilientClientOptions {
  timeoutConfig?: Partial<TimeoutConfig>;
  priorityHeader?: string;
  priorityHeaderValue?: string;
}

/**
 * Enhanced Resilient API Client with advanced timeout handling
 */
export class EnhancedResilientApiClient extends ResilientApiClient {
  private timeoutController: TimeoutController;
  private priorityHeader?: string;
  private priorityHeaderValue?: string;
  
  /**
   * Create a new instance of the Enhanced Resilient API Client
   */
  constructor(options: EnhancedClientOptions) {
    // Initialize the base client
    super(options);
    
    // Initialize timeout controller
    this.timeoutController = getTimeoutController(options.timeoutConfig);
    
    // Set priority headers if provided
    this.priorityHeader = options.priorityHeader;
    this.priorityHeaderValue = options.priorityHeaderValue;
    
    // Add timeout controller interceptors
    this.addTimeoutInterceptors();
  }
  
  /**
   * Add timeout controller interceptors to the API client
   * @private
   */
  private addTimeoutInterceptors(): void {
    // Get the axios instance from the parent class
    const axiosInstance = this.getAxiosInstance();
    
    if (axiosInstance) {
      // Add request interceptor for timeout handling
      axiosInstance.interceptors.request.use(
        this.timeoutController.createRequestInterceptor(),
        (error) => Promise.reject(error)
      );
      
      // Add response interceptor for tracking request completion
      const responseInterceptor = this.timeoutController.createResponseInterceptor();
      axiosInstance.interceptors.response.use(
        responseInterceptor.onSuccess,
        responseInterceptor.onError
      );
    }
  }
  
  /**
   * Exposed method to get the axios instance (protected in parent class)
   */
  getAxiosInstance(): AxiosInstance | null {
    // This is a hack to access the protected apiClient from the parent class
    // In a real implementation, you would modify the parent class to expose this
    return (this as any).apiClient || null;
  }
  
  /**
   * Set timeout for a specific endpoint
   */
  setEndpointTimeout(endpoint: string, timeout: number): void {
    this.timeoutController.setEndpointTimeout(endpoint, timeout);
  }
  
  /**
   * Update timeout configuration
   */
  updateTimeoutConfig(config: Partial<TimeoutConfig>): void {
    this.timeoutController.updateConfig(config);
  }
  
  /**
   * Get active requests status
   */
  getActiveRequests() {
    return this.timeoutController.getActiveRequests();
  }
  
  /**
   * Cancel a specific request
   */
  cancelRequest(requestId: string): void {
    this.timeoutController.cancelRequest(requestId);
  }
  
  /**
   * Cancel all active requests
   */
  cancelAllRequests(): void {
    this.timeoutController.cancelAllRequests();
  }
  
  /**
   * Register a listener for request status changes
   */
  onRequestStatusChange(
    listener: (requests: ReturnType<TimeoutController['getActiveRequests']>) => void
  ): () => void {
    return this.timeoutController.onRequestStatusChange((requestMap) => {
      // Convert the request map to the format expected by the listener
      listener(this.timeoutController.getActiveRequests());
    });
  }
  
  /**
   * Override GET request to add priority headers if needed
   */
  async get<T = any>(url: string, config?: AxiosRequestConfig): Promise<T> {
    const enhancedConfig = this.addPriorityHeaders(config);
    return super.get<T>(url, enhancedConfig);
  }
  
  /**
   * Override POST request to add priority headers if needed
   */
  async post<T = any>(url: string, data?: any, config?: AxiosRequestConfig): Promise<T> {
    const enhancedConfig = this.addPriorityHeaders(config);
    return super.post<T>(url, data, enhancedConfig);
  }
  
  /**
   * Override PUT request to add priority headers if needed
   */
  async put<T = any>(url: string, data?: any, config?: AxiosRequestConfig): Promise<T> {
    const enhancedConfig = this.addPriorityHeaders(config);
    return super.put<T>(url, data, enhancedConfig);
  }
  
  /**
   * Override PATCH request to add priority headers if needed
   */
  async patch<T = any>(url: string, data?: any, config?: AxiosRequestConfig): Promise<T> {
    const enhancedConfig = this.addPriorityHeaders(config);
    return super.patch<T>(url, data, enhancedConfig);
  }
  
  /**
   * Override DELETE request to add priority headers if needed
   */
  async delete<T = any>(url: string, config?: AxiosRequestConfig): Promise<T> {
    const enhancedConfig = this.addPriorityHeaders(config);
    return super.delete<T>(url, enhancedConfig);
  }
  
  /**
   * Add priority headers to request config if configured
   * @private
   */
  private addPriorityHeaders(config?: AxiosRequestConfig): AxiosRequestConfig | undefined {
    if (!this.priorityHeader || !this.priorityHeaderValue) {
      return config;
    }
    
    return {
      ...config,
      headers: {
        ...(config?.headers || {}),
        [this.priorityHeader]: this.priorityHeaderValue,
      },
    };
  }
}

// Create a factory function for easy client creation
export function createEnhancedApiClient(options: EnhancedClientOptions): EnhancedResilientApiClient {
  return new EnhancedResilientApiClient(options);
}
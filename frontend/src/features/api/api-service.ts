/**
 * API Service layer for handling connections to backend services
 * with circuit breaker pattern and retry mechanisms
 */

// Constants for API configuration
const API_TIMEOUT = 15000; // 15 seconds
const MAX_RETRIES = 3;
const RETRY_DELAY = 1000; // 1 second base delay
const CIRCUIT_BREAKER_THRESHOLD = 5; // Number of failures before opening the circuit
const CIRCUIT_BREAKER_RESET_TIMEOUT = 30000; // 30 seconds before resetting

// Circuit breaker state
type CircuitState = 'CLOSED' | 'OPEN' | 'HALF-OPEN';

interface CircuitBreaker {
  state: CircuitState;
  failureCount: number;
  lastFailureTime: number | null;
  services: Record<string, {
    state: CircuitState;
    failureCount: number;
    lastFailureTime: number | null;
  }>;
}

// Initialize circuit breaker
const circuitBreaker: CircuitBreaker = {
  state: 'CLOSED',
  failureCount: 0,
  lastFailureTime: null,
  services: {
    heroku: {
      state: 'CLOSED',
      failureCount: 0,
      lastFailureTime: null
    },
    rdkit: {
      state: 'CLOSED',
      failureCount: 0,
      lastFailureTime: null
    },
    convex: {
      state: 'CLOSED',
      failureCount: 0,
      lastFailureTime: null
    }
  }
};

// API Service configuration
interface ApiConfig {
  baseUrl: string;
  timeout?: number;
  headers?: Record<string, string>;
  service: 'heroku' | 'rdkit' | 'convex';
}

// Default configurations for different backend services
const apiConfigs: Record<string, ApiConfig> = {
  heroku: {
    baseUrl: process.env.NEXT_PUBLIC_API_URL || 'https://cryoprotect-8030e4025428.herokuapp.com/api',
    timeout: API_TIMEOUT,
    service: 'heroku',
    headers: {
      'Content-Type': 'application/json'
    }
  },
  rdkit: {
    baseUrl: process.env.NEXT_PUBLIC_RDKIT_URL || 'https://cryoprotect-rdkit.fly.dev',
    timeout: API_TIMEOUT * 2, // Longer timeout for computational services
    service: 'rdkit',
    headers: {
      'Content-Type': 'application/json'
    }
  },
  convex: {
    baseUrl: process.env.NEXT_PUBLIC_CONVEX_URL || 'https://adaptable-crocodile-452.convex.cloud',
    timeout: API_TIMEOUT,
    service: 'convex',
    headers: {
      'Content-Type': 'application/json'
    }
  }
};

// Request options interface
interface RequestOptions {
  method?: 'GET' | 'POST' | 'PUT' | 'DELETE' | 'PATCH';
  headers?: Record<string, string>;
  body?: any;
  timeoutMs?: number;
  retries?: number;
}

/**
 * Check if circuit is open for a specific service
 */
const isCircuitOpen = (service: string): boolean => {
  const serviceBreaker = circuitBreaker.services[service];
  if (!serviceBreaker) return false;
  
  if (serviceBreaker.state === 'OPEN') {
    // Check if it's time to try again (half-open state)
    if (serviceBreaker.lastFailureTime && 
        Date.now() - serviceBreaker.lastFailureTime > CIRCUIT_BREAKER_RESET_TIMEOUT) {
      serviceBreaker.state = 'HALF-OPEN';
      return false;
    }
    return true;
  }
  
  return false;
};

/**
 * Record a success and reset the circuit breaker if in half-open state
 */
const recordSuccess = (service: string): void => {
  const serviceBreaker = circuitBreaker.services[service];
  if (!serviceBreaker) return;
  
  if (serviceBreaker.state === 'HALF-OPEN') {
    serviceBreaker.state = 'CLOSED';
    serviceBreaker.failureCount = 0;
    serviceBreaker.lastFailureTime = null;
  }
};

/**
 * Record a failure and potentially open the circuit
 */
const recordFailure = (service: string): void => {
  const serviceBreaker = circuitBreaker.services[service];
  if (!serviceBreaker) return;
  
  serviceBreaker.failureCount++;
  serviceBreaker.lastFailureTime = Date.now();
  
  if (serviceBreaker.failureCount >= CIRCUIT_BREAKER_THRESHOLD) {
    serviceBreaker.state = 'OPEN';
    console.warn(`Circuit breaker opened for service: ${service}`);
  }
};

/**
 * Get circuit breaker status for all services
 */
export const getCircuitStatus = (): Record<string, CircuitState> => {
  const status: Record<string, CircuitState> = {};
  
  Object.keys(circuitBreaker.services).forEach(service => {
    status[service] = circuitBreaker.services[service].state;
  });
  
  return status;
};

/**
 * Make an API request with circuit breaker and retry mechanism
 */
export const apiRequest = async <T>(
  endpoint: string,
  options: RequestOptions = {},
  serviceKey: 'heroku' | 'rdkit' | 'convex' = 'heroku'
): Promise<T> => {
  const config = apiConfigs[serviceKey];
  
  if (!config) {
    throw new Error(`Service configuration not found for ${serviceKey}`);
  }
  
  // Check if circuit is open
  if (isCircuitOpen(serviceKey)) {
    throw new Error(`Circuit is open for ${serviceKey} service. Please try again later.`);
  }
  
  const {
    method = 'GET',
    headers = {},
    body,
    timeoutMs = config.timeout,
    retries = MAX_RETRIES
  } = options;
  
  const url = `${config.baseUrl}/${endpoint}`.replace(/\/+/g, '/').replace(':/', '://');
  
  // Merge headers
  const mergedHeaders = { ...config.headers, ...headers };
  
  // Create request options
  const fetchOptions: RequestInit = {
    method,
    headers: mergedHeaders,
    body: body ? JSON.stringify(body) : undefined
  };
  
  // Create AbortController for timeout
  const controller = new AbortController();
  const timeoutId = setTimeout(() => controller.abort(), timeoutMs);
  
  // Add signal to fetch options
  fetchOptions.signal = controller.signal;
  
  // Initialize retry counter
  let currentRetry = 0;
  
  while (true) {
    try {
      const response = await fetch(url, fetchOptions);
      
      // Clear timeout
      clearTimeout(timeoutId);
      
      // Handle response
      if (!response.ok) {
        const errorData = await response.json().catch(() => null);
        
        // Decide whether to retry based on status code
        if (currentRetry < retries && (response.status >= 500 || response.status === 429)) {
          currentRetry++;
          const delay = RETRY_DELAY * Math.pow(2, currentRetry - 1); // Exponential backoff
          await new Promise(resolve => setTimeout(resolve, delay));
          continue;
        }
        
        // Record failure after retries are exhausted
        recordFailure(serviceKey);
        
        throw {
          status: response.status,
          statusText: response.statusText,
          data: errorData
        };
      }
      
      // Parse response
      const data = await response.json();
      
      // Record success
      recordSuccess(serviceKey);
      
      return data as T;
    } catch (error: any) {
      // Clear timeout
      clearTimeout(timeoutId);
      
      // Handle timeout/abort errors
      if (error.name === 'AbortError') {
        if (currentRetry < retries) {
          currentRetry++;
          const delay = RETRY_DELAY * Math.pow(2, currentRetry - 1);
          await new Promise(resolve => setTimeout(resolve, delay));
          continue;
        }
        
        recordFailure(serviceKey);
        throw new Error(`Request timed out after ${timeoutMs}ms`);
      }
      
      // If we've exhausted retries, record failure
      if (currentRetry >= retries) {
        recordFailure(serviceKey);
      }
      
      // For network errors that should be retried
      if (error instanceof TypeError && error.message.includes('fetch')) {
        if (currentRetry < retries) {
          currentRetry++;
          const delay = RETRY_DELAY * Math.pow(2, currentRetry - 1);
          await new Promise(resolve => setTimeout(resolve, delay));
          continue;
        }
      }
      
      // Rethrow error
      throw error;
    }
  }
};

/**
 * API Service for experimental data
 */
export const experimentalDataApi = {
  // Experiments
  getExperiments: () => apiRequest<any[]>('experiments'),
  getExperiment: (id: string) => apiRequest<any>(`experiments/${id}`),
  createExperiment: (data: any) => apiRequest<any>('experiments', { method: 'POST', body: data }),
  updateExperiment: (id: string, data: any) => apiRequest<any>(`experiments/${id}`, { method: 'PUT', body: data }),
  deleteExperiment: (id: string) => apiRequest<void>(`experiments/${id}`, { method: 'DELETE' }),
  
  // Protocols
  getProtocols: () => apiRequest<any[]>('protocols'),
  getProtocol: (id: string) => apiRequest<any>(`protocols/${id}`),
  createProtocol: (data: any) => apiRequest<any>('protocols', { method: 'POST', body: data }),
  updateProtocol: (id: string, data: any) => apiRequest<any>(`protocols/${id}`, { method: 'PUT', body: data }),
  deleteProtocol: (id: string) => apiRequest<void>(`protocols/${id}`, { method: 'DELETE' }),
  getProtocolVersions: (id: string) => apiRequest<any[]>(`protocols/${id}/versions`),
  
  // Results
  getExperimentResults: (id: string) => apiRequest<any[]>(`experiments/${id}/results`),
  addExperimentResult: (id: string, data: any) => apiRequest<any>(`experiments/${id}/results`, { method: 'POST', body: data }),
  
  // RDKit Services - uses the rdkit service
  getMolecularProperties: (smiles: string) => apiRequest<any>(`properties?smiles=${encodeURIComponent(smiles)}`, {}, 'rdkit'),
  getVisualizationData: (smiles: string, format: string = 'svg') => 
    apiRequest<any>(`visualization?smiles=${encodeURIComponent(smiles)}&format=${format}`, {}, 'rdkit'),
};

/**
 * API Service for backup direct connection to Supabase
 * Used when Heroku or Convex services are unavailable
 */
export const backupDirectApi = {
  // Implements the same interface as experimentalDataApi but uses direct connections
  // This would be implemented with direct database queries in a real application
};

/**
 * Get appropriate API service based on circuit status
 */
export const getApiService = () => {
  const status = getCircuitStatus();
  
  // Use the main API service if Heroku is available
  if (status.heroku !== 'OPEN') {
    return experimentalDataApi;
  }
  
  // Fall back to direct connection if Heroku is unavailable
  console.warn('Using backup direct API connection due to Heroku service unavailability');
  return backupDirectApi;
};
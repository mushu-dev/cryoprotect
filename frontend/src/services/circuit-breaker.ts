/**
 * Circuit Breaker Service
 * 
 * This service implements the Circuit Breaker pattern for API calls,
 * preventing cascading failures and allowing graceful degradation.
 * 
 * The circuit breaker has three states:
 * - CLOSED: All requests are allowed through (normal operation)
 * - OPEN: All requests are rejected immediately (failure state)
 * - HALF_OPEN: A limited number of test requests are allowed to check if the service has recovered
 */

// Circuit breaker states
export enum CircuitState {
  CLOSED = 'CLOSED',
  OPEN = 'OPEN',
  HALF_OPEN = 'HALF_OPEN',
}

// Options for circuit breaker configuration
export interface CircuitBreakerOptions {
  /** Failure threshold to trip the circuit (default: 5) */
  failureThreshold?: number;
  
  /** Time window for failure threshold in milliseconds (default: 60000) */
  failureThresholdTimeWindow?: number;
  
  /** Reset timeout in milliseconds (default: 30000) */
  resetTimeout?: number;
  
  /** Number of successful test requests required to close the circuit (default: 2) */
  successThreshold?: number;
  
  /** Maximum number of concurrent requests in half-open state (default: 1) */
  halfOpenMaxConcurrent?: number;
  
  /** Maximum queue size for delayed execution (default: 10) */
  maxQueueSize?: number;
  
  /** Whether to enable fallbacks for circuit-broken requests (default: true) */
  enableFallbacks?: boolean;
}

// Circuit breaker state snapshot
export interface CircuitBreakerSnapshot {
  state: CircuitState;
  failures: number;
  successes: number;
  lastFailure: number | null;
  lastSuccess: number | null;
  lastStateChange: number;
  halfOpenRequests: number;
  startTime: number;
  totalRequests: number;
  totalFailures: number;
  totalSuccesses: number;
  queueSize: number;
}

// Interface for queue entry
interface QueueEntry<T> {
  timestamp: number;
  resolve: (value: T) => void;
  reject: (error: Error) => void;
  operation: () => Promise<T>;
  fallback?: () => Promise<T>;
}

/**
 * Circuit Breaker implementation
 */
export class CircuitBreaker {
  // Circuit state
  private state: CircuitState = CircuitState.CLOSED;
  
  // Circuit configuration
  private failureThreshold: number;
  private failureThresholdTimeWindow: number;
  private resetTimeout: number;
  private successThreshold: number;
  private halfOpenMaxConcurrent: number;
  private maxQueueSize: number;
  private enableFallbacks: boolean;
  
  // Circuit metrics
  private failures: number = 0;
  private successes: number = 0;
  private lastFailure: number | null = null;
  private lastSuccess: number | null = null;
  private lastStateChange: number = Date.now();
  private halfOpenRequests: number = 0;
  private startTime: number = Date.now();
  private totalRequests: number = 0;
  private totalFailures: number = 0;
  private totalSuccesses: number = 0;
  
  // Request queue for when circuit is open
  private queue: Array<QueueEntry<any>> = [];
  
  // State change listeners
  private stateChangeListeners: Array<(state: CircuitState, snapshot: CircuitBreakerSnapshot) => void> = [];
  
  // Fallback functions for different operations
  private fallbacks: Map<string, (error: Error) => Promise<any>> = new Map();
  
  /**
   * Create a new circuit breaker
   */
  constructor(name: string, options: CircuitBreakerOptions = {}) {
    this.failureThreshold = options.failureThreshold || 5;
    this.failureThresholdTimeWindow = options.failureThresholdTimeWindow || 60000; // 1 minute
    this.resetTimeout = options.resetTimeout || 30000; // 30 seconds
    this.successThreshold = options.successThreshold || 2;
    this.halfOpenMaxConcurrent = options.halfOpenMaxConcurrent || 1;
    this.maxQueueSize = options.maxQueueSize || 10;
    this.enableFallbacks = options.enableFallbacks !== false; // Default to true
    
    // Start the health check timer
    this.startHealthCheck();
  }
  
  /**
   * Execute a function with circuit breaker protection
   * @param operationKey A key to identify this operation (for fallbacks)
   * @param operation The function to execute
   * @param fallback Optional fallback function to execute if the circuit is open
   */
  async execute<T>(operationKey: string, operation: () => Promise<T>, fallback?: () => Promise<T>): Promise<T> {
    this.totalRequests++;
    
    // Check circuit state
    switch (this.state) {
      case CircuitState.OPEN:
        return this.handleOpenCircuit(operationKey, operation, fallback);
        
      case CircuitState.HALF_OPEN:
        return this.handleHalfOpenCircuit(operationKey, operation, fallback);
        
      case CircuitState.CLOSED:
      default:
        return this.handleClosedCircuit(operationKey, operation, fallback);
    }
  }
  
  /**
   * Register a fallback function for an operation
   */
  registerFallback<T>(operationKey: string, fallback: (error: Error) => Promise<T>): void {
    this.fallbacks.set(operationKey, fallback);
  }
  
  /**
   * Add a listener for state changes
   */
  onStateChange(listener: (state: CircuitState, snapshot: CircuitBreakerSnapshot) => void): () => void {
    this.stateChangeListeners.push(listener);
    
    // Call immediately with current state
    listener(this.state, this.getSnapshot());
    
    // Return a function to remove the listener
    return () => {
      const index = this.stateChangeListeners.indexOf(listener);
      if (index !== -1) {
        this.stateChangeListeners.splice(index, 1);
      }
    };
  }
  
  /**
   * Get a snapshot of the current circuit state
   */
  getSnapshot(): CircuitBreakerSnapshot {
    return {
      state: this.state,
      failures: this.failures,
      successes: this.successes,
      lastFailure: this.lastFailure,
      lastSuccess: this.lastSuccess,
      lastStateChange: this.lastStateChange,
      halfOpenRequests: this.halfOpenRequests,
      startTime: this.startTime,
      totalRequests: this.totalRequests,
      totalFailures: this.totalFailures,
      totalSuccesses: this.totalSuccesses,
      queueSize: this.queue.length,
    };
  }
  
  /**
   * Reset the circuit breaker to closed state
   */
  reset(): void {
    const previousState = this.state;
    
    this.state = CircuitState.CLOSED;
    this.failures = 0;
    this.successes = 0;
    this.halfOpenRequests = 0;
    this.lastStateChange = Date.now();
    
    // Process any queued requests
    this.processQueue();
    
    // Notify listeners if state changed
    if (previousState !== CircuitState.CLOSED) {
      this.notifyStateChange();
    }
  }
  
  /**
   * Trip the circuit breaker (force open)
   */
  trip(): void {
    if (this.state !== CircuitState.OPEN) {
      this.state = CircuitState.OPEN;
      this.lastStateChange = Date.now();
      this.notifyStateChange();
    }
  }
  
  /**
   * Handle a closed circuit request
   * @private
   */
  private async handleClosedCircuit<T>(
    operationKey: string, 
    operation: () => Promise<T>, 
    fallback?: () => Promise<T>
  ): Promise<T> {
    try {
      // Execute the operation
      const result = await operation();
      
      // Record success
      this.recordSuccess();
      return result;
    } catch (error) {
      // Record failure
      this.recordFailure();
      
      // Check if we should trip the circuit
      if (this.shouldTripCircuit()) {
        this.tripCircuit();
      }
      
      // Try to execute fallback if available
      return this.tryFallback(operationKey, error as Error, fallback);
    }
  }
  
  /**
   * Handle a half-open circuit request
   * @private
   */
  private async handleHalfOpenCircuit<T>(
    operationKey: string, 
    operation: () => Promise<T>, 
    fallback?: () => Promise<T>
  ): Promise<T> {
    // Check if we can make a test request
    if (this.halfOpenRequests >= this.halfOpenMaxConcurrent) {
      // Too many concurrent requests, queue this one if enabled
      if (this.queue.length < this.maxQueueSize) {
        return this.queueRequest(operation, fallback);
      }
      
      // Queue is full, try fallback
      return this.tryFallback(operationKey, new Error('Circuit is half-open and queue is full'), fallback);
    }
    
    // Increment half-open requests counter
    this.halfOpenRequests++;
    
    try {
      // Execute test request
      const result = await operation();
      
      // Record success
      this.recordSuccess();
      
      // Check if we should close the circuit
      if (this.successes >= this.successThreshold) {
        this.closeCircuit();
      }
      
      return result;
    } catch (error) {
      // Record failure and trip the circuit again
      this.recordFailure();
      this.tripCircuit();
      
      // Try fallback
      return this.tryFallback(operationKey, error as Error, fallback);
    } finally {
      // Decrement half-open requests counter
      this.halfOpenRequests--;
    }
  }
  
  /**
   * Handle an open circuit request
   * @private
   */
  private async handleOpenCircuit<T>(
    operationKey: string, 
    operation: () => Promise<T>, 
    fallback?: () => Promise<T>
  ): Promise<T> {
    // Check if the reset timeout has elapsed
    if (this.shouldAttemptReset()) {
      // Transition to half-open state
      this.state = CircuitState.HALF_OPEN;
      this.lastStateChange = Date.now();
      this.successes = 0;
      this.failures = 0;
      this.notifyStateChange();
      
      // Process this request as a half-open request
      return this.handleHalfOpenCircuit(operationKey, operation, fallback);
    }
    
    // If the queue isn't full, queue the request
    if (this.queue.length < this.maxQueueSize) {
      return this.queueRequest(operation, fallback);
    }
    
    // Circuit is open and queue is full, try fallback
    return this.tryFallback(
      operationKey, 
      new Error(`Circuit is open for another ${Math.round((this.lastStateChange + this.resetTimeout - Date.now()) / 1000)}s`), 
      fallback
    );
  }
  
  /**
   * Record a successful operation
   * @private
   */
  private recordSuccess(): void {
    this.successes++;
    this.totalSuccesses++;
    this.lastSuccess = Date.now();
  }
  
  /**
   * Record a failed operation
   * @private
   */
  private recordFailure(): void {
    this.failures++;
    this.totalFailures++;
    this.lastFailure = Date.now();
  }
  
  /**
   * Check if the circuit should be tripped
   * @private
   */
  private shouldTripCircuit(): boolean {
    // Check if we've reached the failure threshold
    if (this.failures < this.failureThreshold) {
      return false;
    }
    
    // Check if the failures occurred within the time window
    const now = Date.now();
    const windowStart = now - this.failureThresholdTimeWindow;
    
    // If the first failure is within the window, trip the circuit
    return this.lastFailure !== null && this.lastFailure >= windowStart;
  }
  
  /**
   * Trip the circuit (open it)
   * @private
   */
  private tripCircuit(): void {
    if (this.state !== CircuitState.OPEN) {
      this.state = CircuitState.OPEN;
      this.lastStateChange = Date.now();
      this.notifyStateChange();
    }
  }
  
  /**
   * Close the circuit
   * @private
   */
  private closeCircuit(): void {
    if (this.state !== CircuitState.CLOSED) {
      this.state = CircuitState.CLOSED;
      this.lastStateChange = Date.now();
      this.failures = 0;
      this.successes = 0;
      this.halfOpenRequests = 0;
      this.notifyStateChange();
      
      // Process any queued requests
      this.processQueue();
    }
  }
  
  /**
   * Check if we should attempt to reset the circuit
   * @private
   */
  private shouldAttemptReset(): boolean {
    if (this.state !== CircuitState.OPEN) {
      return false;
    }
    
    const now = Date.now();
    return (now - this.lastStateChange) >= this.resetTimeout;
  }
  
  /**
   * Queue a request for later execution
   * @private
   */
  private queueRequest<T>(operation: () => Promise<T>, fallback?: () => Promise<T>): Promise<T> {
    return new Promise<T>((resolve, reject) => {
      // Add to queue
      this.queue.push({
        timestamp: Date.now(),
        resolve,
        reject,
        operation,
        fallback,
      });
    });
  }
  
  /**
   * Process queued requests
   * @private
   */
  private processQueue(): void {
    // Process all requests in the queue
    const requests = [...this.queue];
    this.queue = [];
    
    for (const request of requests) {
      // Execute the operation
      this.execute('queued-request', request.operation, request.fallback)
        .then(request.resolve)
        .catch(request.reject);
    }
  }
  
  /**
   * Try to execute a fallback
   * @private
   */
  private async tryFallback<T>(
    operationKey: string, 
    error: Error, 
    providedFallback?: () => Promise<T>
  ): Promise<T> {
    // Check if fallbacks are enabled
    if (!this.enableFallbacks) {
      throw error;
    }
    
    try {
      // Try the provided fallback first
      if (providedFallback) {
        return await providedFallback();
      }
      
      // Try a registered fallback
      const registeredFallback = this.fallbacks.get(operationKey);
      if (registeredFallback) {
        return await registeredFallback(error);
      }
      
      // No fallback available
      throw error;
    } catch (fallbackError) {
      // Fallback also failed
      if (fallbackError !== error) {
        console.error('Fallback failed:', fallbackError);
      }
      
      throw error;
    }
  }
  
  /**
   * Notify listeners of state changes
   * @private
   */
  private notifyStateChange(): void {
    const snapshot = this.getSnapshot();
    
    for (const listener of this.stateChangeListeners) {
      try {
        listener(this.state, snapshot);
      } catch (error) {
        console.error('Error in circuit breaker state change listener:', error);
      }
    }
  }
  
  /**
   * Start the health check timer
   * @private
   */
  private startHealthCheck(): void {
    // Check every resetTimeout / 2 milliseconds
    const healthCheckInterval = Math.max(1000, Math.floor(this.resetTimeout / 2));
    
    setInterval(() => {
      // If the circuit is open and reset timeout has elapsed, transition to half-open
      if (this.state === CircuitState.OPEN && this.shouldAttemptReset()) {
        this.state = CircuitState.HALF_OPEN;
        this.lastStateChange = Date.now();
        this.successes = 0;
        this.failures = 0;
        this.notifyStateChange();
      }
      
      // Clear expired entries from the request queue
      const now = Date.now();
      const maxAge = this.resetTimeout * 2;
      
      this.queue = this.queue.filter(entry => {
        const age = now - entry.timestamp;
        
        if (age > maxAge) {
          // Reject the request with a timeout error
          entry.reject(new Error('Request timeout - circuit was open too long'));
          return false;
        }
        
        return true;
      });
    }, healthCheckInterval);
  }
}

// Map of circuit breakers by name
const circuitBreakers = new Map<string, CircuitBreaker>();

/**
 * Get or create a circuit breaker instance by name
 */
export function getCircuitBreaker(name: string, options?: CircuitBreakerOptions): CircuitBreaker {
  if (!circuitBreakers.has(name)) {
    circuitBreakers.set(name, new CircuitBreaker(name, options));
  }
  
  return circuitBreakers.get(name)!;
}

/**
 * Reset all circuit breakers
 */
export function resetAllCircuitBreakers(): void {
  for (const breaker of circuitBreakers.values()) {
    breaker.reset();
  }
}

/**
 * Get all circuit breakers
 */
export function getAllCircuitBreakers(): Map<string, CircuitBreaker> {
  return new Map(circuitBreakers);
}
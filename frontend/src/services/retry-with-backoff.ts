/**
 * Retry utility with exponential backoff
 * 
 * This module provides a utility for retrying asynchronous operations
 * with exponential backoff to handle transient failures.
 */

export interface RetryOptions {
  /** Maximum number of retry attempts (default: 3) */
  maxRetries?: number;
  
  /** Initial delay in milliseconds (default: 300) */
  initialDelay?: number;
  
  /** Maximum delay in milliseconds (default: 10000) */
  maxDelay?: number;
  
  /** Backoff factor to multiply delay by after each attempt (default: 2) */
  backoffFactor?: number;
  
  /** Jitter factor to randomize delay (0-1, default: 0.1) */
  jitterFactor?: number;
  
  /** Function to determine if an error is retryable (default: all errors are retryable) */
  isRetryable?: (error: Error) => boolean;
  
  /** Function called before each retry attempt */
  onRetry?: (error: Error, attempt: number, delay: number) => void;
  
  /** Whether to abort retries when navigator.onLine is false (default: true) */
  respectNavigatorOnline?: boolean;
}

/**
 * Retry an async operation with exponential backoff
 * 
 * @param operation - The async function to retry
 * @param options - Retry options
 * @returns Result of the operation
 * @throws Last encountered error if all retries fail
 */
export async function retryWithBackoff<T>(
  operation: () => Promise<T>,
  options: RetryOptions = {}
): Promise<T> {
  // Set default options
  const {
    maxRetries = 3,
    initialDelay = 300,
    maxDelay = 10000,
    backoffFactor = 2,
    jitterFactor = 0.1,
    isRetryable = () => true,
    onRetry,
    respectNavigatorOnline = true
  } = options;
  
  let lastError: Error;
  
  // Try the operation up to maxRetries + 1 times (initial + retries)
  for (let attempt = 0; attempt <= maxRetries; attempt++) {
    try {
      // Execute the operation
      return await operation();
    } catch (error) {
      // Cast to Error type
      lastError = error instanceof Error ? error : new Error(String(error));
      
      // Check if we should abort retries
      const isLastAttempt = attempt === maxRetries;
      const isNetworkOffline = respectNavigatorOnline && typeof navigator !== 'undefined' && !navigator.onLine;
      const shouldNotRetry = !isRetryable(lastError);
      
      if (isLastAttempt || isNetworkOffline || shouldNotRetry) {
        throw lastError;
      }
      
      // Calculate delay with exponential backoff and jitter
      const baseDelay = Math.min(initialDelay * Math.pow(backoffFactor, attempt), maxDelay);
      const jitter = 1 - jitterFactor + (Math.random() * jitterFactor * 2);
      const delay = Math.floor(baseDelay * jitter);
      
      // Notify about retry
      if (onRetry) {
        onRetry(lastError, attempt + 1, delay);
      }
      
      // Wait before retrying
      await new Promise(resolve => setTimeout(resolve, delay));
    }
  }
  
  // This should never be reached due to the throw in the loop
  throw lastError!;
}

/**
 * Decorator that adds retry with backoff to an async function
 * 
 * @param options - Retry options
 * @returns Decorated function with retry behavior
 */
export function withRetry<T extends (...args: any[]) => Promise<any>>(
  fn: T,
  options: RetryOptions = {}
): T {
  return (async (...args: Parameters<T>): Promise<ReturnType<T>> => {
    return retryWithBackoff(() => fn(...args), options);
  }) as T;
}

export default retryWithBackoff;
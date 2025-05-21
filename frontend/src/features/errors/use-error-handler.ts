/**
 * Hook for standardized error handling throughout the application
 * 
 * Provides consistent error handling with logging and user feedback.
 */

import { useCallback } from 'react';
import { ConvexError } from 'convex/values';
import { useToast } from '../ui/use-toast';

// Common error types
export type ErrorWithCode = Error & { code?: string };
export type ErrorWithStatus = Error & { status?: number };
export type ErrorWithCause = Error & { cause?: Error };
export type EnhancedError = ErrorWithCode & ErrorWithStatus & ErrorWithCause;

/**
 * Hook for standardized error handling
 */
export function useErrorHandler() {
  const { toast } = useToast();
  
  /**
   * Process an error with logging and user feedback
   */
  const handleError = useCallback((
    error: unknown,
    context: string,
    options?: {
      showToast?: boolean;
      logToConsole?: boolean;
      reportToErrorService?: boolean;
    }
  ) => {
    const showToast = options?.showToast ?? true;
    const logToConsole = options?.logToConsole ?? true;
    const reportToErrorService = options?.reportToErrorService ?? true;
    
    // Default error message
    let errorMessage = context || 'An error occurred';
    let errorDetails = '';
    let errorCode: string | undefined;
    let statusCode: number | undefined;
    
    // Extract information from different error types
    if (error instanceof ConvexError) {
      // Convex-specific errors
      errorMessage = `${context}: ${error.message}`;
      errorDetails = error.data ? JSON.stringify(error.data) : '';
    } else if (error instanceof Error) {
      // Standard Error object
      errorMessage = `${context}: ${error.message}`;
      
      // Extract additional information from enhanced errors
      const enhancedError = error as EnhancedError;
      errorCode = enhancedError.code;
      statusCode = enhancedError.status;
      
      // Include cause if available
      if (enhancedError.cause) {
        errorDetails = `Caused by: ${enhancedError.cause.message}`;
      }
      
      // Include stack trace in development
      if (process.env.NODE_ENV === 'development') {
        errorDetails += `\n${error.stack}`;
      }
    } else if (typeof error === 'string') {
      // String error message
      errorMessage = `${context}: ${error}`;
    } else {
      // Unknown error type
      errorMessage = `${context}: Unknown error`;
      errorDetails = error ? JSON.stringify(error) : '';
    }
    
    // Console logging
    if (logToConsole) {
      console.error(`[ERROR] ${errorMessage}`, error);
      if (errorDetails) {
        console.error(`Details: ${errorDetails}`);
      }
    }
    
    // Show toast notification to user
    if (showToast) {
      toast({
        title: "Error",
        description: errorMessage,
        variant: "destructive",
      });
    }
    
    // Report to error monitoring service (if configured)
    if (reportToErrorService && window.errorReportingService) {
      window.errorReportingService.captureException(error, {
        tags: {
          context,
          errorCode,
          statusCode
        },
        extra: {
          details: errorDetails
        }
      });
    }
    
    // Return the processed error information
    return {
      message: errorMessage,
      details: errorDetails,
      code: errorCode,
      status: statusCode,
      originalError: error
    };
  }, [toast]);
  
  return { handleError };
}

// Mock error reporting service interface for typesafety
declare global {
  interface Window {
    errorReportingService?: {
      captureException: (error: unknown, options?: any) => void;
    };
  }
}
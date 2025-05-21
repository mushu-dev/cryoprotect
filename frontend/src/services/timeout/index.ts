/**
 * Timeout Handling Module
 * 
 * This module exports all components, hooks, and services related to
 * enhanced timeout handling in API requests.
 */

import { TimeoutController, TimeoutConfig, getTimeoutController } from '../timeout-controller';
import { EnhancedResilientApiClient, createEnhancedApiClient } from '../enhanced-resilient-api-client';
import RequestTimeoutMonitor from '@/components/request-timeout-monitor';
import useApiTimeout from '@/hooks/use-api-timeout';

// Re-export everything
export {
  // Services
  TimeoutController,
  getTimeoutController,
  EnhancedResilientApiClient,
  createEnhancedApiClient,
  
  // Components
  RequestTimeoutMonitor,
  
  // Hooks
  useApiTimeout,
};

// Export types
export type { TimeoutConfig };

// Default export with all components
const TimeoutHandling = {
  TimeoutController,
  getTimeoutController,
  EnhancedResilientApiClient,
  createEnhancedApiClient,
  RequestTimeoutMonitor,
  useApiTimeout,
};

export default TimeoutHandling;
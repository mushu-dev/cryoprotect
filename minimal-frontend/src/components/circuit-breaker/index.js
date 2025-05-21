/**
 * Circuit Breaker UI Components
 * 
 * This module exports UI components for displaying circuit breaker status
 * and interacting with circuit breakers in the application.
 */

export { ContextAwareIndicator } from './context-aware-indicator';
export { 
  CircuitBreakerProvider,
  useCircuitBreaker,
  CIRCUIT_STATE
} from './circuit-provider';
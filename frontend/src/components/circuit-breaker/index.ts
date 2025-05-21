/**
 * Circuit Breaker UI Components
 * 
 * This module exports UI components for displaying circuit breaker status
 * and interacting with circuit breakers in the application.
 */

export { CircuitBreakerStatus } from '../circuit-breaker-status';
export { CircuitBreakerDashboard } from '../circuit-breaker-dashboard';
export { CircuitBreakerIndicator } from '../circuit-breaker-indicator';
export { ContextAwareIndicator } from './context-aware-indicator';
export { 
  CircuitBreakerProvider,
  useCircuitBreaker as useCircuitBreakers,
  CIRCUIT_STATE
} from './circuit-provider';
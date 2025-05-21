/**
 * Context-Aware Circuit Breaker Indicator
 * 
 * A visual indicator that displays the current state of a service's circuit breaker.
 */
import React from 'react';
import { useCircuitBreaker, CIRCUIT_STATE } from './circuit-provider';

/**
 * Context-aware circuit breaker indicator component
 * 
 * @param {Object} props - Component props
 * @param {string} props.serviceName - Name of the service to monitor
 * @param {boolean} props.showLabel - Whether to show a text label alongside the indicator
 */
export function ContextAwareIndicator({ serviceName, showLabel = false }) {
  const { getServiceState } = useCircuitBreaker();
  const state = getServiceState(serviceName);
  
  // Determine color based on circuit state
  const getColor = () => {
    switch (state) {
      case CIRCUIT_STATE.CLOSED:
        return 'bg-green-500';
      case CIRCUIT_STATE.HALF_OPEN:
        return 'bg-yellow-500';
      case CIRCUIT_STATE.OPEN:
        return 'bg-red-500';
      default:
        return 'bg-gray-500';
    }
  };
  
  // Determine label text based on state
  const getLabel = () => {
    switch (state) {
      case CIRCUIT_STATE.CLOSED:
        return 'Healthy';
      case CIRCUIT_STATE.HALF_OPEN:
        return 'Testing';
      case CIRCUIT_STATE.OPEN:
        return 'Failing';
      default:
        return 'Unknown';
    }
  };
  
  return (
    <div className="flex items-center">
      <div className={`w-3 h-3 rounded-full ${getColor()}`} />
      {showLabel && (
        <span className="ml-2 text-sm font-medium">
          {getLabel()}
        </span>
      )}
    </div>
  );
}
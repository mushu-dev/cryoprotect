import React from 'react';
import { useCircuitBreaker, CIRCUIT_STATE } from './circuit-breaker/circuit-provider';

export function CircuitBreakerIndicator({ 
  serviceName, 
  showLabel = false, 
  showFailures = false,
  variant = 'pill'
}) {
  const { services } = useCircuitBreaker();
  const service = services[serviceName];
  
  if (!service) {
    return null;
  }
  
  const { state, failures } = service;
  
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
  
  // Determine label based on state
  const getLabel = () => {
    switch (state) {
      case CIRCUIT_STATE.CLOSED:
        return 'Healthy';
      case CIRCUIT_STATE.HALF_OPEN:
        return 'Testing';
      case CIRCUIT_STATE.OPEN:
        return 'Failed';
      default:
        return 'Unknown';
    }
  };
  
  // Pill variant (small circle indicator)
  if (variant === 'pill') {
    return (
      <div className="flex items-center">
        <div className={`w-3 h-3 rounded-full ${getColor()}`} />
        {showLabel && (
          <span className="ml-2 text-sm font-medium">
            {getLabel()}
            {showFailures && failures > 0 && ` (${failures})`}
          </span>
        )}
      </div>
    );
  }
  
  // Badge variant (full badge)
  if (variant === 'badge') {
    const bgColor = 
      state === CIRCUIT_STATE.CLOSED ? 'bg-green-100 text-green-800 dark:bg-green-900 dark:text-green-300' :
      state === CIRCUIT_STATE.HALF_OPEN ? 'bg-yellow-100 text-yellow-800 dark:bg-yellow-900 dark:text-yellow-300' :
      'bg-red-100 text-red-800 dark:bg-red-900 dark:text-red-300';
    
    return (
      <span className={`inline-flex items-center rounded-full px-2.5 py-0.5 text-xs font-medium ${bgColor}`}>
        <span className={`mr-1.5 h-2 w-2 rounded-full ${getColor()}`}></span>
        {getLabel()}
        {showFailures && failures > 0 && ` (${failures})`}
      </span>
    );
  }
  
  return null;
}
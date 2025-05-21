import React from 'react';
import { useCircuitBreaker, CIRCUIT_STATE } from './circuit-provider';

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
  
  return (
    <div className="flex items-center">
      <div className={`w-3 h-3 rounded-full ${getColor()}`} />
      {showLabel && (
        <span className="ml-2 text-sm font-medium">
          {state === CIRCUIT_STATE.CLOSED ? 'Healthy' : 
           state === CIRCUIT_STATE.HALF_OPEN ? 'Testing' : 'Failing'}
        </span>
      )}
    </div>
  );
}
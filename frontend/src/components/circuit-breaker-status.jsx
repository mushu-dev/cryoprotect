import React from 'react';
import { useCircuitBreaker, CIRCUIT_STATE } from './circuit-breaker/circuit-provider';

export function CircuitBreakerStatus() {
  const { services, resetService } = useCircuitBreaker();
  
  // Get all service names
  const serviceNames = Object.keys(services);
  
  if (serviceNames.length === 0) {
    return (
      <div className="p-4 bg-white dark:bg-gray-800 rounded-lg shadow">
        <h3 className="text-lg font-medium">Circuit Breaker Status</h3>
        <p className="text-gray-500 dark:text-gray-400 mt-2">No services monitored yet</p>
      </div>
    );
  }
  
  return (
    <div className="p-4 bg-white dark:bg-gray-800 rounded-lg shadow">
      <h3 className="text-lg font-medium">Circuit Breaker Status</h3>
      
      <div className="mt-4 space-y-3">
        {serviceNames.map(serviceName => {
          const service = services[serviceName];
          const { state, failures } = service;
          
          return (
            <div 
              key={serviceName}
              className="flex items-center justify-between p-3 bg-gray-50 dark:bg-gray-700 rounded-md"
            >
              <div>
                <p className="font-medium">{serviceName}</p>
                <div className="flex items-center mt-1">
                  <div 
                    className={`w-3 h-3 rounded-full mr-2 ${
                      state === CIRCUIT_STATE.CLOSED ? 'bg-green-500' :
                      state === CIRCUIT_STATE.HALF_OPEN ? 'bg-yellow-500' :
                      'bg-red-500'
                    }`} 
                  />
                  <p className="text-sm text-gray-600 dark:text-gray-300">
                    {state === CIRCUIT_STATE.CLOSED ? 'Closed' :
                     state === CIRCUIT_STATE.HALF_OPEN ? 'Half-Open' :
                     'Open'} 
                    {failures > 0 && ` (${failures} failures)`}
                  </p>
                </div>
              </div>
              
              {state !== CIRCUIT_STATE.CLOSED && (
                <button
                  onClick={() => resetService(serviceName)}
                  className="px-3 py-1.5 text-xs font-medium text-white bg-blue-600 rounded-md hover:bg-blue-700 transition-colors"
                >
                  Reset
                </button>
              )}
            </div>
          );
        })}
      </div>
    </div>
  );
}
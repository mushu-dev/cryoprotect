import React, { useState } from 'react';
import { useCircuitBreaker, CIRCUIT_STATE } from './circuit-breaker/circuit-provider';

export function CircuitBreakerDashboard() {
  const { services, resetService } = useCircuitBreaker();
  const [filter, setFilter] = useState('all');
  
  // Get all service names
  const serviceNames = Object.keys(services);
  
  if (serviceNames.length === 0) {
    return (
      <div className="p-4 bg-white dark:bg-gray-800 rounded-lg shadow">
        <h3 className="text-lg font-medium">Service Health Dashboard</h3>
        <p className="text-gray-500 dark:text-gray-400 mt-2">No services monitored yet</p>
      </div>
    );
  }
  
  // Filter services based on the current filter
  const filteredServices = serviceNames.filter(serviceName => {
    if (filter === 'all') return true;
    if (filter === 'closed') return services[serviceName].state === CIRCUIT_STATE.CLOSED;
    if (filter === 'half-open') return services[serviceName].state === CIRCUIT_STATE.HALF_OPEN;
    if (filter === 'open') return services[serviceName].state === CIRCUIT_STATE.OPEN;
    return true;
  });
  
  // Count services by state
  const counts = serviceNames.reduce((acc, serviceName) => {
    const state = services[serviceName].state;
    acc[state] = (acc[state] || 0) + 1;
    return acc;
  }, {});
  
  return (
    <div className="p-4 bg-white dark:bg-gray-800 rounded-lg shadow">
      <div className="flex flex-col sm:flex-row sm:items-center sm:justify-between mb-4">
        <h3 className="text-lg font-medium">Service Health Dashboard</h3>
        
        <div className="mt-3 sm:mt-0 flex space-x-2">
          <button
            onClick={() => setFilter('all')}
            className={`px-3 py-1 text-xs rounded-md ${
              filter === 'all' 
                ? 'bg-gray-200 dark:bg-gray-600' 
                : 'hover:bg-gray-100 dark:hover:bg-gray-700'
            }`}
          >
            All ({serviceNames.length})
          </button>
          <button
            onClick={() => setFilter('closed')}
            className={`px-3 py-1 text-xs rounded-md ${
              filter === 'closed' 
                ? 'bg-green-100 text-green-800 dark:bg-green-900 dark:text-green-300' 
                : 'hover:bg-gray-100 dark:hover:bg-gray-700'
            }`}
          >
            Healthy ({counts[CIRCUIT_STATE.CLOSED] || 0})
          </button>
          <button
            onClick={() => setFilter('half-open')}
            className={`px-3 py-1 text-xs rounded-md ${
              filter === 'half-open' 
                ? 'bg-yellow-100 text-yellow-800 dark:bg-yellow-900 dark:text-yellow-300' 
                : 'hover:bg-gray-100 dark:hover:bg-gray-700'
            }`}
          >
            Testing ({counts[CIRCUIT_STATE.HALF_OPEN] || 0})
          </button>
          <button
            onClick={() => setFilter('open')}
            className={`px-3 py-1 text-xs rounded-md ${
              filter === 'open' 
                ? 'bg-red-100 text-red-800 dark:bg-red-900 dark:text-red-300' 
                : 'hover:bg-gray-100 dark:hover:bg-gray-700'
            }`}
          >
            Failed ({counts[CIRCUIT_STATE.OPEN] || 0})
          </button>
        </div>
      </div>
      
      <div className="mt-4 space-y-3">
        {filteredServices.length === 0 ? (
          <p className="text-gray-500 dark:text-gray-400">No services match the selected filter</p>
        ) : (
          filteredServices.map(serviceName => {
            const service = services[serviceName];
            const { state, failures, lastFailure } = service;
            
            // Format last failure time if available
            const lastFailureTime = lastFailure 
              ? new Date(lastFailure).toLocaleTimeString() 
              : null;
            
            return (
              <div 
                key={serviceName}
                className="flex flex-col sm:flex-row sm:items-center sm:justify-between p-3 bg-gray-50 dark:bg-gray-700 rounded-md"
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
                      {state === CIRCUIT_STATE.CLOSED ? 'Healthy' :
                       state === CIRCUIT_STATE.HALF_OPEN ? 'Testing' :
                       'Failed'} 
                      {failures > 0 && ` (${failures} failures)`}
                      {lastFailureTime && ` - Last failure: ${lastFailureTime}`}
                    </p>
                  </div>
                </div>
                
                {state !== CIRCUIT_STATE.CLOSED && (
                  <button
                    onClick={() => resetService(serviceName)}
                    className="mt-2 sm:mt-0 px-3 py-1.5 text-xs font-medium text-white bg-blue-600 rounded-md hover:bg-blue-700 transition-colors"
                  >
                    Reset
                  </button>
                )}
              </div>
            );
          })
        )}
      </div>
    </div>
  );
}
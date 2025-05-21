import React, { createContext, useContext, useState, useEffect } from 'react';

// Circuit breaker states
export const CIRCUIT_STATE = {
  CLOSED: 'closed',
  OPEN: 'open',
  HALF_OPEN: 'half-open'
};

// Default fallback values
const defaultConfig = {
  // Maximum number of failures before opening the circuit
  failureThreshold: 5,
  // Time in milliseconds to wait before transitioning to half-open
  resetTimeout: 30000,
  // Maximum number of requests to allow in half-open state
  halfOpenRequestLimit: 3
};

// Create context with default values
const CircuitBreakerContext = createContext({
  // Current state of the circuit
  state: CIRCUIT_STATE.CLOSED,
  // Service health status
  services: {},
  // Function to report a failure for a service
  reportFailure: () => {},
  // Function to report a success for a service
  reportSuccess: () => {},
  // Get the state of a specific service
  getServiceState: () => CIRCUIT_STATE.CLOSED,
  // Reset the circuit breaker for a service
  resetService: () => {}
});

export const useCircuitBreaker = () => useContext(CircuitBreakerContext);

export function CircuitBreakerProvider({ children, config = {} }) {
  // Merge default config with provided config
  const circuitConfig = { ...defaultConfig, ...config };
  
  // Track services and their states
  const [services, setServices] = useState({});
  
  // Report a failure for a service
  const reportFailure = (serviceName) => {
    setServices(prev => {
      const service = prev[serviceName] || { 
        failures: 0, 
        state: CIRCUIT_STATE.CLOSED,
        lastFailure: null
      };
      
      const failures = service.state === CIRCUIT_STATE.CLOSED ? 
        service.failures + 1 : service.failures;
      
      const now = Date.now();
      
      // Open circuit if threshold is reached
      const newState = failures >= circuitConfig.failureThreshold ?
        CIRCUIT_STATE.OPEN : service.state;
      
      return {
        ...prev,
        [serviceName]: {
          ...service,
          failures,
          state: newState,
          lastFailure: now
        }
      };
    });
  };
  
  // Report a success for a service
  const reportSuccess = (serviceName) => {
    setServices(prev => {
      const service = prev[serviceName];
      
      // If service doesn't exist or is already closed, do nothing
      if (!service || service.state === CIRCUIT_STATE.CLOSED) {
        return prev;
      }
      
      // If in half-open state, close the circuit on success
      if (service.state === CIRCUIT_STATE.HALF_OPEN) {
        return {
          ...prev,
          [serviceName]: {
            ...service,
            failures: 0,
            state: CIRCUIT_STATE.CLOSED,
            lastFailure: null
          }
        };
      }
      
      return prev;
    });
  };
  
  // Get the state of a specific service
  const getServiceState = (serviceName) => {
    return services[serviceName]?.state || CIRCUIT_STATE.CLOSED;
  };
  
  // Reset a service to closed state
  const resetService = (serviceName) => {
    setServices(prev => ({
      ...prev,
      [serviceName]: {
        failures: 0,
        state: CIRCUIT_STATE.CLOSED,
        lastFailure: null
      }
    }));
  };
  
  // Check for services that need to transition from open to half-open
  useEffect(() => {
    const checkCircuits = () => {
      const now = Date.now();
      
      setServices(prev => {
        const updated = { ...prev };
        let hasChanges = false;
        
        Object.keys(updated).forEach(serviceName => {
          const service = updated[serviceName];
          
          // If service is open and reset timeout has elapsed, transition to half-open
          if (
            service.state === CIRCUIT_STATE.OPEN && 
            service.lastFailure && 
            (now - service.lastFailure) >= circuitConfig.resetTimeout
          ) {
            updated[serviceName] = {
              ...service,
              state: CIRCUIT_STATE.HALF_OPEN,
              requestCount: 0
            };
            hasChanges = true;
          }
        });
        
        return hasChanges ? updated : prev;
      });
    };
    
    // Check circuit states periodically
    const intervalId = setInterval(checkCircuits, 1000);
    
    return () => clearInterval(intervalId);
  }, [circuitConfig.resetTimeout]);
  
  // Create context value
  const contextValue = {
    services,
    reportFailure,
    reportSuccess,
    getServiceState,
    resetService
  };
  
  return (
    <CircuitBreakerContext.Provider value={contextValue}>
      {children}
    </CircuitBreakerContext.Provider>
  );
}
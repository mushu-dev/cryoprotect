// Helper to check if we should use mock services
const shouldUseMockServices = () => {
  // Check for environment variable - server-side
  if (typeof process !== 'undefined' && process.env?.NEXT_PUBLIC_USE_MOCK_DATA === 'true') {
    return true;
  }
  
  // Check for URL param - client-side only
  if (typeof window !== 'undefined') {
    // Allow using mock data via URL parameter for testing
    const urlParams = new URLSearchParams(window.location.search);
    if (urlParams.get('mock') === 'true') {
      return true;
    }
    
    // If the API calls are failing, we can fallback to mock data
    try {
      const apiStatus = localStorage.getItem('api_connection_status');
      if (apiStatus === 'failed') {
        return true;
      }
    } catch (error) {
      // Ignore localStorage errors
    }
  }
  
  return false;
};

// Get the appropriate molecule service
export const getMoleculeService = () => {
  // In real implementation, return actual services based on environment
  return {
    getMolecules: async () => [],
    getMoleculeById: async (id) => null,
    searchMolecules: async (query) => []
  };
};

// Get the appropriate mixture service
export const getMixtureService = () => {
  // In real implementation, return actual services based on environment
  return {
    getMixtures: async () => [],
    getMixtureById: async (id) => null,
    searchMixtures: async (query) => []
  };
};

// Get the appropriate protocol service
export const getProtocolService = () => {
  // In real implementation, return actual services based on environment
  return {
    getProtocols: async () => [],
    getProtocolById: async (id) => null, 
    createProtocol: async (data) => ({}),
    updateProtocol: async (id, data) => ({}),
    deleteProtocol: async (id) => true
  };
};

// Get a context-based protocol service
export const getContextProtocolService = (apiClient) => {
  // In real implementation, return actual services based on environment
  return {
    getProtocols: async () => [],
    getProtocolById: async (id) => null,
    createProtocol: async (data) => ({}),
    updateProtocol: async (id, data) => ({}),
    deleteProtocol: async (id) => true
  };
};

// Mark API connection as failed - will cause fallback to mock data
export const markApiConnectionFailed = () => {
  if (typeof window !== 'undefined') {
    try {
      localStorage.setItem('api_connection_status', 'failed');
    } catch (error) {
      // Ignore localStorage errors
    }
  }
};

// Reset API connection status
export const resetApiConnectionStatus = () => {
  if (typeof window !== 'undefined') {
    try {
      localStorage.removeItem('api_connection_status');
    } catch (error) {
      // Ignore localStorage errors
    }
  }
};
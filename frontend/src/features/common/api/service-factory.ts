import { moleculeService } from '@/features/molecules/services/molecule-service';
import { mockMoleculeService } from '@/features/molecules/services/mock-molecule-service';
import { mixtureService } from '@/features/mixtures/services/mixture-service';
import { mockMixtureService } from '@/features/mixtures/services/mock-mixture-service';
import { ProtocolService, MockProtocolService } from '@/features/protocols/services/protocol-service';
import { ApiProtocolService } from '@/features/protocols/services/protocol-service-api';
import { ContextProtocolService } from '@/features/protocols/services/protocol-service-with-context';
import { useApi } from '@/contexts/api-context';

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
  return shouldUseMockServices() ? mockMoleculeService : moleculeService;
};

// Get the appropriate mixture service
export const getMixtureService = () => {
  return shouldUseMockServices() ? mockMixtureService : mixtureService;
};

// Get the appropriate protocol service
export const getProtocolService = (): ProtocolService => {
  return shouldUseMockServices() 
    ? new MockProtocolService() 
    : new ApiProtocolService();
};

// Get a context-based protocol service
export const getContextProtocolService = (apiClient: any): ProtocolService => {
  return shouldUseMockServices()
    ? new MockProtocolService()
    : new ContextProtocolService(apiClient);
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
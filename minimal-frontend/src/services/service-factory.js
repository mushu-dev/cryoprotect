/**
 * Service factory
 * Provides factory methods to get the appropriate service implementations
 * Supports fallback to mock services when API is unavailable
 */
import { moleculeService } from './molecules/molecule-service';
import { mixtureService } from './mixtures/mixture-service';
import { enhancedConvexMoleculeService } from './enhanced-convex-service';

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
  // Check for enhanced Convex service preference
  if (typeof window !== 'undefined') {
    const urlParams = new URLSearchParams(window.location.search);
    if (urlParams.get('enhanced') === 'true' || window.location.pathname.includes('enhanced')) {
      return enhancedConvexMoleculeService;
    }
  }
  
  // For now, we don't have mock services in minimal-frontend
  // This can be extended later to include mock services
  return moleculeService;
};

// Get the enhanced molecule service directly
export const getEnhancedMoleculeService = () => {
  return enhancedConvexMoleculeService;
};

// Get the appropriate mixture service
export const getMixtureService = () => {
  // For now, we don't have mock services in minimal-frontend
  // This can be extended later to include mock services
  return mixtureService;
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
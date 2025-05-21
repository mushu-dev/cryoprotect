/**
 * Resilient Protocols Hook
 * 
 * This hook provides access to the resilient protocol service,
 * offering enhanced error handling and offline capabilities.
 */

import { useState, useEffect, useCallback } from 'react';
import { ResilientProtocolService } from '../services/protocol-service-resilient';
import { Protocol, ProtocolListParams } from '../services/protocol-service';
import { UUID } from 'crypto';

// Create a singleton instance of the resilient service
const protocolService = new ResilientProtocolService();

/**
 * Hook for accessing protocols with enhanced resilience
 */
export function useResilientProtocols() {
  const [loading, setLoading] = useState<boolean>(false);
  const [error, setError] = useState<Error | null>(null);
  const [protocols, setProtocols] = useState<Protocol[]>([]);
  const [totalProtocols, setTotalProtocols] = useState<number>(0);
  const [connectionStatus, setConnectionStatus] = useState<'connected' | 'degraded' | 'offline'>('connected');

  // Monitor connection status
  useEffect(() => {
    // Check localStorage for API connection status (updated by the resilient service)
    const checkConnectionStatus = () => {
      if (typeof window !== 'undefined') {
        const status = localStorage.getItem('api_connection_status');
        if (status === 'failed') {
          setConnectionStatus('offline');
        }
      }
    };
    
    // Check initially
    checkConnectionStatus();
    
    // Set up an interval to check periodically
    const interval = setInterval(checkConnectionStatus, 10000);
    
    return () => clearInterval(interval);
  }, []);

  /**
   * Fetch protocols with optional filtering and pagination
   */
  const fetchProtocols = useCallback(async (params?: ProtocolListParams) => {
    setLoading(true);
    setError(null);
    
    try {
      const result = await protocolService.getProtocols(params);
      setProtocols(result.data);
      setTotalProtocols(result.total);
      
      // If we're here, we must be connected
      if (connectionStatus !== 'connected') {
        setConnectionStatus('connected');
      }
      
      return result;
    } catch (err) {
      const error = err instanceof Error ? err : new Error('An unknown error occurred');
      setError(error);
      
      // Update connection status based on error
      if (error.message.includes('Circuit breaker is open')) {
        setConnectionStatus('degraded');
      } else if (
        error.message.includes('Network Error') || 
        error.message.includes('ECONNREFUSED') ||
        error.message.includes('ENOTFOUND')
      ) {
        setConnectionStatus('offline');
      }
      
      return {
        data: [],
        total: 0,
        page: 1,
        per_page: 10,
        total_pages: 0
      };
    } finally {
      setLoading(false);
    }
  }, [connectionStatus]);

  /**
   * Get a single protocol by ID
   */
  const getProtocol = useCallback(async (id: UUID) => {
    setLoading(true);
    setError(null);
    
    try {
      const result = await protocolService.getProtocol(id);
      
      // If we're here, we must be connected
      if (connectionStatus !== 'connected') {
        setConnectionStatus('connected');
      }
      
      return result;
    } catch (err) {
      const error = err instanceof Error ? err : new Error('An unknown error occurred');
      setError(error);
      
      // Update connection status based on error
      if (error.message.includes('Circuit breaker is open')) {
        setConnectionStatus('degraded');
      } else if (
        error.message.includes('Network Error') || 
        error.message.includes('ECONNREFUSED') ||
        error.message.includes('ENOTFOUND')
      ) {
        setConnectionStatus('offline');
      }
      
      // Return null for failed protocol fetch
      return null;
    } finally {
      setLoading(false);
    }
  }, [connectionStatus]);

  /**
   * Create a new protocol
   */
  const createProtocol = useCallback(async (protocol: Omit<Protocol, 'id' | 'provenance'>) => {
    setLoading(true);
    setError(null);
    
    try {
      const result = await protocolService.createProtocol(protocol);
      
      // If we're here, we must be connected
      if (connectionStatus !== 'connected') {
        setConnectionStatus('connected');
      }
      
      return result;
    } catch (err) {
      const error = err instanceof Error ? err : new Error('An unknown error occurred');
      setError(error);
      
      // Update connection status based on error
      if (error.message.includes('Circuit breaker is open')) {
        setConnectionStatus('degraded');
      } else if (
        error.message.includes('Network Error') || 
        error.message.includes('ECONNREFUSED') ||
        error.message.includes('ENOTFOUND')
      ) {
        setConnectionStatus('offline');
      }
      
      // Re-throw for the caller to handle
      throw error;
    } finally {
      setLoading(false);
    }
  }, [connectionStatus]);

  /**
   * Update an existing protocol
   */
  const updateProtocol = useCallback(async (id: UUID, protocol: Partial<Protocol>) => {
    setLoading(true);
    setError(null);
    
    try {
      const result = await protocolService.updateProtocol(id, protocol);
      
      // If we're here, we must be connected
      if (connectionStatus !== 'connected') {
        setConnectionStatus('connected');
      }
      
      return result;
    } catch (err) {
      const error = err instanceof Error ? err : new Error('An unknown error occurred');
      setError(error);
      
      // Update connection status based on error
      if (error.message.includes('Circuit breaker is open')) {
        setConnectionStatus('degraded');
      } else if (
        error.message.includes('Network Error') || 
        error.message.includes('ECONNREFUSED') ||
        error.message.includes('ENOTFOUND')
      ) {
        setConnectionStatus('offline');
      }
      
      // Re-throw for the caller to handle
      throw error;
    } finally {
      setLoading(false);
    }
  }, [connectionStatus]);

  /**
   * Delete a protocol
   */
  const deleteProtocol = useCallback(async (id: UUID) => {
    setLoading(true);
    setError(null);
    
    try {
      await protocolService.deleteProtocol(id);
      
      // If we're here, we must be connected
      if (connectionStatus !== 'connected') {
        setConnectionStatus('connected');
      }
      
      return true;
    } catch (err) {
      const error = err instanceof Error ? err : new Error('An unknown error occurred');
      setError(error);
      
      // Update connection status based on error
      if (error.message.includes('Circuit breaker is open')) {
        setConnectionStatus('degraded');
      } else if (
        error.message.includes('Network Error') || 
        error.message.includes('ECONNREFUSED') ||
        error.message.includes('ENOTFOUND')
      ) {
        setConnectionStatus('offline');
      }
      
      // Return false for failed deletion
      return false;
    } finally {
      setLoading(false);
    }
  }, [connectionStatus]);

  /**
   * Get protocol templates
   */
  const getTemplates = useCallback(async (params?: ProtocolListParams) => {
    setLoading(true);
    setError(null);
    
    try {
      const result = await protocolService.getTemplates(params);
      
      // If we're here, we must be connected
      if (connectionStatus !== 'connected') {
        setConnectionStatus('connected');
      }
      
      return result;
    } catch (err) {
      const error = err instanceof Error ? err : new Error('An unknown error occurred');
      setError(error);
      
      // Update connection status based on error
      if (error.message.includes('Circuit breaker is open')) {
        setConnectionStatus('degraded');
      } else if (
        error.message.includes('Network Error') || 
        error.message.includes('ECONNREFUSED') ||
        error.message.includes('ENOTFOUND')
      ) {
        setConnectionStatus('offline');
      }
      
      return {
        data: [],
        total: 0,
        page: 1,
        per_page: 10,
        total_pages: 0
      };
    } finally {
      setLoading(false);
    }
  }, [connectionStatus]);

  /**
   * Create from template
   */
  const createFromTemplate = useCallback(async (templateId: UUID, customizations?: Partial<Protocol>) => {
    setLoading(true);
    setError(null);
    
    try {
      const result = await protocolService.createFromTemplate(templateId, customizations);
      
      // If we're here, we must be connected
      if (connectionStatus !== 'connected') {
        setConnectionStatus('connected');
      }
      
      return result;
    } catch (err) {
      const error = err instanceof Error ? err : new Error('An unknown error occurred');
      setError(error);
      
      // Update connection status based on error
      if (error.message.includes('Circuit breaker is open')) {
        setConnectionStatus('degraded');
      } else if (
        error.message.includes('Network Error') || 
        error.message.includes('ECONNREFUSED') ||
        error.message.includes('ENOTFOUND')
      ) {
        setConnectionStatus('offline');
      }
      
      // Re-throw for the caller to handle
      throw error;
    } finally {
      setLoading(false);
    }
  }, [connectionStatus]);

  // Return the hook's API
  return {
    // State
    protocols,
    totalProtocols,
    loading,
    error,
    connectionStatus,
    
    // Methods
    fetchProtocols,
    getProtocol,
    createProtocol,
    updateProtocol,
    deleteProtocol,
    getTemplates,
    createFromTemplate,
    
    // Direct access to the service for advanced usage
    service: protocolService,
  };
}
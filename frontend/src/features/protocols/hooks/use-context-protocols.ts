/**
 * Context Protocols Hook
 * 
 * This hook provides access to protocol operations using the API context
 * for enhanced resilience and app-wide connection status tracking.
 */

import { useState, useEffect, useCallback } from 'react';
import { useApi } from '@/contexts/api-context';
import { ContextProtocolService } from '../services/protocol-service-with-context';
import { Protocol, ProtocolListParams } from '../services/protocol-service';
import { UUID } from 'crypto';
import { getContextProtocolService } from '@/features/common/api/service-factory';

/**
 * Hook for accessing protocols with enhanced resilience via context
 */
export function useContextProtocols() {
  const { apiClient, connectionStatus, refreshConnection } = useApi();
  const [protocolService] = useState(() => getContextProtocolService(apiClient));
  const [loading, setLoading] = useState<boolean>(false);
  const [error, setError] = useState<Error | null>(null);
  const [protocols, setProtocols] = useState<Protocol[]>([]);
  const [totalProtocols, setTotalProtocols] = useState<number>(0);

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
      
      return result;
    } catch (err) {
      const error = err instanceof Error ? err : new Error('An unknown error occurred');
      setError(error);
      
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
  }, [protocolService]);

  /**
   * Get a single protocol by ID
   */
  const getProtocol = useCallback(async (id: UUID) => {
    setLoading(true);
    setError(null);
    
    try {
      const result = await protocolService.getProtocol(id);
      return result;
    } catch (err) {
      const error = err instanceof Error ? err : new Error('An unknown error occurred');
      setError(error);
      return null;
    } finally {
      setLoading(false);
    }
  }, [protocolService]);

  /**
   * Create a new protocol
   */
  const createProtocol = useCallback(async (protocol: Omit<Protocol, 'id' | 'provenance'>) => {
    setLoading(true);
    setError(null);
    
    try {
      const result = await protocolService.createProtocol(protocol);
      return result;
    } catch (err) {
      const error = err instanceof Error ? err : new Error('An unknown error occurred');
      setError(error);
      throw error;
    } finally {
      setLoading(false);
    }
  }, [protocolService]);

  /**
   * Update an existing protocol
   */
  const updateProtocol = useCallback(async (id: UUID, protocol: Partial<Protocol>) => {
    setLoading(true);
    setError(null);
    
    try {
      const result = await protocolService.updateProtocol(id, protocol);
      return result;
    } catch (err) {
      const error = err instanceof Error ? err : new Error('An unknown error occurred');
      setError(error);
      throw error;
    } finally {
      setLoading(false);
    }
  }, [protocolService]);

  /**
   * Delete a protocol
   */
  const deleteProtocol = useCallback(async (id: UUID) => {
    setLoading(true);
    setError(null);
    
    try {
      await protocolService.deleteProtocol(id);
      return true;
    } catch (err) {
      const error = err instanceof Error ? err : new Error('An unknown error occurred');
      setError(error);
      return false;
    } finally {
      setLoading(false);
    }
  }, [protocolService]);

  /**
   * Get protocol templates
   */
  const getTemplates = useCallback(async (params?: ProtocolListParams) => {
    setLoading(true);
    setError(null);
    
    try {
      const result = await protocolService.getTemplates(params);
      return result;
    } catch (err) {
      const error = err instanceof Error ? err : new Error('An unknown error occurred');
      setError(error);
      
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
  }, [protocolService]);

  /**
   * Create from template
   */
  const createFromTemplate = useCallback(async (templateId: UUID, customizations?: Partial<Protocol>) => {
    setLoading(true);
    setError(null);
    
    try {
      const result = await protocolService.createFromTemplate(templateId, customizations);
      return result;
    } catch (err) {
      const error = err instanceof Error ? err : new Error('An unknown error occurred');
      setError(error);
      throw error;
    } finally {
      setLoading(false);
    }
  }, [protocolService]);

  /**
   * Compare protocol versions
   */
  const compareProtocolVersions = useCallback(async (protocolId1: UUID, protocolId2: UUID) => {
    setLoading(true);
    setError(null);
    
    try {
      const result = await protocolService.compareProtocolVersions(protocolId1, protocolId2);
      return result;
    } catch (err) {
      const error = err instanceof Error ? err : new Error('An unknown error occurred');
      setError(error);
      return null;
    } finally {
      setLoading(false);
    }
  }, [protocolService]);

  /**
   * Validate a protocol
   */
  const validateProtocol = useCallback(async (protocol: Protocol) => {
    setLoading(true);
    setError(null);
    
    try {
      const result = await protocolService.validateProtocol(protocol);
      return result;
    } catch (err) {
      const error = err instanceof Error ? err : new Error('An unknown error occurred');
      setError(error);
      return {
        valid: false,
        errors: [{ message: error.message, severity: 'error' as const }],
        warnings: [],
        suggestions: [],
      };
    } finally {
      setLoading(false);
    }
  }, [protocolService]);

  /**
   * Export a protocol
   */
  const exportProtocol = useCallback(async (id: UUID, format: 'json' | 'yaml' | 'pdf' | 'human-readable') => {
    setLoading(true);
    setError(null);
    
    try {
      const result = await protocolService.exportProtocol(id, format);
      return result;
    } catch (err) {
      const error = err instanceof Error ? err : new Error('An unknown error occurred');
      setError(error);
      return null;
    } finally {
      setLoading(false);
    }
  }, [protocolService]);

  // Return the hook's API
  return {
    // State
    protocols,
    totalProtocols,
    loading,
    error,
    connectionStatus,
    
    // Connection management
    refreshConnection,
    
    // Methods
    fetchProtocols,
    getProtocol,
    createProtocol,
    updateProtocol,
    deleteProtocol,
    getTemplates,
    createFromTemplate,
    compareProtocolVersions,
    validateProtocol,
    exportProtocol,
    
    // Direct access to the service for advanced usage
    service: protocolService,
  };
}
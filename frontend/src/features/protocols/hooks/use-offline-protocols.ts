/**
 * Offline-Capable Protocols Hook
 * 
 * This hook provides access to protocols with full offline capabilities,
 * allowing users to work with protocols even when no connection is available.
 * It leverages the SyncService for data synchronization and offline storage.
 */

import { useState, useEffect, useCallback } from 'react';
import { UUID } from 'crypto';
import { OfflineProtocolService } from '../services/protocol-service-offline';
import { Protocol, ProtocolListParams } from '../services/protocol-service';
import { ResilientApiClient } from '@/services/resilient-api-client';
import { getSyncService } from '@/services/sync-service';
import { getOfflineStorage } from '@/services/offline-storage';

// Initialize dependencies
const offlineStorage = getOfflineStorage();
const apiClient = new ResilientApiClient({
  baseURL: '/api/v1',
  timeout: 15000,
  enableCache: true,
  retries: 3,
  logRequests: process.env.NODE_ENV === 'development',
});
const syncService = getSyncService(apiClient);

// Create a singleton instance of the offline service
const protocolService = new OfflineProtocolService(apiClient, syncService);

/**
 * Hook for accessing protocols with full offline capabilities
 */
export function useOfflineProtocols() {
  const [loading, setLoading] = useState<boolean>(false);
  const [error, setError] = useState<Error | null>(null);
  const [protocols, setProtocols] = useState<Protocol[]>([]);
  const [totalProtocols, setTotalProtocols] = useState<number>(0);
  const [connectionStatus, setConnectionStatus] = useState<'connected' | 'degraded' | 'offline'>('connected');
  const [syncStatus, setSyncStatus] = useState<{
    pendingChanges: number;
    lastSyncTime: number | null;
  }>({
    pendingChanges: 0,
    lastSyncTime: null,
  });

  // Monitor connection status and sync status
  useEffect(() => {
    // Check connection status
    const connectionStatusListener = apiClient.onConnectionStatusChange((status) => {
      setConnectionStatus(status);
    });

    // Get initial sync status
    const updateSyncStatus = async () => {
      try {
        const pendingChanges = await syncService.getPendingChangesCount();
        const lastSyncTime = await syncService.getLastSyncTime();
        
        setSyncStatus({
          pendingChanges,
          lastSyncTime,
        });
      } catch (error) {
        console.error('Failed to get sync status:', error);
      }
    };
    
    // Register for sync events
    const syncSuccessListener = syncService.onSyncSuccess(() => {
      updateSyncStatus();
    });
    
    // Handle browser online/offline events
    const handleOnline = () => {
      console.log('Browser online event detected');
      // When we come back online, try to sync
      if (navigator.onLine) {
        syncService.sync().catch(console.error);
      }
    };
    
    const handleOffline = () => {
      console.log('Browser offline event detected');
      setConnectionStatus('offline');
    };
    
    // Add event listeners
    window.addEventListener('online', handleOnline);
    window.addEventListener('offline', handleOffline);
    
    // Update initially
    updateSyncStatus();
    
    // Set up interval to update periodically
    const interval = setInterval(updateSyncStatus, 30000);
    
    return () => {
      clearInterval(interval);
      connectionStatusListener();
      syncSuccessListener();
      window.removeEventListener('online', handleOnline);
      window.removeEventListener('offline', handleOffline);
    };
  }, []);

  /**
   * Trigger manual sync
   */
  const synchronize = useCallback(async () => {
    setLoading(true);
    
    try {
      const result = await syncService.sync({ forceSync: true });
      
      // Update connection status based on sync result
      if (result.success) {
        setConnectionStatus('connected');
      } else if (result.failed > 0 && navigator.onLine) {
        setConnectionStatus('degraded');
      } else if (!navigator.onLine) {
        setConnectionStatus('offline');
      }
      
      // Update sync status
      const pendingChanges = await syncService.getPendingChangesCount();
      const lastSyncTime = await syncService.getLastSyncTime();
      
      setSyncStatus({
        pendingChanges,
        lastSyncTime,
      });
      
      return result;
    } catch (error) {
      console.error('Sync error:', error);
      return {
        success: false,
        synced: 0,
        failed: 1,
        conflicts: 0,
      };
    } finally {
      setLoading(false);
    }
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
  }, []);

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
  }, []);

  /**
   * Create a new protocol
   */
  const createProtocol = useCallback(async (protocol: Omit<Protocol, 'id' | 'provenance'>) => {
    setLoading(true);
    setError(null);
    
    try {
      const result = await protocolService.createProtocol(protocol);
      
      // Update sync status after creation
      const pendingChanges = await syncService.getPendingChangesCount();
      setSyncStatus(prev => ({
        ...prev,
        pendingChanges,
      }));
      
      return result;
    } catch (err) {
      const error = err instanceof Error ? err : new Error('An unknown error occurred');
      setError(error);
      throw error;
    } finally {
      setLoading(false);
    }
  }, []);

  /**
   * Update an existing protocol
   */
  const updateProtocol = useCallback(async (id: UUID, protocol: Partial<Protocol>) => {
    setLoading(true);
    setError(null);
    
    try {
      const result = await protocolService.updateProtocol(id, protocol);
      
      // Update sync status after update
      const pendingChanges = await syncService.getPendingChangesCount();
      setSyncStatus(prev => ({
        ...prev,
        pendingChanges,
      }));
      
      return result;
    } catch (err) {
      const error = err instanceof Error ? err : new Error('An unknown error occurred');
      setError(error);
      throw error;
    } finally {
      setLoading(false);
    }
  }, []);

  /**
   * Delete a protocol
   */
  const deleteProtocol = useCallback(async (id: UUID) => {
    setLoading(true);
    setError(null);
    
    try {
      await protocolService.deleteProtocol(id);
      
      // Update sync status after deletion
      const pendingChanges = await syncService.getPendingChangesCount();
      setSyncStatus(prev => ({
        ...prev,
        pendingChanges,
      }));
      
      return true;
    } catch (err) {
      const error = err instanceof Error ? err : new Error('An unknown error occurred');
      setError(error);
      return false;
    } finally {
      setLoading(false);
    }
  }, []);

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
  }, []);

  /**
   * Create from template
   */
  const createFromTemplate = useCallback(async (templateId: UUID, customizations?: Partial<Protocol>) => {
    setLoading(true);
    setError(null);
    
    try {
      const result = await protocolService.createFromTemplate(templateId, customizations);
      
      // Update sync status after creation
      const pendingChanges = await syncService.getPendingChangesCount();
      setSyncStatus(prev => ({
        ...prev,
        pendingChanges,
      }));
      
      return result;
    } catch (err) {
      const error = err instanceof Error ? err : new Error('An unknown error occurred');
      setError(error);
      throw error;
    } finally {
      setLoading(false);
    }
  }, []);

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
      throw error;
    } finally {
      setLoading(false);
    }
  }, []);

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
      throw error;
    } finally {
      setLoading(false);
    }
  }, []);

  /**
   * Check if device is online
   */
  const isOnline = useCallback(() => {
    return navigator.onLine;
  }, []);
  
  /**
   * Check if there are pending changes
   */
  const hasPendingChanges = useCallback(() => {
    return syncStatus.pendingChanges > 0;
  }, [syncStatus.pendingChanges]);

  /**
   * Get information about temporary (not yet synced) protocols
   */
  const getTemporaryProtocols = useCallback(async () => {
    try {
      const allProtocols = await syncService.getAllEntities<Protocol>('protocols');
      return allProtocols.filter(protocol => protocol._temp_id === true);
    } catch (error) {
      console.error('Failed to get temporary protocols:', error);
      return [];
    }
  }, []);

  // Return the hook's API
  return {
    // State
    protocols,
    totalProtocols,
    loading,
    error,
    connectionStatus,
    syncStatus,
    
    // Sync methods
    synchronize,
    isOnline,
    hasPendingChanges,
    getTemporaryProtocols,
    
    // Protocol methods
    fetchProtocols,
    getProtocol,
    createProtocol,
    updateProtocol,
    deleteProtocol,
    getTemplates,
    createFromTemplate,
    compareProtocolVersions,
    validateProtocol,
    
    // Direct access to the service for advanced usage
    service: protocolService,
  };
}
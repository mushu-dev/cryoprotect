/**
 * Hook for accessing offline mode status and functions
 * 
 * This hook provides access to offline mode detection, synchronization
 * and pending changes information throughout the application.
 */

import { useState, useEffect, useCallback } from 'react';
import { ResilientApiClient, ConnectionStatus } from '@/services/resilient-api-client';
import { getSyncService } from '@/services/sync-service';

interface UseOfflineModeOptions {
  apiClient?: ResilientApiClient;
}

/**
 * Hook for working with offline status
 */
export function useOfflineMode(options: UseOfflineModeOptions = {}) {
  const [isOnline, setIsOnline] = useState(
    typeof navigator !== 'undefined' ? navigator.onLine : true
  );
  
  const [connectionStatus, setConnectionStatus] = useState<ConnectionStatus>('connected');
  
  const [syncStatus, setSyncStatus] = useState<{
    pendingChanges: number;
    lastSyncTime: number | null;
    isSyncing: boolean;
  }>({
    pendingChanges: 0,
    lastSyncTime: null,
    isSyncing: false,
  });

  // Get or create API client instance
  const apiClient = options.apiClient || new ResilientApiClient({
    baseURL: '/api/v1',
    enableCache: true,
    retries: 3,
  });

  // Get sync service
  const syncService = getSyncService(apiClient);

  // Update sync status
  const updateSyncStatus = useCallback(async () => {
    try {
      const pendingChanges = await syncService.getPendingChangesCount();
      const lastSyncTime = await syncService.getLastSyncTime();
      
      setSyncStatus(prev => ({
        ...prev,
        pendingChanges,
        lastSyncTime,
      }));
    } catch (error) {
      console.error('Failed to get sync status:', error);
    }
  }, [syncService]);
  
  // Sync method
  const synchronize = useCallback(async (options = { forceSync: true }) => {
    if (syncStatus.isSyncing) return { success: false, message: 'Sync already in progress' };
    
    setSyncStatus(prev => ({ ...prev, isSyncing: true }));
    
    try {
      const result = await syncService.sync(options);
      await updateSyncStatus();
      return { success: true, result };
    } catch (error) {
      console.error('Sync error:', error);
      return { 
        success: false, 
        message: error instanceof Error ? error.message : 'Unknown error' 
      };
    } finally {
      setSyncStatus(prev => ({ ...prev, isSyncing: false }));
    }
  }, [syncService, syncStatus.isSyncing, updateSyncStatus]);
  
  // Format last sync time
  const formatLastSyncTime = useCallback(() => {
    if (!syncStatus.lastSyncTime) return 'Never';
    
    const now = new Date();
    const syncTime = new Date(syncStatus.lastSyncTime);
    const diffMs = now.getTime() - syncTime.getTime();
    
    // If less than a minute ago
    if (diffMs < 60 * 1000) {
      return 'Just now';
    }
    
    // If less than an hour ago
    if (diffMs < 60 * 60 * 1000) {
      const minutes = Math.floor(diffMs / (60 * 1000));
      return `${minutes} minute${minutes !== 1 ? 's' : ''} ago`;
    }
    
    // If less than a day ago
    if (diffMs < 24 * 60 * 60 * 1000) {
      const hours = Math.floor(diffMs / (60 * 60 * 1000));
      return `${hours} hour${hours !== 1 ? 's' : ''} ago`;
    }
    
    // Otherwise, show date
    return syncTime.toLocaleString();
  }, [syncStatus.lastSyncTime]);

  // Set up event listeners
  useEffect(() => {
    // Check online status
    setIsOnline(navigator.onLine);
    
    // Listen for connection status changes
    const unsubscribe = apiClient.onConnectionStatusChange((status) => {
      setConnectionStatus(status);
      setIsOnline(status !== 'offline');
    });
    
    // Listen for sync events
    const syncSuccessListener = syncService.onSyncSuccess(() => {
      updateSyncStatus();
    });
    
    // Handle browser online/offline events
    const handleOnline = () => {
      setIsOnline(true);
      synchronize({ forceSync: false }).catch(console.error);
    };
    
    const handleOffline = () => {
      setIsOnline(false);
    };
    
    // Add event listeners
    window.addEventListener('online', handleOnline);
    window.addEventListener('offline', handleOffline);
    
    // Update sync status initially
    updateSyncStatus();
    
    // Cleanup on unmount
    return () => {
      unsubscribe();
      syncSuccessListener();
      window.removeEventListener('online', handleOnline);
      window.removeEventListener('offline', handleOffline);
    };
  }, [apiClient, syncService, synchronize, updateSyncStatus]);

  return {
    // Basic status
    isOnline,
    connectionStatus,
    
    // Sync status
    hasPendingChanges: syncStatus.pendingChanges > 0,
    pendingChangesCount: syncStatus.pendingChanges,
    lastSyncTime: syncStatus.lastSyncTime,
    formattedLastSyncTime: formatLastSyncTime(),
    isSyncing: syncStatus.isSyncing,
    
    // Actions
    synchronize,
    
    // Access to services
    syncService,
    apiClient,
  };
}

export default useOfflineMode;
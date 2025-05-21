/**
 * Synchronization Service
 * 
 * This service manages data synchronization between the local offline storage
 * and the remote API, ensuring data consistency and handling conflicts.
 */

import { getOfflineStorage, StorageKey } from './offline-storage';
import { ResilientApiClient, ConnectionStatus } from './resilient-api-client';

// Constants
const SYNC_METADATA_KEY = 'sync_metadata';
const SYNC_INTERVAL = 5 * 60 * 1000; // 5 minutes

// Types
export interface SyncOptions {
  forceSync?: boolean;
  showNotifications?: boolean;
}

export interface SyncResult {
  success: boolean;
  synced: number;
  failed: number;
  conflicts: number;
  details?: string[];
}

export interface SyncMetadata {
  lastSyncTime: number;
  pendingChanges: Record<string, boolean>;
  syncErrors: Record<string, string>;
  deviceId: string;
}

export interface EntityMetadata {
  modifiedAt: number;
  modifiedBy: string;
  version: number;
  synchronized: boolean;
  syncAttempts: number;
  syncError?: string;
}

/**
 * Sync Service for handling data synchronization
 */
export class SyncService {
  private apiClient: ResilientApiClient;
  private offlineStorage = getOfflineStorage();
  private syncInterval: number | null = null;
  private deviceId: string;
  private syncInProgress = false;
  private onSyncStartHandlers: (() => void)[] = [];
  private onSyncCompleteHandlers: ((result: SyncResult) => void)[] = [];
  private onSyncErrorHandlers: ((error: Error) => void)[] = [];
  private entityPrefixes: Record<string, string> = {
    protocols: 'protocol_',
    experiments: 'experiment_',
    mixtures: 'mixture_',
    molecules: 'molecule_',
  };
  
  /**
   * Create a new sync service
   */
  constructor(apiClient: ResilientApiClient) {
    this.apiClient = apiClient;
    
    // Generate a device ID if not available
    const storedDeviceId = localStorage.getItem('device_id');
    if (storedDeviceId) {
      this.deviceId = storedDeviceId;
    } else {
      this.deviceId = this.generateDeviceId();
      localStorage.setItem('device_id', this.deviceId);
    }
    
    // Set up listeners for connectivity changes
    apiClient.onConnectionStatusChange(this.handleConnectionStatusChange.bind(this));
    
    // Set up listeners for offline storage sync events
    this.offlineStorage.onSyncSuccess(this.handleSyncSuccess.bind(this));
    this.offlineStorage.onSyncFailure(this.handleSyncFailure.bind(this));
  }
  
  /**
   * Start automatic synchronization
   */
  startAutoSync(): void {
    // Stop any existing sync interval
    this.stopAutoSync();
    
    // Start a new sync interval
    this.syncInterval = window.setInterval(
      () => this.sync({ showNotifications: false }),
      SYNC_INTERVAL
    );
    
    // Add event listeners for online/offline events
    window.addEventListener('online', this.handleOnline.bind(this));
    
    console.log('Automatic synchronization started');
  }
  
  /**
   * Stop automatic synchronization
   */
  stopAutoSync(): void {
    if (this.syncInterval !== null) {
      clearInterval(this.syncInterval);
      this.syncInterval = null;
    }
    
    // Remove event listeners
    window.removeEventListener('online', this.handleOnline);
    
    console.log('Automatic synchronization stopped');
  }
  
  /**
   * Synchronize data with the server
   */
  async sync(options: SyncOptions = {}): Promise<SyncResult> {
    // Skip if sync is already in progress
    if (this.syncInProgress) {
      return {
        success: false,
        synced: 0,
        failed: 0,
        conflicts: 0,
        details: ['Sync already in progress'],
      };
    }
    
    // Skip if offline and not forced
    if (!navigator.onLine && !options.forceSync) {
      return {
        success: false,
        synced: 0,
        failed: 0,
        conflicts: 0,
        details: ['Device is offline'],
      };
    }
    
    this.syncInProgress = true;
    
    // Notify sync start handlers
    this.notifySyncStart();
    
    try {
      // Get sync metadata
      const metadata = await this.getSyncMetadata();
      
      // Get all entities that need syncing
      const pendingEntities = await this.getPendingEntities();
      
      if (pendingEntities.length === 0 && !options.forceSync) {
        // Nothing to sync
        this.syncInProgress = false;
        
        const result: SyncResult = {
          success: true,
          synced: 0,
          failed: 0,
          conflicts: 0,
          details: ['No changes to sync'],
        };
        
        this.notifySyncComplete(result);
        return result;
      }
      
      console.log(`Syncing ${pendingEntities.length} entities with the server`);
      
      // Process each entity
      let syncedCount = 0;
      let failedCount = 0;
      let conflictCount = 0;
      const details: string[] = [];
      
      for (const entityKey of pendingEntities) {
        try {
          const syncResult = await this.syncEntity(entityKey);
          
          if (syncResult.success) {
            syncedCount++;
            details.push(`Synced: ${entityKey}`);
            
            // Mark as synced in metadata
            delete metadata.pendingChanges[entityKey];
            delete metadata.syncErrors[entityKey];
          } else if (syncResult.conflict) {
            conflictCount++;
            details.push(`Conflict: ${entityKey} - ${syncResult.error || 'Version conflict'}`);
            
            // Mark as conflict in metadata
            metadata.syncErrors[entityKey] = syncResult.error || 'Version conflict';
          } else {
            failedCount++;
            details.push(`Failed: ${entityKey} - ${syncResult.error || 'Unknown error'}`);
            
            // Mark as failed in metadata
            metadata.syncErrors[entityKey] = syncResult.error || 'Unknown error';
          }
        } catch (error) {
          failedCount++;
          const errorMessage = error instanceof Error ? error.message : String(error);
          details.push(`Error: ${entityKey} - ${errorMessage}`);
          
          // Mark as failed in metadata
          metadata.syncErrors[entityKey] = errorMessage;
        }
      }
      
      // If force sync, try to download all data too
      if (options.forceSync) {
        await this.downloadAllEntities();
      }
      
      // Update sync metadata
      metadata.lastSyncTime = Date.now();
      await this.setSyncMetadata(metadata);
      
      // Process queued requests
      await this.offlineStorage.processRequestQueue();
      
      const result: SyncResult = {
        success: syncedCount > 0 && failedCount === 0,
        synced: syncedCount,
        failed: failedCount,
        conflicts: conflictCount,
        details,
      };
      
      // Notify sync complete handlers
      this.notifySyncComplete(result);
      
      return result;
    } catch (error) {
      // Handle sync error
      const errorMessage = error instanceof Error ? error.message : String(error);
      console.error('Sync error:', errorMessage);
      
      const result: SyncResult = {
        success: false,
        synced: 0,
        failed: 1,
        conflicts: 0,
        details: [`Sync error: ${errorMessage}`],
      };
      
      // Notify sync error handlers
      this.notifySyncError(error instanceof Error ? error : new Error(errorMessage));
      
      return result;
    } finally {
      this.syncInProgress = false;
    }
  }
  
  /**
   * Sync a specific entity
   */
  private async syncEntity(entityKey: string): Promise<{
    success: boolean;
    conflict: boolean;
    error?: string;
  }> {
    try {
      // Get entity data from offline storage
      const entityData = await this.offlineStorage.get(entityKey);
      
      if (!entityData) {
        return { success: false, conflict: false, error: 'Entity not found in offline storage' };
      }
      
      // Extract entity metadata
      const entityMetadata = entityData._metadata as EntityMetadata;
      
      if (!entityMetadata) {
        return { success: false, conflict: false, error: 'Entity has no metadata' };
      }
      
      // If already synchronized, skip
      if (entityMetadata.synchronized) {
        return { success: true, conflict: false };
      }
      
      // Determine entity type and ID
      const { entityType, entityId } = this.parseEntityKey(entityKey);
      
      if (!entityType || !entityId) {
        return { success: false, conflict: false, error: 'Invalid entity key format' };
      }
      
      // Prepare data for API
      const apiData = { ...entityData };
      delete apiData._metadata;
      
      // Include metadata in headers
      const headers: Record<string, string> = {
        'X-Client-Modified-At': entityMetadata.modifiedAt.toString(),
        'X-Client-Modified-By': entityMetadata.modifiedBy,
        'X-Client-Version': entityMetadata.version.toString(),
        'X-Device-Id': this.deviceId,
      };
      
      // Make API request
      try {
        const response = await this.apiClient.put(
          `/${entityType}/${entityId}`,
          apiData,
          { headers }
        );
        
        // Update entity metadata
        entityMetadata.synchronized = true;
        entityMetadata.syncAttempts += 1;
        delete entityMetadata.syncError;
        
        // Save updated entity with metadata
        await this.offlineStorage.set(entityKey, {
          ...apiData,
          _metadata: entityMetadata,
        });
        
        return { success: true, conflict: false };
      } catch (error) {
        // Handle API error
        if (error instanceof Error && error.message.includes('409')) {
          // Conflict - server has a newer version
          
          // Try to get the server version
          try {
            const serverData = await this.apiClient.get(`/${entityType}/${entityId}`);
            
            // Resolve conflict (for now, just mark as conflict)
            entityMetadata.syncError = 'Version conflict';
            entityMetadata.syncAttempts += 1;
            
            // Save updated entity with metadata
            await this.offlineStorage.set(entityKey, {
              ...entityData,
              _metadata: entityMetadata,
            });
            
            return { success: false, conflict: true, error: 'Version conflict' };
          } catch (getError) {
            // Failed to get server version
            entityMetadata.syncError = 'Could not retrieve server version';
            entityMetadata.syncAttempts += 1;
            
            // Save updated entity with metadata
            await this.offlineStorage.set(entityKey, {
              ...entityData,
              _metadata: entityMetadata,
            });
            
            return { success: false, conflict: true, error: 'Could not retrieve server version' };
          }
        } else {
          // Other API error
          const errorMessage = error instanceof Error ? error.message : String(error);
          
          // Update entity metadata
          entityMetadata.syncError = errorMessage;
          entityMetadata.syncAttempts += 1;
          
          // Save updated entity with metadata
          await this.offlineStorage.set(entityKey, {
            ...entityData,
            _metadata: entityMetadata,
          });
          
          return { success: false, conflict: false, error: errorMessage };
        }
      }
    } catch (error) {
      // Handle any other errors
      const errorMessage = error instanceof Error ? error.message : String(error);
      return { success: false, conflict: false, error: errorMessage };
    }
  }
  
  /**
   * Download all entities from the server
   */
  private async downloadAllEntities(): Promise<void> {
    try {
      // For each entity type, download all entities
      for (const [entityType, prefix] of Object.entries(this.entityPrefixes)) {
        try {
          // Get all entities of this type
          const response = await this.apiClient.get(`/${entityType}`);
          
          if (Array.isArray(response.data)) {
            console.log(`Downloaded ${response.data.length} ${entityType} from server`);
            
            // Process each entity
            for (const entity of response.data) {
              if (entity && entity.id) {
                const entityKey = `${prefix}${entity.id}`;
                
                // Check if we already have this entity
                const existingEntity = await this.offlineStorage.get(entityKey);
                
                if (existingEntity) {
                  const existingMetadata = existingEntity._metadata as EntityMetadata;
                  
                  // Only replace if server version is newer or we don't have metadata
                  if (!existingMetadata || 
                      (entity.updated_at && new Date(entity.updated_at).getTime() > existingMetadata.modifiedAt)) {
                    // Create metadata for server entity
                    const serverMetadata: EntityMetadata = {
                      modifiedAt: entity.updated_at ? new Date(entity.updated_at).getTime() : Date.now(),
                      modifiedBy: entity.updated_by || 'server',
                      version: entity.version || 1,
                      synchronized: true,
                      syncAttempts: 0,
                    };
                    
                    // Save server entity with metadata
                    await this.offlineStorage.set(entityKey, {
                      ...entity,
                      _metadata: serverMetadata,
                    });
                  }
                } else {
                  // Create metadata for server entity
                  const serverMetadata: EntityMetadata = {
                    modifiedAt: entity.updated_at ? new Date(entity.updated_at).getTime() : Date.now(),
                    modifiedBy: entity.updated_by || 'server',
                    version: entity.version || 1,
                    synchronized: true,
                    syncAttempts: 0,
                  };
                  
                  // Save server entity with metadata
                  await this.offlineStorage.set(entityKey, {
                    ...entity,
                    _metadata: serverMetadata,
                  });
                }
              }
            }
          }
        } catch (error) {
          console.error(`Failed to download ${entityType}:`, error);
        }
      }
    } catch (error) {
      console.error('Failed to download entities:', error);
    }
  }
  
  /**
   * Save an entity to offline storage and mark it for sync
   */
  async saveEntity(entityType: string, entityId: string, data: any): Promise<void> {
    const prefix = this.entityPrefixes[entityType];
    
    if (!prefix) {
      throw new Error(`Unknown entity type: ${entityType}`);
    }
    
    const entityKey = `${prefix}${entityId}`;
    
    // Get existing entity if available
    const existingEntity = await this.offlineStorage.get(entityKey);
    const existingMetadata = existingEntity?._metadata as EntityMetadata | undefined;
    
    // Create or update metadata
    const metadata: EntityMetadata = {
      modifiedAt: Date.now(),
      modifiedBy: this.deviceId,
      version: existingMetadata ? existingMetadata.version + 1 : 1,
      synchronized: false,
      syncAttempts: 0,
    };
    
    // Save entity with metadata
    await this.offlineStorage.set(entityKey, {
      ...data,
      _metadata: metadata,
    });
    
    // Mark as pending in sync metadata
    const syncMetadata = await this.getSyncMetadata();
    syncMetadata.pendingChanges[entityKey] = true;
    await this.setSyncMetadata(syncMetadata);
    
    // If we're online, try to sync immediately
    if (navigator.onLine) {
      this.sync({ showNotifications: false }).catch(error => {
        console.error('Failed to sync after save:', error);
      });
    }
  }
  
  /**
   * Get an entity from offline storage
   */
  async getEntity<T = any>(entityType: string, entityId: string): Promise<T | null> {
    const prefix = this.entityPrefixes[entityType];
    
    if (!prefix) {
      throw new Error(`Unknown entity type: ${entityType}`);
    }
    
    const entityKey = `${prefix}${entityId}`;
    
    // Get entity from offline storage
    const entity = await this.offlineStorage.get<T & { _metadata?: EntityMetadata }>(entityKey);
    
    if (!entity) {
      return null;
    }
    
    // Clone entity and remove metadata
    const result = { ...entity };
    delete result._metadata;
    
    return result as T;
  }
  
  /**
   * Get all entities of a specific type
   */
  async getAllEntities<T = any>(entityType: string): Promise<T[]> {
    const prefix = this.entityPrefixes[entityType];
    
    if (!prefix) {
      throw new Error(`Unknown entity type: ${entityType}`);
    }
    
    // Get all keys with this prefix
    const keys = await this.offlineStorage.getKeysByPrefix(prefix);
    
    // Get all entities
    const entities: T[] = [];
    
    for (const key of keys) {
      const entity = await this.offlineStorage.get<T & { _metadata?: EntityMetadata }>(key);
      
      if (entity) {
        // Clone entity and remove metadata
        const result = { ...entity };
        delete result._metadata;
        
        entities.push(result as T);
      }
    }
    
    return entities;
  }
  
  /**
   * Delete an entity
   */
  async deleteEntity(entityType: string, entityId: string): Promise<void> {
    const prefix = this.entityPrefixes[entityType];
    
    if (!prefix) {
      throw new Error(`Unknown entity type: ${entityType}`);
    }
    
    const entityKey = `${prefix}${entityId}`;
    
    // Delete from offline storage
    await this.offlineStorage.remove(entityKey);
    
    // Queue delete request
    await this.offlineStorage.queueRequest(
      `/${entityType}/${entityId}`,
      'DELETE',
      undefined,
      { 'X-Device-Id': this.deviceId }
    );
    
    // Update sync metadata
    const metadata = await this.getSyncMetadata();
    delete metadata.pendingChanges[entityKey];
    delete metadata.syncErrors[entityKey];
    await this.setSyncMetadata(metadata);
    
    // If we're online, try to sync immediately
    if (navigator.onLine) {
      this.offlineStorage.processRequestQueue().catch(error => {
        console.error('Failed to process request queue after delete:', error);
      });
    }
  }
  
  /**
   * Register a callback for sync start
   */
  onSyncStart(callback: () => void): () => void {
    this.onSyncStartHandlers.push(callback);
    
    // Return a function to unregister the callback
    return () => {
      const index = this.onSyncStartHandlers.indexOf(callback);
      if (index !== -1) {
        this.onSyncStartHandlers.splice(index, 1);
      }
    };
  }
  
  /**
   * Register a callback for sync completion
   */
  onSyncComplete(callback: (result: SyncResult) => void): () => void {
    this.onSyncCompleteHandlers.push(callback);
    
    // Return a function to unregister the callback
    return () => {
      const index = this.onSyncCompleteHandlers.indexOf(callback);
      if (index !== -1) {
        this.onSyncCompleteHandlers.splice(index, 1);
      }
    };
  }
  
  /**
   * Register a callback for sync errors
   */
  onSyncError(callback: (error: Error) => void): () => void {
    this.onSyncErrorHandlers.push(callback);
    
    // Return a function to unregister the callback
    return () => {
      const index = this.onSyncErrorHandlers.indexOf(callback);
      if (index !== -1) {
        this.onSyncErrorHandlers.splice(index, 1);
      }
    };
  }
  
  /**
   * Handle connection status changes
   */
  private handleConnectionStatusChange(status: ConnectionStatus): void {
    if (status === 'connected') {
      // We're back online, try to sync
      this.sync({ showNotifications: false }).catch(error => {
        console.error('Failed to sync after connection status change:', error);
      });
    }
  }
  
  /**
   * Handle device coming back online
   */
  private handleOnline(): void {
    console.log('Device is back online, syncing...');
    
    // Sync with a slight delay to ensure the connection is stable
    setTimeout(() => {
      this.sync({ showNotifications: true }).catch(error => {
        console.error('Failed to sync after coming back online:', error);
      });
    }, 2000);
  }
  
  /**
   * Handle successful sync from offline storage
   */
  private handleSyncSuccess(keys: StorageKey[]): void {
    console.log(`${keys.length} items synced successfully from request queue`);
  }
  
  /**
   * Handle sync failure from offline storage
   */
  private handleSyncFailure(error: Error): void {
    console.error('Failed to sync from request queue:', error);
  }
  
  /**
   * Notify sync start handlers
   */
  private notifySyncStart(): void {
    for (const handler of this.onSyncStartHandlers) {
      try {
        handler();
      } catch (error) {
        console.error('Error in sync start handler:', error);
      }
    }
  }
  
  /**
   * Notify sync complete handlers
   */
  private notifySyncComplete(result: SyncResult): void {
    for (const handler of this.onSyncCompleteHandlers) {
      try {
        handler(result);
      } catch (error) {
        console.error('Error in sync complete handler:', error);
      }
    }
  }
  
  /**
   * Notify sync error handlers
   */
  private notifySyncError(error: Error): void {
    for (const handler of this.onSyncErrorHandlers) {
      try {
        handler(error);
      } catch (handlerError) {
        console.error('Error in sync error handler:', handlerError);
      }
    }
  }
  
  /**
   * Generate a unique device ID
   */
  private generateDeviceId(): string {
    const navigator = window.navigator;
    const screen = window.screen;
    let deviceId = '';
    
    // Use available device information to generate a unique ID
    const components = [
      navigator.userAgent,
      screen.width.toString(),
      screen.height.toString(),
      navigator.language,
      new Date().getTimezoneOffset().toString(),
    ];
    
    // Generate a hash of the components
    deviceId = components.join('|');
    
    // Use simple hashing algorithm
    let hash = 0;
    for (let i = 0; i < deviceId.length; i++) {
      const char = deviceId.charCodeAt(i);
      hash = ((hash << 5) - hash) + char;
      hash = hash & hash; // Convert to 32bit integer
    }
    
    // Add a timestamp for uniqueness
    const timestamp = Date.now().toString(36);
    
    // Format the device ID
    return `device_${hash.toString(36)}_${timestamp}`;
  }
  
  /**
   * Get sync metadata from offline storage
   */
  private async getSyncMetadata(): Promise<SyncMetadata> {
    const metadata = await this.offlineStorage.get<SyncMetadata>(SYNC_METADATA_KEY);
    
    if (metadata) {
      return metadata;
    }
    
    // Create new metadata
    return {
      lastSyncTime: 0,
      pendingChanges: {},
      syncErrors: {},
      deviceId: this.deviceId,
    };
  }
  
  /**
   * Set sync metadata in offline storage
   */
  private async setSyncMetadata(metadata: SyncMetadata): Promise<void> {
    await this.offlineStorage.set(SYNC_METADATA_KEY, metadata);
  }
  
  /**
   * Get all entities that need syncing
   */
  private async getPendingEntities(): Promise<string[]> {
    const metadata = await this.getSyncMetadata();
    return Object.keys(metadata.pendingChanges);
  }
  
  /**
   * Parse an entity key into type and ID
   */
  private parseEntityKey(entityKey: string): { entityType: string; entityId: string } {
    for (const [entityType, prefix] of Object.entries(this.entityPrefixes)) {
      if (entityKey.startsWith(prefix)) {
        return {
          entityType,
          entityId: entityKey.substring(prefix.length),
        };
      }
    }
    
    return { entityType: '', entityId: '' };
  }
  
  /**
   * Get all sync errors
   */
  async getSyncErrors(): Promise<Record<string, string>> {
    const metadata = await this.getSyncMetadata();
    return metadata.syncErrors;
  }
  
  /**
   * Get the last sync time
   */
  async getLastSyncTime(): Promise<number> {
    const metadata = await this.getSyncMetadata();
    return metadata.lastSyncTime;
  }
  
  /**
   * Get the number of pending changes
   */
  async getPendingChangesCount(): Promise<number> {
    const metadata = await this.getSyncMetadata();
    return Object.keys(metadata.pendingChanges).length;
  }
  
  /**
   * Get the device ID
   */
  getDeviceId(): string {
    return this.deviceId;
  }
}

// Create a singleton instance
let syncServiceInstance: SyncService | null = null;

/**
 * Get the sync service instance
 */
export function getSyncService(apiClient: ResilientApiClient): SyncService {
  if (!syncServiceInstance) {
    syncServiceInstance = new SyncService(apiClient);
  }
  
  return syncServiceInstance;
}
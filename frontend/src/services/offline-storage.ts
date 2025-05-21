/**
 * Offline Storage Service
 * 
 * This service provides a robust offline storage solution with:
 * - Structured data storage using IndexedDB
 * - Fallback to localStorage when IndexedDB is not available
 * - Request queue for offline mutations
 * - Data expiration policies
 * - Storage limit management
 * - Data encryption for sensitive information
 */

import { v4 as uuidv4 } from 'uuid';

// Types
export type StorageKey = string;
export type StorageValue = any;
export type ExpiryTime = number | null; // null means no expiry
export type StorageArea = 'app' | 'user' | 'cache' | 'sync';

export interface StorageOptions {
  expiry?: ExpiryTime;
  encrypt?: boolean;
  area?: StorageArea;
}

export interface StorageEntry {
  key: StorageKey;
  value: StorageValue;
  timestamp: number;
  expiry: ExpiryTime;
  area: StorageArea;
  encrypted: boolean;
}

export interface QueuedRequest {
  id: string;
  url: string;
  method: string;
  body?: any;
  headers?: Record<string, string>;
  timestamp: number;
  retryCount: number;
  status: 'pending' | 'processing' | 'failed';
  errorMessage?: string;
}

// Constants
const DB_NAME = 'cryoprotect_offline_storage';
const DB_VERSION = 1;
const STORAGE_STORE = 'storage';
const REQUEST_QUEUE_STORE = 'request_queue';
const METADATA_STORE = 'metadata';
const DEFAULT_EXPIRY = 7 * 24 * 60 * 60 * 1000; // 7 days
const MAX_STORAGE_SIZE = 50 * 1024 * 1024; // 50MB
const ENCRYPTION_KEY_STORAGE = 'cryoprotect_encryption_key';
const LOCALSTORAGE_PREFIX = 'cryoprotect_offline_';

/**
 * Offline Storage Service
 */
export class OfflineStorage {
  private db: IDBDatabase | null = null;
  private dbPromise: Promise<IDBDatabase> | null = null;
  private encryptionKey: CryptoKey | null = null;
  private ready = false;
  private useLocalStorageFallback = false;
  private queueProcessingInterval: number | null = null;
  private onSyncSuccessHandlers: ((keys: StorageKey[]) => void)[] = [];
  private onSyncFailureHandlers: ((error: Error) => void)[] = [];
  private onStorageChangeHandlers: ((key: StorageKey, newValue: StorageValue | null) => void)[] = [];
  
  /**
   * Initialize the offline storage
   */
  async init(): Promise<void> {
    try {
      // Try to initialize IndexedDB
      this.db = await this.initIndexedDB();
      this.ready = true;
      
      // Try to initialize encryption
      await this.initEncryption();
      
      console.log('Offline storage initialized with IndexedDB');
    } catch (error) {
      console.warn('Failed to initialize IndexedDB, falling back to localStorage', error);
      this.useLocalStorageFallback = true;
      this.ready = true;
    }
    
    // Initialize request queue processing
    this.startQueueProcessing();
    
    // Clean up expired items
    this.cleanupExpiredItems();
    
    return Promise.resolve();
  }
  
  /**
   * Initialize IndexedDB
   */
  private initIndexedDB(): Promise<IDBDatabase> {
    if (this.dbPromise) {
      return this.dbPromise;
    }
    
    this.dbPromise = new Promise((resolve, reject) => {
      if (!window.indexedDB) {
        reject(new Error('IndexedDB is not supported'));
        return;
      }
      
      const request = indexedDB.open(DB_NAME, DB_VERSION);
      
      request.onerror = () => {
        reject(new Error('Failed to open IndexedDB'));
      };
      
      request.onsuccess = (event) => {
        const db = (event.target as IDBOpenDBRequest).result;
        resolve(db);
      };
      
      request.onupgradeneeded = (event) => {
        const db = (event.target as IDBOpenDBRequest).result;
        
        // Create storage object store
        if (!db.objectStoreNames.contains(STORAGE_STORE)) {
          const store = db.createObjectStore(STORAGE_STORE, { keyPath: 'key' });
          store.createIndex('area', 'area', { unique: false });
          store.createIndex('timestamp', 'timestamp', { unique: false });
          store.createIndex('expiry', 'expiry', { unique: false });
        }
        
        // Create request queue object store
        if (!db.objectStoreNames.contains(REQUEST_QUEUE_STORE)) {
          const store = db.createObjectStore(REQUEST_QUEUE_STORE, { keyPath: 'id' });
          store.createIndex('status', 'status', { unique: false });
          store.createIndex('timestamp', 'timestamp', { unique: false });
        }
        
        // Create metadata object store
        if (!db.objectStoreNames.contains(METADATA_STORE)) {
          db.createObjectStore(METADATA_STORE, { keyPath: 'key' });
        }
      };
    });
    
    return this.dbPromise;
  }
  
  /**
   * Initialize encryption for sensitive data
   */
  private async initEncryption(): Promise<void> {
    try {
      // Check if crypto is available
      if (!window.crypto || !window.crypto.subtle) {
        console.warn('Web Crypto API is not supported, encryption will be disabled');
        return;
      }
      
      // Check if we already have an encryption key
      let storedKey = localStorage.getItem(ENCRYPTION_KEY_STORAGE);
      
      if (storedKey) {
        // Import existing key
        const keyData = JSON.parse(storedKey);
        const keyBuffer = new Uint8Array(keyData).buffer;
        
        this.encryptionKey = await window.crypto.subtle.importKey(
          'raw',
          keyBuffer,
          { name: 'AES-GCM' },
          false,
          ['encrypt', 'decrypt']
        );
      } else {
        // Generate a new key
        this.encryptionKey = await window.crypto.subtle.generateKey(
          {
            name: 'AES-GCM',
            length: 256,
          },
          true,
          ['encrypt', 'decrypt']
        );
        
        // Export and store the key
        const keyBuffer = await window.crypto.subtle.exportKey('raw', this.encryptionKey);
        const keyArray = Array.from(new Uint8Array(keyBuffer));
        localStorage.setItem(ENCRYPTION_KEY_STORAGE, JSON.stringify(keyArray));
      }
      
      console.log('Encryption initialized successfully');
    } catch (error) {
      console.warn('Failed to initialize encryption', error);
      this.encryptionKey = null;
    }
  }
  
  /**
   * Encrypt sensitive data
   */
  private async encrypt(data: any): Promise<string> {
    if (!this.encryptionKey) {
      console.warn('Encryption is not available, storing data unencrypted');
      return JSON.stringify(data);
    }
    
    try {
      // Generate a random initialization vector
      const iv = window.crypto.getRandomValues(new Uint8Array(12));
      
      // Convert data to string and then to ArrayBuffer
      const dataString = JSON.stringify(data);
      const encoder = new TextEncoder();
      const dataBuffer = encoder.encode(dataString);
      
      // Encrypt the data
      const encryptedBuffer = await window.crypto.subtle.encrypt(
        {
          name: 'AES-GCM',
          iv,
        },
        this.encryptionKey,
        dataBuffer
      );
      
      // Combine IV and encrypted data
      const encryptedArray = new Uint8Array(encryptedBuffer);
      const result = new Uint8Array(iv.length + encryptedArray.length);
      result.set(iv);
      result.set(encryptedArray, iv.length);
      
      // Convert to Base64 for storage
      return btoa(String.fromCharCode.apply(null, Array.from(result)));
    } catch (error) {
      console.error('Encryption failed', error);
      return JSON.stringify(data);
    }
  }
  
  /**
   * Decrypt sensitive data
   */
  private async decrypt(encryptedData: string): Promise<any> {
    if (!this.encryptionKey) {
      try {
        return JSON.parse(encryptedData);
      } catch (e) {
        return encryptedData;
      }
    }
    
    try {
      // Convert from Base64 to Uint8Array
      const buffer = Uint8Array.from(atob(encryptedData), c => c.charCodeAt(0));
      
      // Extract IV and encrypted data
      const iv = buffer.slice(0, 12);
      const data = buffer.slice(12);
      
      // Decrypt the data
      const decryptedBuffer = await window.crypto.subtle.decrypt(
        {
          name: 'AES-GCM',
          iv,
        },
        this.encryptionKey,
        data
      );
      
      // Convert ArrayBuffer to string and parse JSON
      const decoder = new TextDecoder();
      const decryptedString = decoder.decode(decryptedBuffer);
      return JSON.parse(decryptedString);
    } catch (error) {
      console.error('Decryption failed', error);
      try {
        return JSON.parse(encryptedData);
      } catch (e) {
        return encryptedData;
      }
    }
  }
  
  /**
   * Store an item in offline storage
   */
  async set(key: StorageKey, value: StorageValue, options: StorageOptions = {}): Promise<void> {
    if (!this.ready) {
      await this.init();
    }
    
    const timestamp = Date.now();
    const expiry = options.expiry !== undefined ? options.expiry : DEFAULT_EXPIRY;
    const area = options.area || 'app';
    const shouldEncrypt = options.encrypt || false;
    
    const entry: StorageEntry = {
      key,
      value: shouldEncrypt ? await this.encrypt(value) : value,
      timestamp,
      expiry,
      area,
      encrypted: shouldEncrypt,
    };
    
    if (this.useLocalStorageFallback) {
      // Use localStorage fallback
      try {
        localStorage.setItem(
          `${LOCALSTORAGE_PREFIX}${key}`,
          JSON.stringify(entry)
        );
      } catch (error) {
        console.error('Failed to set item in localStorage', error);
        throw new Error('Failed to store data offline');
      }
    } else {
      // Use IndexedDB
      try {
        const db = await this.initIndexedDB();
        const transaction = db.transaction([STORAGE_STORE], 'readwrite');
        const store = transaction.objectStore(STORAGE_STORE);
        
        await new Promise<void>((resolve, reject) => {
          const request = store.put(entry);
          request.onsuccess = () => resolve();
          request.onerror = () => reject(new Error('Failed to store data in IndexedDB'));
        });
      } catch (error) {
        console.error('Failed to set item in IndexedDB', error);
        throw new Error('Failed to store data offline');
      }
    }
    
    // Notify storage change listeners
    this.notifyStorageChange(key, value);
    
    // Check storage limits
    this.checkStorageLimits();
  }
  
  /**
   * Retrieve an item from offline storage
   */
  async get<T = any>(key: StorageKey): Promise<T | null> {
    if (!this.ready) {
      await this.init();
    }
    
    if (this.useLocalStorageFallback) {
      // Use localStorage fallback
      const item = localStorage.getItem(`${LOCALSTORAGE_PREFIX}${key}`);
      
      if (!item) {
        return null;
      }
      
      try {
        const entry: StorageEntry = JSON.parse(item);
        
        // Check if item is expired
        if (entry.expiry !== null && Date.now() > entry.timestamp + entry.expiry) {
          localStorage.removeItem(`${LOCALSTORAGE_PREFIX}${key}`);
          return null;
        }
        
        // Decrypt if necessary
        if (entry.encrypted) {
          return this.decrypt(entry.value);
        }
        
        return entry.value;
      } catch (error) {
        console.error('Failed to parse item from localStorage', error);
        return null;
      }
    } else {
      // Use IndexedDB
      try {
        const db = await this.initIndexedDB();
        const transaction = db.transaction([STORAGE_STORE], 'readonly');
        const store = transaction.objectStore(STORAGE_STORE);
        
        const entry = await new Promise<StorageEntry | undefined>((resolve, reject) => {
          const request = store.get(key);
          request.onsuccess = () => resolve(request.result);
          request.onerror = () => reject(new Error('Failed to get data from IndexedDB'));
        });
        
        if (!entry) {
          return null;
        }
        
        // Check if item is expired
        if (entry.expiry !== null && Date.now() > entry.timestamp + entry.expiry) {
          await this.remove(key);
          return null;
        }
        
        // Decrypt if necessary
        if (entry.encrypted) {
          return this.decrypt(entry.value);
        }
        
        return entry.value;
      } catch (error) {
        console.error('Failed to get item from IndexedDB', error);
        return null;
      }
    }
  }
  
  /**
   * Remove an item from offline storage
   */
  async remove(key: StorageKey): Promise<void> {
    if (!this.ready) {
      await this.init();
    }
    
    if (this.useLocalStorageFallback) {
      // Use localStorage fallback
      localStorage.removeItem(`${LOCALSTORAGE_PREFIX}${key}`);
    } else {
      // Use IndexedDB
      try {
        const db = await this.initIndexedDB();
        const transaction = db.transaction([STORAGE_STORE], 'readwrite');
        const store = transaction.objectStore(STORAGE_STORE);
        
        await new Promise<void>((resolve, reject) => {
          const request = store.delete(key);
          request.onsuccess = () => resolve();
          request.onerror = () => reject(new Error('Failed to remove data from IndexedDB'));
        });
      } catch (error) {
        console.error('Failed to remove item from IndexedDB', error);
        throw new Error('Failed to remove data offline');
      }
    }
    
    // Notify storage change listeners
    this.notifyStorageChange(key, null);
  }
  
  /**
   * Check if an item exists in offline storage
   */
  async has(key: StorageKey): Promise<boolean> {
    if (!this.ready) {
      await this.init();
    }
    
    if (this.useLocalStorageFallback) {
      // Use localStorage fallback
      return localStorage.getItem(`${LOCALSTORAGE_PREFIX}${key}`) !== null;
    } else {
      // Use IndexedDB
      try {
        const db = await this.initIndexedDB();
        const transaction = db.transaction([STORAGE_STORE], 'readonly');
        const store = transaction.objectStore(STORAGE_STORE);
        
        const entry = await new Promise<StorageEntry | undefined>((resolve, reject) => {
          const request = store.get(key);
          request.onsuccess = () => resolve(request.result);
          request.onerror = () => reject(new Error('Failed to check data in IndexedDB'));
        });
        
        return !!entry;
      } catch (error) {
        console.error('Failed to check item in IndexedDB', error);
        return false;
      }
    }
  }
  
  /**
   * Queue a request to be processed when online
   */
  async queueRequest(
    url: string,
    method: string,
    body?: any,
    headers?: Record<string, string>
  ): Promise<string> {
    if (!this.ready) {
      await this.init();
    }
    
    const requestId = uuidv4();
    
    const queuedRequest: QueuedRequest = {
      id: requestId,
      url,
      method,
      body,
      headers,
      timestamp: Date.now(),
      retryCount: 0,
      status: 'pending',
    };
    
    if (this.useLocalStorageFallback) {
      // Use localStorage fallback
      const queueKey = `${LOCALSTORAGE_PREFIX}queue_${requestId}`;
      localStorage.setItem(queueKey, JSON.stringify(queuedRequest));
      
      // Update queue indices
      const queueIds = this.getQueueIdsFromLocalStorage();
      queueIds.push(requestId);
      localStorage.setItem(`${LOCALSTORAGE_PREFIX}queue_ids`, JSON.stringify(queueIds));
    } else {
      // Use IndexedDB
      try {
        const db = await this.initIndexedDB();
        const transaction = db.transaction([REQUEST_QUEUE_STORE], 'readwrite');
        const store = transaction.objectStore(REQUEST_QUEUE_STORE);
        
        await new Promise<void>((resolve, reject) => {
          const request = store.add(queuedRequest);
          request.onsuccess = () => resolve();
          request.onerror = () => reject(new Error('Failed to queue request in IndexedDB'));
        });
      } catch (error) {
        console.error('Failed to queue request in IndexedDB', error);
        throw new Error('Failed to queue request for offline processing');
      }
    }
    
    // Try to process the queue immediately if we might be online
    this.processRequestQueue();
    
    return requestId;
  }
  
  /**
   * Get all queued requests
   */
  async getQueuedRequests(): Promise<QueuedRequest[]> {
    if (!this.ready) {
      await this.init();
    }
    
    if (this.useLocalStorageFallback) {
      // Use localStorage fallback
      const queueIds = this.getQueueIdsFromLocalStorage();
      const requests: QueuedRequest[] = [];
      
      for (const id of queueIds) {
        const queueKey = `${LOCALSTORAGE_PREFIX}queue_${id}`;
        const requestJson = localStorage.getItem(queueKey);
        
        if (requestJson) {
          try {
            requests.push(JSON.parse(requestJson));
          } catch (error) {
            console.error('Failed to parse queued request from localStorage', error);
          }
        }
      }
      
      return requests;
    } else {
      // Use IndexedDB
      try {
        const db = await this.initIndexedDB();
        const transaction = db.transaction([REQUEST_QUEUE_STORE], 'readonly');
        const store = transaction.objectStore(REQUEST_QUEUE_STORE);
        
        return new Promise<QueuedRequest[]>((resolve, reject) => {
          const request = store.getAll();
          request.onsuccess = () => resolve(request.result);
          request.onerror = () => reject(new Error('Failed to get queued requests from IndexedDB'));
        });
      } catch (error) {
        console.error('Failed to get queued requests from IndexedDB', error);
        return [];
      }
    }
  }
  
  /**
   * Get pending request count
   */
  async getPendingRequestCount(): Promise<number> {
    const requests = await this.getQueuedRequests();
    return requests.filter(r => r.status === 'pending').length;
  }
  
  /**
   * Start processing the request queue periodically
   */
  private startQueueProcessing(): void {
    // Process queue initially
    this.processRequestQueue();
    
    // Set up interval to process queue
    if (this.queueProcessingInterval === null) {
      this.queueProcessingInterval = window.setInterval(
        () => this.processRequestQueue(),
        30000 // Check every 30 seconds
      );
    }
    
    // Add online event listener to process queue when we go back online
    window.addEventListener('online', () => {
      console.log('Back online, processing request queue');
      this.processRequestQueue();
    });
  }
  
  /**
   * Stop processing the request queue
   */
  stopQueueProcessing(): void {
    if (this.queueProcessingInterval !== null) {
      clearInterval(this.queueProcessingInterval);
      this.queueProcessingInterval = null;
    }
  }
  
  /**
   * Process the request queue
   */
  async processRequestQueue(): Promise<void> {
    // Skip processing if offline
    if (!navigator.onLine) {
      return;
    }
    
    const requests = await this.getQueuedRequests();
    const pendingRequests = requests.filter(r => r.status === 'pending');
    
    if (pendingRequests.length === 0) {
      return;
    }
    
    console.log(`Processing ${pendingRequests.length} pending requests`);
    
    // Process each request
    const successfulRequests: string[] = [];
    
    for (const request of pendingRequests) {
      try {
        // Mark as processing
        await this.updateRequestStatus(request.id, 'processing');
        
        // Execute the request
        const response = await fetch(request.url, {
          method: request.method,
          headers: {
            'Content-Type': 'application/json',
            ...(request.headers || {}),
          },
          body: request.body ? JSON.stringify(request.body) : undefined,
        });
        
        if (response.ok) {
          // Request succeeded, remove from queue
          await this.removeRequest(request.id);
          successfulRequests.push(request.id);
        } else {
          // Request failed, update status and retry count
          const errorText = await response.text();
          
          await this.updateRequest(request.id, {
            status: 'failed',
            retryCount: request.retryCount + 1,
            errorMessage: `HTTP ${response.status}: ${errorText}`,
          });
        }
      } catch (error) {
        // Network error, update status and retry count
        await this.updateRequest(request.id, {
          status: 'failed',
          retryCount: request.retryCount + 1,
          errorMessage: error instanceof Error ? error.message : String(error),
        });
      }
    }
    
    // Notify sync success handlers if any requests succeeded
    if (successfulRequests.length > 0) {
      this.notifySyncSuccess(successfulRequests);
    }
  }
  
  /**
   * Update a request's status
   */
  private async updateRequestStatus(id: string, status: 'pending' | 'processing' | 'failed'): Promise<void> {
    return this.updateRequest(id, { status });
  }
  
  /**
   * Update a request
   */
  private async updateRequest(id: string, updates: Partial<QueuedRequest>): Promise<void> {
    if (this.useLocalStorageFallback) {
      // Use localStorage fallback
      const queueKey = `${LOCALSTORAGE_PREFIX}queue_${id}`;
      const requestJson = localStorage.getItem(queueKey);
      
      if (requestJson) {
        try {
          const request = JSON.parse(requestJson);
          const updatedRequest = { ...request, ...updates };
          localStorage.setItem(queueKey, JSON.stringify(updatedRequest));
        } catch (error) {
          console.error('Failed to update queued request in localStorage', error);
        }
      }
    } else {
      // Use IndexedDB
      try {
        const db = await this.initIndexedDB();
        const transaction = db.transaction([REQUEST_QUEUE_STORE], 'readwrite');
        const store = transaction.objectStore(REQUEST_QUEUE_STORE);
        
        // Get the current request
        const request = await new Promise<QueuedRequest | undefined>((resolve, reject) => {
          const getRequest = store.get(id);
          getRequest.onsuccess = () => resolve(getRequest.result);
          getRequest.onerror = () => reject(new Error('Failed to get request from IndexedDB'));
        });
        
        if (request) {
          // Update the request
          const updatedRequest = { ...request, ...updates };
          
          await new Promise<void>((resolve, reject) => {
            const putRequest = store.put(updatedRequest);
            putRequest.onsuccess = () => resolve();
            putRequest.onerror = () => reject(new Error('Failed to update request in IndexedDB'));
          });
        }
      } catch (error) {
        console.error('Failed to update request in IndexedDB', error);
      }
    }
  }
  
  /**
   * Remove a request from the queue
   */
  private async removeRequest(id: string): Promise<void> {
    if (this.useLocalStorageFallback) {
      // Use localStorage fallback
      const queueKey = `${LOCALSTORAGE_PREFIX}queue_${id}`;
      localStorage.removeItem(queueKey);
      
      // Update queue indices
      const queueIds = this.getQueueIdsFromLocalStorage();
      const updatedQueueIds = queueIds.filter(queueId => queueId !== id);
      localStorage.setItem(`${LOCALSTORAGE_PREFIX}queue_ids`, JSON.stringify(updatedQueueIds));
    } else {
      // Use IndexedDB
      try {
        const db = await this.initIndexedDB();
        const transaction = db.transaction([REQUEST_QUEUE_STORE], 'readwrite');
        const store = transaction.objectStore(REQUEST_QUEUE_STORE);
        
        await new Promise<void>((resolve, reject) => {
          const request = store.delete(id);
          request.onsuccess = () => resolve();
          request.onerror = () => reject(new Error('Failed to remove request from IndexedDB'));
        });
      } catch (error) {
        console.error('Failed to remove request from IndexedDB', error);
      }
    }
  }
  
  /**
   * Get queue IDs from localStorage
   */
  private getQueueIdsFromLocalStorage(): string[] {
    const queueIdsJson = localStorage.getItem(`${LOCALSTORAGE_PREFIX}queue_ids`);
    
    if (queueIdsJson) {
      try {
        return JSON.parse(queueIdsJson);
      } catch (error) {
        console.error('Failed to parse queue IDs from localStorage', error);
      }
    }
    
    return [];
  }
  
  /**
   * Clean up expired items
   */
  private async cleanupExpiredItems(): Promise<void> {
    if (!this.ready) {
      await this.init();
    }
    
    const now = Date.now();
    
    if (this.useLocalStorageFallback) {
      // Use localStorage fallback
      const keysToRemove: string[] = [];
      
      for (let i = 0; i < localStorage.length; i++) {
        const key = localStorage.key(i);
        
        if (key && key.startsWith(LOCALSTORAGE_PREFIX) && !key.startsWith(`${LOCALSTORAGE_PREFIX}queue_`)) {
          const item = localStorage.getItem(key);
          
          if (item) {
            try {
              const entry: StorageEntry = JSON.parse(item);
              
              if (entry.expiry !== null && now > entry.timestamp + entry.expiry) {
                keysToRemove.push(key);
              }
            } catch (error) {
              // If we can't parse it, remove it
              keysToRemove.push(key);
            }
          }
        }
      }
      
      // Remove expired items
      for (const key of keysToRemove) {
        localStorage.removeItem(key);
      }
      
      if (keysToRemove.length > 0) {
        console.log(`Cleaned up ${keysToRemove.length} expired items from localStorage`);
      }
    } else {
      // Use IndexedDB
      try {
        const db = await this.initIndexedDB();
        const transaction = db.transaction([STORAGE_STORE], 'readwrite');
        const store = transaction.objectStore(STORAGE_STORE);
        const index = store.index('expiry');
        
        // Get all items with expiry (not null)
        const expiryRange = IDBKeyRange.lowerBound(0);
        const request = index.openCursor(expiryRange);
        
        let removeCount = 0;
        
        await new Promise<void>((resolve) => {
          request.onsuccess = (event) => {
            const cursor = (event.target as IDBRequest).result as IDBCursorWithValue;
            
            if (cursor) {
              const entry = cursor.value as StorageEntry;
              
              // Check if the item is expired
              if (entry.expiry !== null && now > entry.timestamp + entry.expiry) {
                cursor.delete();
                removeCount++;
              }
              
              cursor.continue();
            } else {
              if (removeCount > 0) {
                console.log(`Cleaned up ${removeCount} expired items from IndexedDB`);
              }
              resolve();
            }
          };
          
          request.onerror = () => {
            console.error('Failed to cleanup expired items in IndexedDB');
            resolve();
          };
        });
      } catch (error) {
        console.error('Failed to cleanup expired items', error);
      }
    }
    
    // Schedule next cleanup
    setTimeout(() => this.cleanupExpiredItems(), 60 * 60 * 1000); // Run every hour
  }
  
  /**
   * Check and enforce storage limits
   */
  private async checkStorageLimits(): Promise<void> {
    // Only check periodically to avoid performance issues
    const metadata = await this.getMetadata('lastStorageCheck');
    const now = Date.now();
    
    if (metadata && now - metadata.timestamp < 60 * 60 * 1000) {
      // Last check was less than an hour ago
      return;
    }
    
    // Update last check timestamp
    await this.setMetadata('lastStorageCheck', { timestamp: now });
    
    if (this.useLocalStorageFallback) {
      // Check localStorage usage
      try {
        let totalSize = 0;
        const prefix = LOCALSTORAGE_PREFIX;
        
        for (let i = 0; i < localStorage.length; i++) {
          const key = localStorage.key(i);
          
          if (key && key.startsWith(prefix)) {
            const value = localStorage.getItem(key) || '';
            totalSize += key.length + value.length;
          }
        }
        
        // Convert to MB
        const usedMB = totalSize / (1024 * 1024);
        
        // If we're approaching the limit, clean up old items
        if (usedMB > MAX_STORAGE_SIZE * 0.8) {
          console.log(`Storage usage high (${usedMB.toFixed(2)}MB), cleaning up oldest items`);
          await this.cleanupOldestItems();
        }
      } catch (error) {
        console.error('Failed to check localStorage limits', error);
      }
    } else {
      // Check IndexedDB usage
      try {
        // Use the experimental storage manager API if available
        if (navigator.storage && navigator.storage.estimate) {
          const estimate = await navigator.storage.estimate();
          
          if (estimate.usage && estimate.quota) {
            const usedMB = estimate.usage / (1024 * 1024);
            const quotaMB = estimate.quota / (1024 * 1024);
            const usagePercentage = (estimate.usage / estimate.quota) * 100;
            
            console.log(`Storage usage: ${usedMB.toFixed(2)}MB of ${quotaMB.toFixed(2)}MB (${usagePercentage.toFixed(2)}%)`);
            
            // If we're approaching the limit, clean up old items
            if (usagePercentage > 80) {
              console.log('Storage usage high, cleaning up oldest items');
              await this.cleanupOldestItems();
            }
          }
        }
      } catch (error) {
        console.error('Failed to check IndexedDB limits', error);
      }
    }
  }
  
  /**
   * Clean up oldest items to free up space
   */
  private async cleanupOldestItems(): Promise<void> {
    if (this.useLocalStorageFallback) {
      // Use localStorage fallback
      const items: { key: string; timestamp: number }[] = [];
      const prefix = LOCALSTORAGE_PREFIX;
      
      // Collect all items with their timestamps
      for (let i = 0; i < localStorage.length; i++) {
        const key = localStorage.key(i);
        
        if (key && key.startsWith(prefix) && !key.startsWith(`${prefix}queue_`)) {
          const item = localStorage.getItem(key);
          
          if (item) {
            try {
              const entry: StorageEntry = JSON.parse(item);
              items.push({ key, timestamp: entry.timestamp });
            } catch (error) {
              // If we can't parse it, assume it's old
              items.push({ key, timestamp: 0 });
            }
          }
        }
      }
      
      // Sort by timestamp (oldest first)
      items.sort((a, b) => a.timestamp - b.timestamp);
      
      // Remove oldest 20% of items
      const itemsToRemove = Math.ceil(items.length * 0.2);
      
      for (let i = 0; i < itemsToRemove && i < items.length; i++) {
        localStorage.removeItem(items[i].key);
      }
      
      console.log(`Removed ${itemsToRemove} oldest items from localStorage`);
    } else {
      // Use IndexedDB
      try {
        const db = await this.initIndexedDB();
        const transaction = db.transaction([STORAGE_STORE], 'readwrite');
        const store = transaction.objectStore(STORAGE_STORE);
        const index = store.index('timestamp');
        
        // Get all items sorted by timestamp
        const request = index.openCursor();
        
        // Count total items to determine how many to remove
        const countRequest = store.count();
        const totalItems = await new Promise<number>((resolve) => {
          countRequest.onsuccess = () => resolve(countRequest.result);
          countRequest.onerror = () => resolve(0);
        });
        
        const itemsToRemove = Math.ceil(totalItems * 0.2); // Remove oldest 20%
        let removedCount = 0;
        
        await new Promise<void>((resolve) => {
          request.onsuccess = (event) => {
            const cursor = (event.target as IDBRequest).result as IDBCursorWithValue;
            
            if (cursor && removedCount < itemsToRemove) {
              // Skip queue-related items
              const entry = cursor.value as StorageEntry;
              const key = entry.key as string;
              
              if (!key.startsWith('queue_')) {
                cursor.delete();
                removedCount++;
              }
              
              cursor.continue();
            } else {
              console.log(`Removed ${removedCount} oldest items from IndexedDB`);
              resolve();
            }
          };
          
          request.onerror = () => {
            console.error('Failed to cleanup oldest items in IndexedDB');
            resolve();
          };
        });
      } catch (error) {
        console.error('Failed to cleanup oldest items', error);
      }
    }
  }
  
  /**
   * Get metadata from storage
   */
  private async getMetadata(key: string): Promise<any> {
    if (!this.ready) {
      await this.init();
    }
    
    if (this.useLocalStorageFallback) {
      // Use localStorage fallback
      const metadataKey = `${LOCALSTORAGE_PREFIX}metadata_${key}`;
      const metadataJson = localStorage.getItem(metadataKey);
      
      if (metadataJson) {
        try {
          return JSON.parse(metadataJson);
        } catch (error) {
          console.error('Failed to parse metadata from localStorage', error);
        }
      }
      
      return null;
    } else {
      // Use IndexedDB
      try {
        const db = await this.initIndexedDB();
        const transaction = db.transaction([METADATA_STORE], 'readonly');
        const store = transaction.objectStore(METADATA_STORE);
        
        return new Promise<any>((resolve, reject) => {
          const request = store.get(key);
          request.onsuccess = () => resolve(request.result?.value);
          request.onerror = () => reject(new Error('Failed to get metadata from IndexedDB'));
        });
      } catch (error) {
        console.error('Failed to get metadata from IndexedDB', error);
        return null;
      }
    }
  }
  
  /**
   * Set metadata in storage
   */
  private async setMetadata(key: string, value: any): Promise<void> {
    if (!this.ready) {
      await this.init();
    }
    
    if (this.useLocalStorageFallback) {
      // Use localStorage fallback
      const metadataKey = `${LOCALSTORAGE_PREFIX}metadata_${key}`;
      
      try {
        localStorage.setItem(metadataKey, JSON.stringify(value));
      } catch (error) {
        console.error('Failed to set metadata in localStorage', error);
      }
    } else {
      // Use IndexedDB
      try {
        const db = await this.initIndexedDB();
        const transaction = db.transaction([METADATA_STORE], 'readwrite');
        const store = transaction.objectStore(METADATA_STORE);
        
        await new Promise<void>((resolve, reject) => {
          const request = store.put({ key, value });
          request.onsuccess = () => resolve();
          request.onerror = () => reject(new Error('Failed to set metadata in IndexedDB'));
        });
      } catch (error) {
        console.error('Failed to set metadata in IndexedDB', error);
      }
    }
  }
  
  /**
   * Register a callback for successful synchronization
   */
  onSyncSuccess(callback: (keys: StorageKey[]) => void): () => void {
    this.onSyncSuccessHandlers.push(callback);
    
    // Return a function to unregister the callback
    return () => {
      const index = this.onSyncSuccessHandlers.indexOf(callback);
      if (index !== -1) {
        this.onSyncSuccessHandlers.splice(index, 1);
      }
    };
  }
  
  /**
   * Register a callback for sync failures
   */
  onSyncFailure(callback: (error: Error) => void): () => void {
    this.onSyncFailureHandlers.push(callback);
    
    // Return a function to unregister the callback
    return () => {
      const index = this.onSyncFailureHandlers.indexOf(callback);
      if (index !== -1) {
        this.onSyncFailureHandlers.splice(index, 1);
      }
    };
  }
  
  /**
   * Register a callback for storage changes
   */
  onStorageChange(callback: (key: StorageKey, newValue: StorageValue | null) => void): () => void {
    this.onStorageChangeHandlers.push(callback);
    
    // Return a function to unregister the callback
    return () => {
      const index = this.onStorageChangeHandlers.indexOf(callback);
      if (index !== -1) {
        this.onStorageChangeHandlers.splice(index, 1);
      }
    };
  }
  
  /**
   * Notify sync success handlers
   */
  private notifySyncSuccess(keys: StorageKey[]): void {
    for (const handler of this.onSyncSuccessHandlers) {
      try {
        handler(keys);
      } catch (error) {
        console.error('Error in sync success handler', error);
      }
    }
  }
  
  /**
   * Notify sync failure handlers
   */
  private notifySyncFailure(error: Error): void {
    for (const handler of this.onSyncFailureHandlers) {
      try {
        handler(error);
      } catch (handlerError) {
        console.error('Error in sync failure handler', handlerError);
      }
    }
  }
  
  /**
   * Notify storage change handlers
   */
  private notifyStorageChange(key: StorageKey, newValue: StorageValue | null): void {
    for (const handler of this.onStorageChangeHandlers) {
      try {
        handler(key, newValue);
      } catch (error) {
        console.error('Error in storage change handler', error);
      }
    }
  }
  
  /**
   * Get all keys with a specific prefix
   */
  async getKeysByPrefix(prefix: string): Promise<StorageKey[]> {
    if (!this.ready) {
      await this.init();
    }
    
    if (this.useLocalStorageFallback) {
      // Use localStorage fallback
      const keys: StorageKey[] = [];
      const storagePrefix = LOCALSTORAGE_PREFIX;
      
      for (let i = 0; i < localStorage.length; i++) {
        const key = localStorage.key(i);
        
        if (key && key.startsWith(storagePrefix)) {
          const actualKey = key.substring(storagePrefix.length);
          
          if (actualKey.startsWith(prefix)) {
            keys.push(actualKey);
          }
        }
      }
      
      return keys;
    } else {
      // Use IndexedDB
      try {
        const db = await this.initIndexedDB();
        const transaction = db.transaction([STORAGE_STORE], 'readonly');
        const store = transaction.objectStore(STORAGE_STORE);
        
        const allKeys = await new Promise<StorageKey[]>((resolve, reject) => {
          const request = store.getAllKeys();
          request.onsuccess = () => resolve(request.result as StorageKey[]);
          request.onerror = () => reject(new Error('Failed to get keys from IndexedDB'));
        });
        
        // Filter keys by prefix
        return allKeys.filter(key => 
          typeof key === 'string' && key.startsWith(prefix)
        );
      } catch (error) {
        console.error('Failed to get keys by prefix from IndexedDB', error);
        return [];
      }
    }
  }
  
  /**
   * Clear all stored data
   */
  async clear(): Promise<void> {
    if (!this.ready) {
      await this.init();
    }
    
    if (this.useLocalStorageFallback) {
      // Use localStorage fallback
      const keysToRemove: string[] = [];
      const prefix = LOCALSTORAGE_PREFIX;
      
      for (let i = 0; i < localStorage.length; i++) {
        const key = localStorage.key(i);
        
        if (key && key.startsWith(prefix)) {
          keysToRemove.push(key);
        }
      }
      
      // Remove all items
      for (const key of keysToRemove) {
        localStorage.removeItem(key);
      }
      
      console.log(`Cleared ${keysToRemove.length} items from localStorage`);
    } else {
      // Use IndexedDB
      try {
        const db = await this.initIndexedDB();
        
        // Clear all object stores
        const transaction = db.transaction([STORAGE_STORE, REQUEST_QUEUE_STORE, METADATA_STORE], 'readwrite');
        
        await Promise.all([
          new Promise<void>((resolve, reject) => {
            const request = transaction.objectStore(STORAGE_STORE).clear();
            request.onsuccess = () => resolve();
            request.onerror = () => reject(new Error('Failed to clear storage store'));
          }),
          new Promise<void>((resolve, reject) => {
            const request = transaction.objectStore(REQUEST_QUEUE_STORE).clear();
            request.onsuccess = () => resolve();
            request.onerror = () => reject(new Error('Failed to clear request queue store'));
          }),
          new Promise<void>((resolve, reject) => {
            const request = transaction.objectStore(METADATA_STORE).clear();
            request.onsuccess = () => resolve();
            request.onerror = () => reject(new Error('Failed to clear metadata store'));
          }),
        ]);
        
        console.log('Cleared all data from IndexedDB');
      } catch (error) {
        console.error('Failed to clear IndexedDB data', error);
        throw new Error('Failed to clear offline storage');
      }
    }
  }
  
  /**
   * Close the offline storage
   */
  close(): void {
    // Stop queue processing
    this.stopQueueProcessing();
    
    // Close IndexedDB connection if open
    if (this.db) {
      this.db.close();
      this.db = null;
      this.dbPromise = null;
    }
    
    this.ready = false;
    console.log('Offline storage closed');
  }
}

// Create a singleton instance
let offlineStorageInstance: OfflineStorage | null = null;

/**
 * Get the offline storage instance
 */
export function getOfflineStorage(): OfflineStorage {
  if (!offlineStorageInstance) {
    offlineStorageInstance = new OfflineStorage();
    
    // Initialize the storage
    offlineStorageInstance.init().catch(error => {
      console.error('Failed to initialize offline storage', error);
    });
  }
  
  return offlineStorageInstance;
}
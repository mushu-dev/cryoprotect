# Offline Components

This module provides a comprehensive solution for adding offline capabilities to your application. It includes components for handling offline states, synchronization, and data persistence.

## Components

### OfflineWrapper

The `OfflineWrapper` component provides offline capabilities to any page or component:

```tsx
<OfflineWrapper
  fallback={<CustomOfflineFallback />}
  requireOnline={false}
  showStatusIndicator={true}
  onConnectionChange={(status) => console.log(`Connection status: ${status}`)}
>
  <YourComponent />
</OfflineWrapper>
```

#### Props

- `children`: The content to render
- `fallback` (optional): Custom fallback UI when offline
- `apiClient` (optional): Custom API client instance
- `requireOnline` (optional): Whether the wrapped component requires an internet connection
- `showStatusIndicator` (optional): Whether to show the offline status indicator
- `onConnectionChange` (optional): Callback when connection status changes

### OfflineStatusIndicator

The `OfflineStatusIndicator` displays the current connection status:

```tsx
<OfflineStatusIndicator className="custom-class" />
```

#### Props

- `className` (optional): Additional CSS classes

## Hook

### useOfflineMode

The `useOfflineMode` hook provides access to offline functionality from any component:

```tsx
const {
  isOnline,
  connectionStatus,
  hasPendingChanges,
  pendingChangesCount,
  isSyncing,
  lastSyncTime,
  formattedLastSyncTime,
  synchronize,
  apiClient,
  syncService
} = useOfflineMode();
```

#### Options

- `apiClient` (optional): Custom API client instance

#### Return Values

- `isOnline`: Boolean indicating if the application is online
- `connectionStatus`: Current connection status ('connected', 'degraded', or 'offline')
- `hasPendingChanges`: Boolean indicating if there are pending changes
- `pendingChangesCount`: Number of pending changes
- `isSyncing`: Boolean indicating if synchronization is in progress
- `lastSyncTime`: Timestamp of the last successful sync
- `formattedLastSyncTime`: Human-readable version of last sync time
- `synchronize`: Function to trigger synchronization
- `apiClient`: Reference to the API client
- `syncService`: Reference to the sync service

## Usage Examples

### Basic Usage

Add offline capability to a page:

```tsx
// pages/example.tsx
import { OfflineWrapper } from '@/components/offline';

const ExamplePage = () => {
  return (
    <OfflineWrapper>
      <YourComponent />
    </OfflineWrapper>
  );
};

export default ExamplePage;
```

### Using the Hook

Access offline functionality in a component:

```tsx
// components/DataEditor.tsx
import { useOfflineMode } from '@/components/offline';

const DataEditor = () => {
  const { isOnline, hasPendingChanges, synchronize } = useOfflineMode();
  
  const saveData = async (data) => {
    // Save data locally
    await localStorage.setItem('saved_data', JSON.stringify(data));
    
    // Try to sync if online
    if (isOnline) {
      await synchronize();
    }
  };
  
  return (
    <div>
      {hasPendingChanges && (
        <button onClick={() => synchronize()}>
          Sync Pending Changes
        </button>
      )}
      
      {/* Rest of your component */}
    </div>
  );
};
```

### Requiring Online Connection

Create a page that requires an internet connection:

```tsx
// pages/admin.tsx
import { OfflineWrapper } from '@/components/offline';

const AdminPage = () => {
  return (
    <OfflineWrapper 
      requireOnline={true}
      fallback={<CustomOfflineFallback />}
    >
      <AdminDashboard />
    </OfflineWrapper>
  );
};

export default AdminPage;
```

## Advanced Usage

### Custom API Client

Use a custom API client with additional configuration:

```tsx
// Create a custom API client
const apiClient = new ResilientApiClient({
  baseURL: '/api/v2',
  timeout: 30000,
  retries: 5,
  authToken: 'your-auth-token',
  enableCache: true
});

// Use it with the OfflineWrapper
<OfflineWrapper apiClient={apiClient}>
  <YourComponent />
</OfflineWrapper>

// Or with the hook
const { isOnline } = useOfflineMode({ apiClient });
```

### Listening for Connection Changes

React to connection status changes:

```tsx
<OfflineWrapper
  onConnectionChange={(status) => {
    if (status === 'connected') {
      showNotification('You are back online!');
    } else if (status === 'offline') {
      showNotification('You are offline. Changes will be saved locally.');
    }
  }}
>
  <YourComponent />
</OfflineWrapper>
```

### Working with Offline Data

Store and retrieve data while offline:

```tsx
import { getOfflineStorage } from '@/services/offline-storage';
import { useOfflineMode } from '@/components/offline';

const DataComponent = () => {
  const [data, setData] = useState(null);
  const { isOnline } = useOfflineMode();
  const offlineStorage = getOfflineStorage();
  
  // Load data
  useEffect(() => {
    const loadData = async () => {
      // Try to get from API if online
      if (isOnline) {
        try {
          const apiData = await fetchFromApi();
          setData(apiData);
          
          // Save to offline storage for later
          await offlineStorage.set('my_data', apiData);
          return;
        } catch (error) {
          console.error('Failed to load from API:', error);
        }
      }
      
      // Fall back to offline storage
      const offlineData = await offlineStorage.get('my_data');
      if (offlineData) {
        setData(offlineData);
      }
    };
    
    loadData();
  }, [isOnline]);
  
  // Save data
  const saveData = async (newData) => {
    // Update state
    setData(newData);
    
    // Save to offline storage
    await offlineStorage.set('my_data', newData);
    
    // Try to save to API if online
    if (isOnline) {
      try {
        await saveToApi(newData);
      } catch (error) {
        console.error('Failed to save to API:', error);
      }
    }
  };
  
  return (
    <div>
      {/* Your component UI */}
    </div>
  );
};
```

## Architecture

The offline components are built on top of several services:

1. **ResilientApiClient**: Provides robust API communication with retry logic, circuit breaking, and error handling.

2. **SyncService**: Manages synchronization of offline changes with the server, including conflict resolution.

3. **OfflineStorage**: Handles persistent storage using IndexedDB with localStorage fallback.

These services work together to provide a seamless offline experience for users.

## Best Practices

1. **Always wrap with OfflineWrapper**: Wrap your pages or top-level components with the OfflineWrapper component.

2. **Store critical data offline**: Use the `offlineStorage.set()` method to persist important data.

3. **React to connection changes**: Use `onConnectionChange` or the status from `useOfflineMode` to adapt your UI.

4. **Provide fallbacks**: Create appropriate fallback UIs for when users are offline.

5. **Consider data size**: Be mindful of storage limits when caching large datasets offline.

6. **Handle conflicts**: Implement proper conflict resolution when synchronizing data.

7. **Show sync status**: Let users know when changes are pending synchronization.
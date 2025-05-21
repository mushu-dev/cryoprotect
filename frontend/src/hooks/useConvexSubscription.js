/**
 * Convex Subscription Hook
 * 
 * This hook provides real-time subscriptions to Convex data for collaborative features
 * using the Convex client's subscription mechanism.
 */
import { useState, useEffect } from 'react';
import { useQuery } from 'convex/react';

/**
 * Custom hook for real-time data subscriptions from Convex
 * 
 * @param {Function} queryFunction - The Convex query function to subscribe to
 * @param {Array} args - Arguments to pass to the query function
 * @param {Object} options - Additional options
 * @param {Function} options.onDataChange - Callback when data changes
 * @returns {Object} Subscription state with data, loading, and error
 */
export function useConvexSubscription(queryFunction, args = [], options = {}) {
  const [state, setState] = useState({
    data: null,
    previousData: null,
    isLoading: true,
    error: null,
    lastUpdated: null
  });

  // Spread args to ensure they're passed correctly to useQuery
  const result = useQuery(queryFunction, ...args);

  useEffect(() => {
    // Check if the data has changed
    const hasChanged = JSON.stringify(result) !== JSON.stringify(state.data);
    
    if (hasChanged) {
      const previousData = state.data;
      const newState = {
        data: result,
        previousData,
        isLoading: false,
        error: null,
        lastUpdated: new Date()
      };
      
      setState(newState);
      
      // Call the onDataChange callback if provided
      if (options.onDataChange && previousData !== null) {
        options.onDataChange({
          previousData,
          currentData: result,
          changes: detectChanges(previousData, result)
        });
      }
    }
  }, [result, state.data, options]);

  return {
    ...state,
    refresh: () => setState({ ...state, isLoading: true })
  };
}

/**
 * Hook specifically for collaborative document editing
 * 
 * @param {String} documentId - The ID of the document being edited
 * @param {String} documentType - The type of document (e.g., 'protocols', 'experiments')
 * @param {Object} options - Additional options
 * @returns {Object} Document state with collaboration features
 */
export function useCollaborativeDocument(documentId, documentType, options = {}) {
  const [collaborators, setCollaborators] = useState([]);
  const [localChanges, setLocalChanges] = useState({});
  const [mergeConflicts, setMergeConflicts] = useState([]);
  const [lastSynced, setLastSynced] = useState(null);

  // Subscribe to the document data
  const { data, isLoading, error, lastUpdated } = useConvexSubscription(
    // This would be the actual Convex query function, e.g., api.documents.get
    `api.${documentType}.get`,
    [documentId],
    {
      onDataChange: ({ previousData, currentData, changes }) => {
        // Handle remote changes
        handleRemoteChanges(previousData, currentData, changes);
        
        // Call the onDataChange callback if provided
        if (options.onDataChange) {
          options.onDataChange({ previousData, currentData, changes });
        }
      }
    }
  );

  // Subscribe to active collaborators
  useConvexSubscription(
    // This would be the actual Convex query function, e.g., api.collaborators.list
    'api.collaborators.list',
    [documentId],
    {
      onDataChange: ({ currentData }) => {
        setCollaborators(currentData || []);
      }
    }
  );

  // Function to handle remote changes
  const handleRemoteChanges = (previousData, currentData, changes) => {
    // Check for conflicts with local changes
    const conflicts = detectConflicts(localChanges, changes);
    
    if (conflicts.length > 0) {
      setMergeConflicts(conflicts);
    } else {
      // No conflicts, clear any local changes that were successfully synced
      const remainingChanges = { ...localChanges };
      
      Object.keys(changes).forEach(key => {
        if (remainingChanges[key] !== undefined) {
          delete remainingChanges[key];
        }
      });
      
      setLocalChanges(remainingChanges);
    }
    
    setLastSynced(new Date());
  };

  // Function to make local changes
  const updateDocument = async (changes) => {
    // Update local state immediately for optimistic UI
    setLocalChanges({ ...localChanges, ...changes });
    
    try {
      // This would be the actual Convex mutation function, e.g., api.documents.update
      // await mutate(`api.${documentType}.update`, { id: documentId, ...changes });
      
      // After successful update, clear these changes from localChanges
      // This will happen automatically when the subscription receives the update
    } catch (error) {
      console.error('Error updating document:', error);
      // Keep the changes in localChanges for potential retry
    }
  };

  // Function to resolve merge conflicts
  const resolveConflicts = (resolutions) => {
    // Apply the conflict resolutions
    const resolvedChanges = {};
    
    resolutions.forEach(({ field, resolution }) => {
      resolvedChanges[field] = resolution;
    });
    
    // Update the document with the resolved conflicts
    updateDocument(resolvedChanges);
    
    // Clear the conflicts
    setMergeConflicts([]);
  };

  return {
    document: data,
    isLoading,
    error,
    lastUpdated,
    collaborators,
    localChanges,
    mergeConflicts,
    lastSynced,
    updateDocument,
    resolveConflicts
  };
}

/**
 * Detect changes between previous and current data
 * 
 * @param {Object} previousData - The previous data object
 * @param {Object} currentData - The current data object
 * @returns {Object} Object containing the changed fields and their new values
 */
function detectChanges(previousData, currentData) {
  if (!previousData || !currentData) return {};

  const changes = {};
  
  // Find fields that have changed
  Object.keys(currentData).forEach(key => {
    if (JSON.stringify(previousData[key]) !== JSON.stringify(currentData[key])) {
      changes[key] = currentData[key];
    }
  });
  
  return changes;
}

/**
 * Detect conflicts between local changes and remote changes
 * 
 * @param {Object} localChanges - Local changes that haven't been synced yet
 * @param {Object} remoteChanges - Changes that came from the server
 * @returns {Array} Array of conflict objects
 */
function detectConflicts(localChanges, remoteChanges) {
  const conflicts = [];
  
  Object.keys(localChanges).forEach(key => {
    if (remoteChanges[key] !== undefined && 
        JSON.stringify(localChanges[key]) !== JSON.stringify(remoteChanges[key])) {
      conflicts.push({
        field: key,
        localValue: localChanges[key],
        remoteValue: remoteChanges[key]
      });
    }
  });
  
  return conflicts;
}
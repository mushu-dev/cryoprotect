/**
 * Collaborative Experiment Editor Component
 * 
 * This component provides real-time collaborative editing for experiment data
 * using Convex's real-time subscriptions.
 */
import React, { useState, useEffect } from 'react';
import { useCollaborativeDocument } from '../hooks/useConvexSubscription';
import { useMutation } from 'convex/react';
import { api } from '../convex/_generated/api';

const CollaborativeExperimentEditor = ({ experimentId, onSave }) => {
  const [formData, setFormData] = useState({});
  const [isSaving, setIsSaving] = useState(false);
  const [isDirty, setIsDirty] = useState(false);
  
  // Use the collaborative document hook
  const {
    document: experiment,
    isLoading,
    error,
    collaborators,
    localChanges,
    mergeConflicts,
    lastSynced,
    updateDocument,
    resolveConflicts
  } = useCollaborativeDocument(experimentId, 'experiments', {
    onDataChange: ({ currentData }) => {
      // When data changes from the server, update our local form state
      // but preserve any unsaved changes
      setFormData(prevData => ({
        ...currentData,
        ...prevData // Keep any unsaved local changes
      }));
    }
  });
  
  // Convex mutation for saving experiment data
  const updateExperiment = useMutation(api.experiments.update);
  
  // Initialize form data when the experiment loads
  useEffect(() => {
    if (experiment && !isLoading && !isDirty) {
      setFormData(experiment);
    }
  }, [experiment, isLoading, isDirty]);
  
  // Handle form field changes
  const handleChange = (field, value) => {
    // Update local form state
    setFormData(prev => ({
      ...prev,
      [field]: value
    }));
    
    // Mark the form as dirty (has unsaved changes)
    setIsDirty(true);
  };
  
  // Save changes to Convex
  const handleSave = async () => {
    if (!isDirty) return;
    
    setIsSaving(true);
    
    try {
      // Prepare changes object with only modified fields
      const changes = {};
      Object.keys(formData).forEach(key => {
        if (JSON.stringify(formData[key]) !== JSON.stringify(experiment?.[key])) {
          changes[key] = formData[key];
        }
      });
      
      // Update document in Convex
      await updateExperiment({ id: experimentId, ...changes });
      
      // Call the real-time update function for optimistic UI
      updateDocument(changes);
      
      // Reset form state
      setIsDirty(false);
      
      // Call the onSave callback if provided
      if (onSave) onSave(formData);
    } catch (error) {
      console.error('Error saving experiment:', error);
    } finally {
      setIsSaving(false);
    }
  };
  
  // Handle conflict resolution
  const handleResolveConflict = (field, resolution) => {
    resolveConflicts([{ field, resolution }]);
    
    // Update form data with the resolution
    setFormData(prev => ({
      ...prev,
      [field]: resolution
    }));
  };
  
  // Render loading state
  if (isLoading) {
    return <div className="loading">Loading experiment data...</div>;
  }
  
  // Render error state
  if (error) {
    return <div className="error">Error loading experiment: {error.message}</div>;
  }
  
  // Render merge conflicts if any
  if (mergeConflicts.length > 0) {
    return (
      <div className="merge-conflicts">
        <h3>Merge Conflicts</h3>
        <p>Another user has made changes to this experiment. Please resolve the conflicts:</p>
        
        {mergeConflicts.map(conflict => (
          <div key={conflict.field} className="conflict">
            <h4>{conflict.field}</h4>
            <div className="conflict-options">
              <div className="local-change">
                <h5>Your changes:</h5>
                <pre>{JSON.stringify(conflict.localValue, null, 2)}</pre>
                <button onClick={() => handleResolveConflict(conflict.field, conflict.localValue)}>
                  Keep my changes
                </button>
              </div>
              
              <div className="remote-change">
                <h5>Their changes:</h5>
                <pre>{JSON.stringify(conflict.remoteValue, null, 2)}</pre>
                <button onClick={() => handleResolveConflict(conflict.field, conflict.remoteValue)}>
                  Use their changes
                </button>
              </div>
            </div>
          </div>
        ))}
      </div>
    );
  }
  
  // Render active collaborators
  const renderCollaborators = () => {
    if (!collaborators || collaborators.length === 0) return null;
    
    return (
      <div className="collaborators">
        <h4>Active Collaborators:</h4>
        <ul>
          {collaborators.map(user => (
            <li key={user.id} className="collaborator">
              <span className="collaborator-indicator" style={{ backgroundColor: user.color }}></span>
              <span className="collaborator-name">{user.name || user.email}</span>
            </li>
          ))}
        </ul>
      </div>
    );
  };
  
  // Render the main form
  return (
    <div className="collaborative-experiment-editor">
      <div className="editor-header">
        <h2>Edit Experiment: {formData.name}</h2>
        
        {renderCollaborators()}
        
        {lastSynced && (
          <div className="last-synced">
            Last synchronized: {new Date(lastSynced).toLocaleString()}
          </div>
        )}
      </div>
      
      <div className="form-group">
        <label htmlFor="name">Experiment Name</label>
        <input
          type="text"
          id="name"
          value={formData.name || ''}
          onChange={e => handleChange('name', e.target.value)}
        />
      </div>
      
      <div className="form-group">
        <label htmlFor="description">Description</label>
        <textarea
          id="description"
          value={formData.description || ''}
          onChange={e => handleChange('description', e.target.value)}
          rows={4}
        />
      </div>
      
      <div className="form-group">
        <label htmlFor="experimentType">Experiment Type</label>
        <select
          id="experimentType"
          value={formData.experimentType || ''}
          onChange={e => handleChange('experimentType', e.target.value)}
        >
          <option value="">Select experiment type...</option>
          <option value="cryopreservation">Cryopreservation</option>
          <option value="vitrification">Vitrification</option>
          <option value="freezing">Freezing</option>
          <option value="thawing">Thawing</option>
          <option value="storage">Storage</option>
          <option value="other">Other</option>
        </select>
      </div>
      
      <div className="form-group">
        <label htmlFor="status">Status</label>
        <select
          id="status"
          value={formData.status || ''}
          onChange={e => handleChange('status', e.target.value)}
        >
          <option value="">Select status...</option>
          <option value="planned">Planned</option>
          <option value="in_progress">In Progress</option>
          <option value="completed">Completed</option>
          <option value="aborted">Aborted</option>
          <option value="failed">Failed</option>
        </select>
      </div>
      
      <div className="form-group">
        <label htmlFor="startDate">Start Date</label>
        <input
          type="date"
          id="startDate"
          value={formData.startDate || ''}
          onChange={e => handleChange('startDate', e.target.value)}
        />
      </div>
      
      <div className="form-group">
        <label htmlFor="endDate">End Date</label>
        <input
          type="date"
          id="endDate"
          value={formData.endDate || ''}
          onChange={e => handleChange('endDate', e.target.value)}
        />
      </div>
      
      <div className="form-group">
        <label htmlFor="researcher">Researcher</label>
        <input
          type="text"
          id="researcher"
          value={formData.researcher || ''}
          onChange={e => handleChange('researcher', e.target.value)}
        />
      </div>
      
      <div className="form-group">
        <label htmlFor="notes">Notes</label>
        <textarea
          id="notes"
          value={formData.notes || ''}
          onChange={e => handleChange('notes', e.target.value)}
          rows={6}
        />
      </div>
      
      {/* Display local changes indicator */}
      {isDirty && (
        <div className="local-changes-indicator">
          <span className="dot"></span> You have unsaved changes
        </div>
      )}
      
      <div className="form-actions">
        <button
          className="secondary"
          onClick={() => setFormData(experiment || {})}
          disabled={isSaving || !isDirty}
        >
          Discard Changes
        </button>
        
        <button
          className="primary"
          onClick={handleSave}
          disabled={isSaving || !isDirty}
        >
          {isSaving ? 'Saving...' : 'Save Changes'}
        </button>
      </div>
    </div>
  );
};

export default CollaborativeExperimentEditor;
/**
 * Collaborative Experiment Editor Wrapper
 * 
 * This component provides a wrapper around the CollaborativeExperimentEditor component
 * to make it easy to integrate into existing pages. It handles conditional rendering
 * based on whether Convex is enabled.
 */
import React, { useState, useEffect } from 'react';
import CollaborativeExperimentEditor from '../../../components/CollaborativeExperimentEditor';
import ExperimentDetail from './ExperimentDetail';
import { useMutation } from 'convex/react';
import { api } from '../../../convex/_generated/api';

/**
 * A wrapper component that conditionally renders either the collaborative editor
 * or the standard experiment detail component based on configuration
 * 
 * @param {Object} props - Component props
 * @param {Object} props.experiment - The experiment data to display or edit
 * @param {Function} props.onSave - Function to call when saving changes
 * @param {boolean} props.readOnly - Whether the experiment should be read-only
 * @returns {React.ReactElement} Rendered component
 */
export default function CollaborativeExperimentEditorWrapper({ 
  experiment, 
  onSave,
  readOnly = false
}) {
  // State for tracking whether Convex is enabled
  const [isConvexEnabled, setIsConvexEnabled] = useState(false);
  const [isEditing, setIsEditing] = useState(false);
  const [isSaving, setIsSaving] = useState(false);
  
  // Convex mutations
  const updateExperiment = useMutation(api.experiments.enhanced_experiments.updateEnhancedExperiment);
  
  // Check if Convex is enabled on mount
  useEffect(() => {
    setIsConvexEnabled(process.env.NEXT_PUBLIC_USE_CONVEX === 'true');
  }, []);
  
  // Handle save callback
  const handleSave = async (data) => {
    setIsSaving(true);
    
    try {
      // If we're using Convex, the data is already saved in real-time
      // But we might need to perform additional updates
      if (isConvexEnabled && experiment._id) {
        // Make any final updates needed
        const updates = {
          // Add any fields that might need to be saved but aren't part of the collaborative editor
        };
        
        if (Object.keys(updates).length > 0) {
          await updateExperiment({ 
            experimentId: experiment._id,
            update: updates
          });
        }
      }
      
      setIsEditing(false);
      
      // Call the parent's onSave callback if provided
      if (onSave) {
        onSave(data);
      }
    } catch (err) {
      console.error('Error saving experiment:', err);
      // You could add error handling UI here
    } finally {
      setIsSaving(false);
    }
  };
  
  // If editing mode is not forced and experiment is read-only, show regular component
  if (readOnly || !isEditing) {
    return (
      <div className="experiment-wrapper">
        <ExperimentDetail 
          experiment={experiment} 
        />
        
        {!readOnly && (
          <div className="mt-6 flex justify-end">
            <button
              onClick={() => setIsEditing(true)}
              className="inline-flex items-center justify-center rounded-md bg-primary px-4 py-2 text-sm font-medium text-primary-foreground shadow transition-colors hover:bg-primary/90"
            >
              {isConvexEnabled ? 'Edit Collaboratively' : 'Edit Experiment'}
            </button>
          </div>
        )}
      </div>
    );
  }
  
  // If Convex is not enabled, we should redirect to the regular edit page
  if (!isConvexEnabled) {
    return (
      <div className="p-6 bg-yellow-50 border border-yellow-200 rounded-md">
        <h3 className="text-lg font-medium text-yellow-800">Real-time Collaboration Not Available</h3>
        <p className="mt-2 text-yellow-700">
          Real-time collaboration features require Convex to be enabled.
          Please edit this experiment using the standard editor instead.
        </p>
        <div className="mt-4">
          <button
            onClick={() => setIsEditing(false)}
            className="inline-flex items-center justify-center rounded-md bg-yellow-100 px-4 py-2 text-sm font-medium text-yellow-800 shadow transition-colors hover:bg-yellow-200"
          >
            Return to View Mode
          </button>
        </div>
      </div>
    );
  }
  
  // If Convex is enabled and we're in edit mode, show the collaborative editor
  return (
    <div className="experiment-collaborative-wrapper">
      {experiment && experiment._id && (
        <>
          <CollaborativeExperimentEditor
            experimentId={experiment._id}
            onSave={handleSave}
          />
          
          <div className="mt-6 flex justify-end">
            <button
              onClick={() => setIsEditing(false)}
              className="mr-2 inline-flex items-center justify-center rounded-md bg-gray-200 px-4 py-2 text-sm font-medium text-gray-800 shadow transition-colors hover:bg-gray-300"
              disabled={isSaving}
            >
              Cancel
            </button>
            
            <button
              onClick={() => handleSave(experiment)}
              className="inline-flex items-center justify-center rounded-md bg-primary px-4 py-2 text-sm font-medium text-primary-foreground shadow transition-colors hover:bg-primary/90"
              disabled={isSaving}
            >
              {isSaving ? 'Saving...' : 'Save Changes'}
            </button>
          </div>
        </>
      )}
    </div>
  );
}
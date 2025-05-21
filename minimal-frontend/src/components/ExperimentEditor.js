/**
 * Experiment Editor Component
 * 
 * This component provides editing for experiment data.
 * Works with both Convex and standard API modes.
 */
import React, { useState, useEffect } from 'react';

const ExperimentEditor = ({ 
  experimentId, 
  initialData, 
  onSave, 
  mode = 'standard' 
}) => {
  const [formData, setFormData] = useState({
    name: '',
    description: '',
    experimentType: '',
    status: '',
    startDate: '',
    endDate: '',
    researcher: '',
    tissueType: '',
    cryoprotectant: '',
    concentration: '',
    temperature: '',
    duration: '',
    notes: ''
  });
  const [isSaving, setIsSaving] = useState(false);
  const [isDirty, setIsDirty] = useState(false);
  const [error, setError] = useState(null);
  
  // Initialize form data
  useEffect(() => {
    if (initialData) {
      setFormData(prevData => ({
        ...prevData,
        ...initialData
      }));
    }
  }, [initialData]);
  
  // Handle form field changes
  const handleChange = (field, value) => {
    setFormData(prev => ({
      ...prev,
      [field]: value
    }));
    setIsDirty(true);
    setError(null);
  };
  
  // Validate form data
  const validateForm = () => {
    const errors = [];
    
    if (!formData.name?.trim()) {
      errors.push('Experiment name is required');
    }
    
    if (!formData.experimentType) {
      errors.push('Experiment type is required');
    }
    
    if (!formData.status) {
      errors.push('Status is required');
    }
    
    if (formData.endDate && formData.startDate && 
        new Date(formData.endDate) < new Date(formData.startDate)) {
      errors.push('End date must be after start date');
    }
    
    return errors;
  };
  
  // Save changes
  const handleSave = async () => {
    const validationErrors = validateForm();
    if (validationErrors.length > 0) {
      setError(validationErrors.join(', '));
      return;
    }
    
    if (!isDirty) return;
    
    setIsSaving(true);
    setError(null);
    
    try {
      // Call the onSave callback with the form data
      if (onSave) {
        await onSave(formData);
      }
      
      // Reset form state
      setIsDirty(false);
    } catch (error) {
      console.error('Error saving experiment:', error);
      setError('Error saving experiment: ' + error.message);
    } finally {
      setIsSaving(false);
    }
  };
  
  // Discard changes
  const handleDiscard = () => {
    if (initialData) {
      setFormData(prevData => ({
        ...prevData,
        ...initialData
      }));
    } else {
      setFormData({
        name: '',
        description: '',
        experimentType: '',
        status: '',
        startDate: '',
        endDate: '',
        researcher: '',
        tissueType: '',
        cryoprotectant: '',
        concentration: '',
        temperature: '',
        duration: '',
        notes: ''
      });
    }
    setIsDirty(false);
    setError(null);
  };
  
  return (
    <div className="experiment-editor">
      <style jsx>{`
        .experiment-editor {
          max-width: 800px;
          margin: 0 auto;
          padding: 2rem;
          background: white;
          border-radius: 12px;
          box-shadow: 0 4px 20px rgba(0, 0, 0, 0.1);
        }
        
        .editor-header {
          margin-bottom: 2rem;
          padding-bottom: 1rem;
          border-bottom: 2px solid #f0f0f0;
        }
        
        .editor-header h2 {
          margin: 0 0 0.5rem 0;
          color: #2c3e50;
          font-size: 1.8rem;
          font-weight: 600;
        }
        
        .form-section {
          margin-bottom: 2rem;
        }
        
        .form-section h3 {
          margin: 0 0 1rem 0;
          color: #34495e;
          font-size: 1.2rem;
          font-weight: 500;
          border-bottom: 1px solid #ecf0f1;
          padding-bottom: 0.5rem;
        }
        
        .form-grid {
          display: grid;
          grid-template-columns: 1fr;
          gap: 1.5rem;
        }
        
        @media (min-width: 768px) {
          .form-grid {
            grid-template-columns: 1fr 1fr;
          }
          
          .form-grid.single-column {
            grid-template-columns: 1fr;
          }
        }
        
        .form-group {
          display: flex;
          flex-direction: column;
        }
        
        .form-group.full-width {
          grid-column: 1 / -1;
        }
        
        .form-group label {
          margin-bottom: 0.5rem;
          font-weight: 500;
          color: #2c3e50;
          font-size: 0.9rem;
        }
        
        .form-group input,
        .form-group select,
        .form-group textarea {
          padding: 0.75rem;
          border: 2px solid #e1e8ed;
          border-radius: 8px;
          font-size: 1rem;
          transition: all 0.2s ease;
          background: white;
        }
        
        .form-group input:focus,
        .form-group select:focus,
        .form-group textarea:focus {
          outline: none;
          border-color: #3498db;
          box-shadow: 0 0 0 3px rgba(52, 152, 219, 0.1);
        }
        
        .form-group textarea {
          resize: vertical;
          min-height: 100px;
        }
        
        .form-group.large-textarea textarea {
          min-height: 150px;
        }
        
        .error-message {
          background: #ffebee;
          color: #c62828;
          padding: 1rem;
          border-radius: 8px;
          margin-bottom: 1rem;
          border-left: 4px solid #e53e3e;
        }
        
        .local-changes-indicator {
          display: flex;
          align-items: center;
          gap: 0.5rem;
          padding: 0.75rem 1rem;
          background: #fff3cd;
          border: 1px solid #ffeaa7;
          border-radius: 8px;
          margin-bottom: 1rem;
          color: #856404;
          font-size: 0.9rem;
        }
        
        .dot {
          width: 8px;
          height: 8px;
          background: #f39c12;
          border-radius: 50%;
          animation: pulse 2s infinite;
        }
        
        @keyframes pulse {
          0% { opacity: 1; }
          50% { opacity: 0.5; }
          100% { opacity: 1; }
        }
        
        .form-actions {
          display: flex;
          justify-content: flex-end;
          gap: 1rem;
          margin-top: 2rem;
          padding-top: 1rem;
          border-top: 2px solid #f0f0f0;
        }
        
        .form-actions button {
          padding: 0.75rem 1.5rem;
          border: none;
          border-radius: 8px;
          font-size: 1rem;
          font-weight: 500;
          cursor: pointer;
          transition: all 0.2s ease;
          min-width: 120px;
        }
        
        .form-actions button:disabled {
          opacity: 0.6;
          cursor: not-allowed;
        }
        
        .form-actions button.secondary {
          background: #ecf0f1;
          color: #2c3e50;
        }
        
        .form-actions button.secondary:hover:not(:disabled) {
          background: #d5dbdb;
        }
        
        .form-actions button.primary {
          background: #3498db;
          color: white;
        }
        
        .form-actions button.primary:hover:not(:disabled) {
          background: #2980b9;
        }
        
        .status-badge {
          display: inline-block;
          padding: 0.25rem 0.75rem;
          border-radius: 20px;
          font-size: 0.8rem;
          font-weight: 500;
          text-transform: uppercase;
          letter-spacing: 0.5px;
        }
        
        .status-planned { background: #e3f2fd; color: #1976d2; }
        .status-in_progress { background: #fff3e0; color: #f57c00; }
        .status-completed { background: #e8f5e8; color: #388e3c; }
        .status-aborted { background: #ffebee; color: #d32f2f; }
        .status-failed { background: #ffebee; color: #d32f2f; }
      `}</style>
      
      <div className="editor-header">
        <h2>
          {experimentId ? 'Edit Experiment' : 'Create New Experiment'}
          {formData.name && `: ${formData.name}`}
        </h2>
        {formData.status && (
          <div className={`status-badge status-${formData.status}`}>
            {formData.status.replace('_', ' ')}
          </div>
        )}
      </div>
      
      {error && (
        <div className="error-message">
          {error}
        </div>
      )}
      
      <div className="form-section">
        <h3>Basic Information</h3>
        <div className="form-grid">
          <div className="form-group">
            <label htmlFor="name">Experiment Name *</label>
            <input
              type="text"
              id="name"
              value={formData.name || ''}
              onChange={e => handleChange('name', e.target.value)}
              placeholder="Enter experiment name"
            />
          </div>
          
          <div className="form-group">
            <label htmlFor="experimentType">Experiment Type *</label>
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
              <option value="viability_study">Viability Study</option>
              <option value="recovery_study">Recovery Study</option>
              <option value="other">Other</option>
            </select>
          </div>
          
          <div className="form-group">
            <label htmlFor="status">Status *</label>
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
            <label htmlFor="researcher">Researcher</label>
            <input
              type="text"
              id="researcher"
              value={formData.researcher || ''}
              onChange={e => handleChange('researcher', e.target.value)}
              placeholder="Lead researcher name"
            />
          </div>
          
          <div className="form-group full-width">
            <label htmlFor="description">Description</label>
            <textarea
              id="description"
              value={formData.description || ''}
              onChange={e => handleChange('description', e.target.value)}
              placeholder="Describe the experiment objectives and methodology"
              rows={3}
            />
          </div>
        </div>
      </div>
      
      <div className="form-section">
        <h3>Timeline</h3>
        <div className="form-grid">
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
        </div>
      </div>
      
      <div className="form-section">
        <h3>Experimental Parameters</h3>
        <div className="form-grid">
          <div className="form-group">
            <label htmlFor="tissueType">Tissue/Cell Type</label>
            <input
              type="text"
              id="tissueType"
              value={formData.tissueType || ''}
              onChange={e => handleChange('tissueType', e.target.value)}
              placeholder="e.g., HeLa cells, liver tissue"
            />
          </div>
          
          <div className="form-group">
            <label htmlFor="cryoprotectant">Cryoprotectant</label>
            <input
              type="text"
              id="cryoprotectant"
              value={formData.cryoprotectant || ''}
              onChange={e => handleChange('cryoprotectant', e.target.value)}
              placeholder="e.g., DMSO, glycerol"
            />
          </div>
          
          <div className="form-group">
            <label htmlFor="concentration">Concentration</label>
            <input
              type="text"
              id="concentration"
              value={formData.concentration || ''}
              onChange={e => handleChange('concentration', e.target.value)}
              placeholder="e.g., 10%, 1M"
            />
          </div>
          
          <div className="form-group">
            <label htmlFor="temperature">Temperature</label>
            <input
              type="text"
              id="temperature"
              value={formData.temperature || ''}
              onChange={e => handleChange('temperature', e.target.value)}
              placeholder="e.g., -80°C, -196°C"
            />
          </div>
          
          <div className="form-group">
            <label htmlFor="duration">Duration</label>
            <input
              type="text"
              id="duration"
              value={formData.duration || ''}
              onChange={e => handleChange('duration', e.target.value)}
              placeholder="e.g., 24 hours, 1 week"
            />
          </div>
        </div>
      </div>
      
      <div className="form-section">
        <h3>Notes</h3>
        <div className="form-grid single-column">
          <div className="form-group large-textarea">
            <label htmlFor="notes">Additional Notes</label>
            <textarea
              id="notes"
              value={formData.notes || ''}
              onChange={e => handleChange('notes', e.target.value)}
              placeholder="Any additional observations, methods, or important notes about this experiment"
              rows={6}
            />
          </div>
        </div>
      </div>
      
      {isDirty && (
        <div className="local-changes-indicator">
          <span className="dot"></span> You have unsaved changes
        </div>
      )}
      
      <div className="form-actions">
        <button
          className="secondary"
          onClick={handleDiscard}
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

export default ExperimentEditor;
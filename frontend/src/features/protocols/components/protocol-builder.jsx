import React, { useState } from 'react';

/**
 * Protocol Builder Component
 * 
 * A form for creating and editing protocols
 */
export default function ProtocolBuilder({ initialData = {}, onSave, onCancel }) {
  const [protocolData, setProtocolData] = useState({
    name: '',
    description: '',
    objective: '',
    materials: [],
    steps: [],
    ...initialData
  });
  
  const [newStep, setNewStep] = useState({
    title: '',
    description: '',
    duration: '',
    temperature: '',
    notes: ''
  });
  
  const [newMaterial, setNewMaterial] = useState('');
  const [errors, setErrors] = useState({});
  
  // Handle basic input changes
  const handleChange = (e) => {
    const { name, value } = e.target;
    setProtocolData(prev => ({
      ...prev,
      [name]: value
    }));
    
    // Clear error for this field if it exists
    if (errors[name]) {
      setErrors(prev => {
        const newErrors = {...prev};
        delete newErrors[name];
        return newErrors;
      });
    }
  };
  
  // Handle step input changes
  const handleStepChange = (e) => {
    const { name, value } = e.target;
    setNewStep(prev => ({
      ...prev,
      [name]: value
    }));
  };
  
  // Add a new step to the protocol
  const addStep = () => {
    // Validate step has at least a title
    if (!newStep.title.trim()) {
      setErrors(prev => ({
        ...prev,
        stepTitle: 'Step title is required'
      }));
      return;
    }
    
    setProtocolData(prev => ({
      ...prev,
      steps: [...prev.steps, { ...newStep, id: Date.now().toString() }]
    }));
    
    // Clear the form
    setNewStep({
      title: '',
      description: '',
      duration: '',
      temperature: '',
      notes: ''
    });
    
    // Clear any step-related errors
    setErrors(prev => {
      const newErrors = {...prev};
      delete newErrors.stepTitle;
      return newErrors;
    });
  };
  
  // Remove a step from the protocol
  const removeStep = (stepId) => {
    setProtocolData(prev => ({
      ...prev,
      steps: prev.steps.filter(step => step.id !== stepId)
    }));
  };
  
  // Add a new material to the protocol
  const addMaterial = () => {
    if (!newMaterial.trim()) return;
    
    setProtocolData(prev => ({
      ...prev,
      materials: [...prev.materials, newMaterial.trim()]
    }));
    
    setNewMaterial('');
  };
  
  // Remove a material from the protocol
  const removeMaterial = (index) => {
    setProtocolData(prev => ({
      ...prev,
      materials: prev.materials.filter((_, i) => i !== index)
    }));
  };
  
  // Validate form before submission
  const validateForm = () => {
    const newErrors = {};
    
    if (!protocolData.name.trim()) {
      newErrors.name = 'Protocol name is required';
    }
    
    if (!protocolData.description.trim()) {
      newErrors.description = 'Description is required';
    }
    
    if (!protocolData.objective.trim()) {
      newErrors.objective = 'Objective is required';
    }
    
    if (protocolData.steps.length === 0) {
      newErrors.steps = 'At least one step is required';
    }
    
    setErrors(newErrors);
    return Object.keys(newErrors).length === 0;
  };
  
  // Handle form submission
  const handleSubmit = (e) => {
    e.preventDefault();
    
    if (!validateForm()) {
      return;
    }
    
    onSave(protocolData);
  };
  
  return (
    <div className="bg-white dark:bg-gray-800 rounded-md shadow p-6">
      <form onSubmit={handleSubmit}>
        <div className="space-y-6">
          {/* Basic Information */}
          <div>
            <h3 className="text-lg font-medium mb-4">Basic Information</h3>
            <div className="grid grid-cols-1 gap-4">
              <div>
                <label htmlFor="name" className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1">
                  Protocol Name
                </label>
                <input
                  type="text"
                  id="name"
                  name="name"
                  value={protocolData.name}
                  onChange={handleChange}
                  className="w-full px-3 py-2 border border-gray-300 dark:border-gray-600 rounded-md shadow-sm focus:outline-none focus:ring-primary focus:border-primary dark:bg-gray-700"
                  placeholder="Enter protocol name"
                />
                {errors.name && <p className="mt-1 text-sm text-red-600 dark:text-red-400">{errors.name}</p>}
              </div>
              
              <div>
                <label htmlFor="description" className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1">
                  Description
                </label>
                <textarea
                  id="description"
                  name="description"
                  value={protocolData.description}
                  onChange={handleChange}
                  rows={3}
                  className="w-full px-3 py-2 border border-gray-300 dark:border-gray-600 rounded-md shadow-sm focus:outline-none focus:ring-primary focus:border-primary dark:bg-gray-700"
                  placeholder="Describe the protocol"
                />
                {errors.description && <p className="mt-1 text-sm text-red-600 dark:text-red-400">{errors.description}</p>}
              </div>
              
              <div>
                <label htmlFor="objective" className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1">
                  Objective
                </label>
                <textarea
                  id="objective"
                  name="objective"
                  value={protocolData.objective}
                  onChange={handleChange}
                  rows={2}
                  className="w-full px-3 py-2 border border-gray-300 dark:border-gray-600 rounded-md shadow-sm focus:outline-none focus:ring-primary focus:border-primary dark:bg-gray-700"
                  placeholder="What is the goal of this protocol?"
                />
                {errors.objective && <p className="mt-1 text-sm text-red-600 dark:text-red-400">{errors.objective}</p>}
              </div>
            </div>
          </div>
          
          {/* Materials */}
          <div>
            <h3 className="text-lg font-medium mb-4">Materials</h3>
            <div className="flex space-x-2 mb-4">
              <input
                type="text"
                value={newMaterial}
                onChange={(e) => setNewMaterial(e.target.value)}
                className="flex-1 px-3 py-2 border border-gray-300 dark:border-gray-600 rounded-md shadow-sm focus:outline-none focus:ring-primary focus:border-primary dark:bg-gray-700"
                placeholder="Add a material"
              />
              <button
                type="button"
                onClick={addMaterial}
                className="px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700"
              >
                Add
              </button>
            </div>
            
            {protocolData.materials.length === 0 ? (
              <p className="text-sm text-gray-500 dark:text-gray-400 italic">No materials added yet</p>
            ) : (
              <ul className="space-y-2">
                {protocolData.materials.map((material, index) => (
                  <li key={index} className="flex justify-between items-center p-2 bg-gray-50 dark:bg-gray-700 rounded-md">
                    <span>{material}</span>
                    <button
                      type="button"
                      onClick={() => removeMaterial(index)}
                      className="text-red-600 dark:text-red-400 hover:text-red-800 dark:hover:text-red-300"
                    >
                      Remove
                    </button>
                  </li>
                ))}
              </ul>
            )}
          </div>
          
          {/* Steps */}
          <div>
            <h3 className="text-lg font-medium mb-4">Protocol Steps</h3>
            
            {/* Step form */}
            <div className="bg-gray-50 dark:bg-gray-700 p-4 rounded-md mb-4">
              <h4 className="font-medium mb-2">Add New Step</h4>
              <div className="grid grid-cols-1 md:grid-cols-2 gap-4 mb-4">
                <div>
                  <label htmlFor="stepTitle" className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1">
                    Step Title
                  </label>
                  <input
                    type="text"
                    id="stepTitle"
                    name="title"
                    value={newStep.title}
                    onChange={handleStepChange}
                    className="w-full px-3 py-2 border border-gray-300 dark:border-gray-600 rounded-md shadow-sm focus:outline-none focus:ring-primary focus:border-primary dark:bg-gray-700"
                    placeholder="E.g., Sample Preparation"
                  />
                  {errors.stepTitle && <p className="mt-1 text-sm text-red-600 dark:text-red-400">{errors.stepTitle}</p>}
                </div>
                
                <div className="grid grid-cols-2 gap-4">
                  <div>
                    <label htmlFor="stepDuration" className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1">
                      Duration
                    </label>
                    <input
                      type="text"
                      id="stepDuration"
                      name="duration"
                      value={newStep.duration}
                      onChange={handleStepChange}
                      className="w-full px-3 py-2 border border-gray-300 dark:border-gray-600 rounded-md shadow-sm focus:outline-none focus:ring-primary focus:border-primary dark:bg-gray-700"
                      placeholder="E.g., 30 min"
                    />
                  </div>
                  
                  <div>
                    <label htmlFor="stepTemperature" className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1">
                      Temperature
                    </label>
                    <input
                      type="text"
                      id="stepTemperature"
                      name="temperature"
                      value={newStep.temperature}
                      onChange={handleStepChange}
                      className="w-full px-3 py-2 border border-gray-300 dark:border-gray-600 rounded-md shadow-sm focus:outline-none focus:ring-primary focus:border-primary dark:bg-gray-700"
                      placeholder="E.g., 4°C"
                    />
                  </div>
                </div>
              </div>
              
              <div className="mb-4">
                <label htmlFor="stepDescription" className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1">
                  Step Description
                </label>
                <textarea
                  id="stepDescription"
                  name="description"
                  value={newStep.description}
                  onChange={handleStepChange}
                  rows={2}
                  className="w-full px-3 py-2 border border-gray-300 dark:border-gray-600 rounded-md shadow-sm focus:outline-none focus:ring-primary focus:border-primary dark:bg-gray-700"
                  placeholder="Describe what needs to be done in this step"
                />
              </div>
              
              <div>
                <label htmlFor="stepNotes" className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1">
                  Notes (Optional)
                </label>
                <textarea
                  id="stepNotes"
                  name="notes"
                  value={newStep.notes}
                  onChange={handleStepChange}
                  rows={2}
                  className="w-full px-3 py-2 border border-gray-300 dark:border-gray-600 rounded-md shadow-sm focus:outline-none focus:ring-primary focus:border-primary dark:bg-gray-700"
                  placeholder="Any additional notes or warnings"
                />
              </div>
              
              <div className="mt-4">
                <button
                  type="button"
                  onClick={addStep}
                  className="w-full py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700"
                >
                  Add Step
                </button>
              </div>
            </div>
            
            {/* Steps list */}
            {errors.steps && <p className="mb-2 text-sm text-red-600 dark:text-red-400">{errors.steps}</p>}
            
            {protocolData.steps.length === 0 ? (
              <p className="text-sm text-gray-500 dark:text-gray-400 italic">No steps added yet</p>
            ) : (
              <ol className="space-y-4">
                {protocolData.steps.map((step, index) => (
                  <li key={step.id} className="p-4 bg-gray-50 dark:bg-gray-700 rounded-md">
                    <div className="flex justify-between items-start mb-2">
                      <div className="flex items-center">
                        <span className="bg-blue-600 text-white w-6 h-6 rounded-full flex items-center justify-center mr-2">
                          {index + 1}
                        </span>
                        <h4 className="font-medium">{step.title}</h4>
                      </div>
                      <button
                        type="button"
                        onClick={() => removeStep(step.id)}
                        className="text-red-600 dark:text-red-400 hover:text-red-800 dark:hover:text-red-300"
                      >
                        Remove
                      </button>
                    </div>
                    
                    <p className="mb-2 text-sm">{step.description}</p>
                    
                    <div className="flex flex-wrap gap-2 text-xs text-gray-500 dark:text-gray-400">
                      {step.duration && <span>Duration: {step.duration}</span>}
                      {step.temperature && <span>Temperature: {step.temperature}</span>}
                    </div>
                    
                    {step.notes && (
                      <div className="mt-2 p-2 bg-yellow-50 dark:bg-yellow-900/20 text-sm rounded">
                        <strong>Notes:</strong> {step.notes}
                      </div>
                    )}
                  </li>
                ))}
              </ol>
            )}
          </div>
          
          {/* Form Actions */}
          <div className="flex justify-end space-x-2 pt-4 border-t">
            <button
              type="button"
              onClick={onCancel}
              className="px-4 py-2 border border-gray-300 dark:border-gray-600 rounded-md shadow-sm text-sm font-medium text-gray-700 dark:text-gray-300 bg-white dark:bg-gray-800 hover:bg-gray-50 dark:hover:bg-gray-700"
            >
              Cancel
            </button>
            <button
              type="submit"
              className="px-4 py-2 bg-primary text-white rounded-md shadow-sm text-sm font-medium hover:bg-primary/90"
            >
              Save Protocol
            </button>
          </div>
        </div>
      </form>
    </div>
  );
}
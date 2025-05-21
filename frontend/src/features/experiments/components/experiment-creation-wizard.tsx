import React, { useState } from 'react';
import { useRouter } from 'next/navigation';
import { createExperiment } from '../actions/experiment-actions';
import { MOCK_PROTOCOLS } from '../../protocols/data/mock-protocols';

/**
 * Component for creating a new experiment with a step-by-step wizard interface
 */
export function ExperimentCreationWizard() {
  const router = useRouter();
  const [currentStep, setCurrentStep] = useState(1);
  const [isSubmitting, setIsSubmitting] = useState(false);
  const [experimentData, setExperimentData] = useState({
    title: '',
    description: '',
    protocol: '',
    cellType: '',
    freezingRate: '',
    storageTemperature: '',
    thawingMethod: '',
    cryoprotectants: [{ id: '', name: '', concentration: '' }]
  });
  
  const totalSteps = 4;
  
  // Handle input changes
  const handleChange = (e: React.ChangeEvent<HTMLInputElement | HTMLTextAreaElement | HTMLSelectElement>) => {
    const { name, value } = e.target;
    setExperimentData(prev => ({
      ...prev,
      [name]: value
    }));
  };
  
  // Handle cryoprotectant changes
  const handleCryoprotectantChange = (index: number, field: string, value: string) => {
    setExperimentData(prev => {
      const updatedCryoprotectants = [...prev.cryoprotectants];
      updatedCryoprotectants[index] = {
        ...updatedCryoprotectants[index],
        [field]: value
      };
      return {
        ...prev,
        cryoprotectants: updatedCryoprotectants
      };
    });
  };
  
  // Add a new cryoprotectant field
  const addCryoprotectant = () => {
    setExperimentData(prev => ({
      ...prev,
      cryoprotectants: [...prev.cryoprotectants, { id: '', name: '', concentration: '' }]
    }));
  };
  
  // Remove a cryoprotectant field
  const removeCryoprotectant = (index: number) => {
    if (experimentData.cryoprotectants.length > 1) {
      setExperimentData(prev => ({
        ...prev,
        cryoprotectants: prev.cryoprotectants.filter((_, i) => i !== index)
      }));
    }
  };
  
  // Handle protocol selection
  const handleProtocolSelect = (protocolId: string) => {
    setExperimentData(prev => ({
      ...prev,
      protocol: protocolId
    }));
  };
  
  // Go to next step
  const handleNextStep = () => {
    setCurrentStep(prev => Math.min(prev + 1, totalSteps));
  };
  
  // Go to previous step
  const handlePrevStep = () => {
    setCurrentStep(prev => Math.max(prev - 1, 1));
  };
  
  // Get protocol data for a given ID
  const getProtocolData = (protocolId: string) => {
    return MOCK_PROTOCOLS.find(p => p.id === protocolId);
  };
  
  // Submit the experiment
  const handleSubmit = async () => {
    try {
      setIsSubmitting(true);
      
      // Get protocol data for the selected protocol
      const protocolData = getProtocolData(experimentData.protocol);
      
      // Create a new experiment object
      const newExperiment = {
        ...experimentData,
        status: 'In Progress',
        protocol: {
          id: experimentData.protocol,
          name: protocolData?.name || 'Unknown Protocol',
          version: protocolData?.version || '1.0.0'
        },
        results: {
          viability: 'Pending',
          recovery: 'Pending',
          functionality: 'Pending',
          notes: ''
        }
      };
      
      // Call the server action to create the experiment
      const result = await createExperiment(newExperiment);
      
      // Navigate to the newly created experiment
      router.push(`/experiments/${result.id}`);
    } catch (error) {
      console.error('Error creating experiment:', error);
      // Show error message
    } finally {
      setIsSubmitting(false);
    }
  };
  
  // Render step 1: Basic Information
  const renderStep1 = () => (
    <div className="space-y-6">
      <div>
        <label htmlFor="title" className="block text-sm font-medium mb-2">
          Experiment Title <span className="text-red-500">*</span>
        </label>
        <input
          type="text"
          id="title"
          name="title"
          value={experimentData.title}
          onChange={handleChange}
          className="w-full rounded-md border px-3 py-2"
          placeholder="e.g., DMSO Concentration Optimization"
          required
        />
      </div>
      
      <div>
        <label htmlFor="description" className="block text-sm font-medium mb-2">
          Description <span className="text-red-500">*</span>
        </label>
        <textarea
          id="description"
          name="description"
          value={experimentData.description}
          onChange={handleChange}
          className="w-full rounded-md border px-3 py-2"
          placeholder="Brief description of the experiment's purpose"
          rows={4}
          required
        />
      </div>
      
      <div>
        <label htmlFor="cellType" className="block text-sm font-medium mb-2">
          Cell Type <span className="text-red-500">*</span>
        </label>
        <input
          type="text"
          id="cellType"
          name="cellType"
          value={experimentData.cellType}
          onChange={handleChange}
          className="w-full rounded-md border px-3 py-2"
          placeholder="e.g., Human HeLa cells"
          required
        />
      </div>
    </div>
  );
  
  // Render step 2: Protocol Selection
  const renderStep2 = () => (
    <div className="space-y-6">
      <h3 className="text-lg font-medium mb-4">Select a Protocol</h3>
      
      <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
        {MOCK_PROTOCOLS.map(protocol => (
          <div 
            key={protocol.id}
            className={`border rounded-lg p-4 cursor-pointer ${
              experimentData.protocol === protocol.id 
                ? 'border-primary bg-primary/5' 
                : 'border-muted hover:border-gray-400'
            }`}
            onClick={() => handleProtocolSelect(protocol.id)}
          >
            <div className="flex items-start gap-2">
              <input 
                type="radio" 
                checked={experimentData.protocol === protocol.id}
                onChange={() => {}}
                className="mt-1"
              />
              <div>
                <h4 className="font-medium">{protocol.name}</h4>
                <p className="text-sm text-muted-foreground">{protocol.description}</p>
                <div className="flex items-center gap-2 mt-2 text-sm">
                  <span className="bg-muted px-2 py-0.5 rounded-full text-xs">v{protocol.version}</span>
                  <span className="text-muted-foreground">{protocol.duration}</span>
                </div>
              </div>
            </div>
          </div>
        ))}
      </div>
      
      {experimentData.protocol && (
        <div className="mt-6 p-4 border rounded-lg bg-muted/10">
          <h4 className="font-medium mb-2">Selected Protocol Details</h4>
          {(() => {
            const protocol = getProtocolData(experimentData.protocol);
            return protocol ? (
              <div className="space-y-3">
                <div className="grid grid-cols-1 md:grid-cols-2 gap-x-4 gap-y-2 text-sm">
                  <div>
                    <span className="font-medium">Steps:</span> {protocol.steps.length}
                  </div>
                  <div>
                    <span className="font-medium">Duration:</span> {protocol.duration}
                  </div>
                  <div>
                    <span className="font-medium">Compatible Cell Types:</span>
                    <div className="mt-1">
                      <ul className="list-disc list-inside">
                        {protocol.cell_types.map((cell, index) => (
                          <li key={index}>{cell}</li>
                        ))}
                      </ul>
                    </div>
                  </div>
                  <div>
                    <span className="font-medium">Default Cryoprotectants:</span>
                    <div className="mt-1">
                      <ul className="list-disc list-inside">
                        {protocol.cryoprotectants.map((cp, index) => (
                          <li key={index}>{cp.name} ({cp.concentration})</li>
                        ))}
                      </ul>
                    </div>
                  </div>
                </div>
              </div>
            ) : null;
          })()}
        </div>
      )}
    </div>
  );
  
  // Render step 3: Cryopreservation Parameters
  const renderStep3 = () => (
    <div className="space-y-6">
      <div>
        <label htmlFor="freezingRate" className="block text-sm font-medium mb-2">
          Freezing Rate <span className="text-red-500">*</span>
        </label>
        <input
          type="text"
          id="freezingRate"
          name="freezingRate"
          value={experimentData.freezingRate}
          onChange={handleChange}
          className="w-full rounded-md border px-3 py-2"
          placeholder="e.g., -1°C/min controlled rate"
          required
        />
      </div>
      
      <div>
        <label htmlFor="storageTemperature" className="block text-sm font-medium mb-2">
          Storage Temperature <span className="text-red-500">*</span>
        </label>
        <input
          type="text"
          id="storageTemperature"
          name="storageTemperature"
          value={experimentData.storageTemperature}
          onChange={handleChange}
          className="w-full rounded-md border px-3 py-2"
          placeholder="e.g., -196°C (liquid nitrogen)"
          required
        />
      </div>
      
      <div>
        <label htmlFor="thawingMethod" className="block text-sm font-medium mb-2">
          Thawing Method <span className="text-red-500">*</span>
        </label>
        <input
          type="text"
          id="thawingMethod"
          name="thawingMethod"
          value={experimentData.thawingMethod}
          onChange={handleChange}
          className="w-full rounded-md border px-3 py-2"
          placeholder="e.g., Rapid thawing in 37°C water bath"
          required
        />
      </div>
      
      <div>
        <div className="flex justify-between items-center mb-2">
          <label className="block text-sm font-medium">
            Cryoprotectants <span className="text-red-500">*</span>
          </label>
          <button
            type="button"
            onClick={addCryoprotectant}
            className="text-sm text-primary hover:underline"
          >
            + Add Cryoprotectant
          </button>
        </div>
        
        {experimentData.cryoprotectants.map((cp, index) => (
          <div key={index} className="flex gap-3 mb-3 items-start">
            <div className="flex-1">
              <input
                type="text"
                value={cp.name}
                onChange={(e) => handleCryoprotectantChange(index, 'name', e.target.value)}
                className="w-full rounded-md border px-3 py-2"
                placeholder="Name (e.g., DMSO)"
                required
              />
            </div>
            <div className="w-1/3">
              <input
                type="text"
                value={cp.concentration}
                onChange={(e) => handleCryoprotectantChange(index, 'concentration', e.target.value)}
                className="w-full rounded-md border px-3 py-2"
                placeholder="Concentration (e.g., 10%)"
                required
              />
            </div>
            {experimentData.cryoprotectants.length > 1 && (
              <button
                type="button"
                onClick={() => removeCryoprotectant(index)}
                className="p-2 text-red-500 hover:text-red-700"
              >
                &times;
              </button>
            )}
          </div>
        ))}
      </div>
    </div>
  );
  
  // Render step 4: Review
  const renderStep4 = () => (
    <div className="space-y-6">
      <h3 className="text-lg font-medium mb-4">Review Experiment Details</h3>
      
      <div className="bg-muted/10 border rounded-lg p-4">
        <dl className="grid grid-cols-1 md:grid-cols-2 gap-x-4 gap-y-3">
          <div>
            <dt className="text-sm font-medium text-muted-foreground">Title</dt>
            <dd className="mt-1">{experimentData.title}</dd>
          </div>
          <div>
            <dt className="text-sm font-medium text-muted-foreground">Cell Type</dt>
            <dd className="mt-1">{experimentData.cellType}</dd>
          </div>
          <div className="md:col-span-2">
            <dt className="text-sm font-medium text-muted-foreground">Description</dt>
            <dd className="mt-1">{experimentData.description}</dd>
          </div>
          <div>
            <dt className="text-sm font-medium text-muted-foreground">Protocol</dt>
            <dd className="mt-1">
              {getProtocolData(experimentData.protocol)?.name || 'No protocol selected'}
            </dd>
          </div>
          <div>
            <dt className="text-sm font-medium text-muted-foreground">Freezing Rate</dt>
            <dd className="mt-1">{experimentData.freezingRate}</dd>
          </div>
          <div>
            <dt className="text-sm font-medium text-muted-foreground">Storage Temperature</dt>
            <dd className="mt-1">{experimentData.storageTemperature}</dd>
          </div>
          <div>
            <dt className="text-sm font-medium text-muted-foreground">Thawing Method</dt>
            <dd className="mt-1">{experimentData.thawingMethod}</dd>
          </div>
          <div className="md:col-span-2">
            <dt className="text-sm font-medium text-muted-foreground">Cryoprotectants</dt>
            <dd className="mt-1">
              <ul className="list-disc list-inside">
                {experimentData.cryoprotectants.map((cp, index) => (
                  <li key={index}>{cp.name} ({cp.concentration})</li>
                ))}
              </ul>
            </dd>
          </div>
        </dl>
      </div>
      
      <div className="p-4 bg-yellow-50 border border-yellow-200 rounded-lg">
        <h4 className="font-medium text-yellow-800 mb-2">Next Steps</h4>
        <p className="text-sm text-yellow-700">
          Once created, you'll be able to:
        </p>
        <ul className="list-disc list-inside text-sm text-yellow-700 mt-2">
          <li>Record experiment results and observations</li>
          <li>Upload attachments like images and data files</li>
          <li>Compare results with other experiments</li>
          <li>Generate comprehensive reports</li>
        </ul>
      </div>
    </div>
  );
  
  // Render the current step
  const renderStep = () => {
    switch (currentStep) {
      case 1:
        return renderStep1();
      case 2:
        return renderStep2();
      case 3:
        return renderStep3();
      case 4:
        return renderStep4();
      default:
        return null;
    }
  };
  
  // Determine if the current step is valid
  const isStepValid = () => {
    switch (currentStep) {
      case 1:
        return experimentData.title && experimentData.description && experimentData.cellType;
      case 2:
        return !!experimentData.protocol;
      case 3:
        return (
          experimentData.freezingRate &&
          experimentData.storageTemperature &&
          experimentData.thawingMethod &&
          experimentData.cryoprotectants.length > 0 &&
          experimentData.cryoprotectants.every(cp => cp.name && cp.concentration)
        );
      case 4:
        return true;
      default:
        return false;
    }
  };
  
  return (
    <div className="bg-card rounded-lg border shadow-sm p-6">
      {/* Progress Steps */}
      <div className="mb-8">
        <div className="flex items-center justify-between">
          {Array.from({ length: totalSteps }).map((_, index) => (
            <React.Fragment key={index}>
              {/* Step circle */}
              <div className="relative">
                <div className={`
                  w-10 h-10 rounded-full flex items-center justify-center text-sm font-medium
                  ${currentStep > index + 1 
                    ? 'bg-primary text-white' 
                    : currentStep === index + 1 
                      ? 'bg-primary text-white' 
                      : 'bg-muted text-muted-foreground'
                  }
                `}>
                  {currentStep > index + 1 ? (
                    <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M5 13l4 4L19 7" />
                    </svg>
                  ) : (
                    index + 1
                  )}
                </div>
                <div className="absolute -bottom-6 left-1/2 transform -translate-x-1/2 w-max text-xs font-medium">
                  {index === 0 && 'Basic Info'}
                  {index === 1 && 'Protocol'}
                  {index === 2 && 'Parameters'}
                  {index === 3 && 'Review'}
                </div>
              </div>
              
              {/* Connector line */}
              {index < totalSteps - 1 && (
                <div className={`flex-1 h-0.5 ${
                  currentStep > index + 1 ? 'bg-primary' : 'bg-muted'
                }`} />
              )}
            </React.Fragment>
          ))}
        </div>
      </div>
      
      {/* Step content */}
      <div className="mt-12 pb-6">
        {renderStep()}
      </div>
      
      {/* Navigation buttons */}
      <div className="border-t pt-6 flex justify-between">
        <button
          type="button"
          onClick={handlePrevStep}
          className={`px-4 py-2 rounded-md text-sm font-medium ${
            currentStep === 1 
              ? 'text-muted-foreground cursor-not-allowed' 
              : 'border text-foreground hover:bg-muted'
          }`}
          disabled={currentStep === 1}
        >
          Previous
        </button>
        
        <div>
          {currentStep < totalSteps ? (
            <button
              type="button"
              onClick={handleNextStep}
              disabled={!isStepValid()}
              className={`px-4 py-2 rounded-md bg-primary text-white text-sm font-medium ${
                isStepValid() 
                  ? 'hover:bg-primary/90' 
                  : 'opacity-50 cursor-not-allowed'
              }`}
            >
              Next
            </button>
          ) : (
            <button
              type="button"
              onClick={handleSubmit}
              disabled={isSubmitting}
              className="px-4 py-2 rounded-md bg-primary text-white text-sm font-medium hover:bg-primary/90 flex items-center"
            >
              {isSubmitting ? (
                <>
                  <div className="mr-2 h-4 w-4 border-2 border-white border-t-transparent rounded-full animate-spin"></div>
                  Creating...
                </>
              ) : (
                'Create Experiment'
              )}
            </button>
          )}
        </div>
      </div>
    </div>
  );
}
import React from 'react';
import Link from 'next/link';
import Head from 'next/head';

export default function CreateProtocolPage() {
  // State for form fields
  const [name, setName] = React.useState('');
  const [description, setDescription] = React.useState('');
  const [cellTypes, setCellTypes] = React.useState([]);
  const [newCellType, setNewCellType] = React.useState('');
  const [freezingRate, setFreezingRate] = React.useState('');
  const [storageTemp, setStorageTemp] = React.useState('');
  const [thawingMethod, setThawingMethod] = React.useState('');
  const [selectedCryoprotectants, setSelectedCryoprotectants] = React.useState([]);
  const [steps, setSteps] = React.useState([
    { order: 1, title: '', description: '', duration: '', isCritical: false }
  ]);
  const [isTemplate, setIsTemplate] = React.useState(false);
  const [notes, setNotes] = React.useState('');
  
  // Mock data
  const availableCryoprotectants = [
    { id: '1', name: 'DMSO', chemical_formula: 'C2H6OS' },
    { id: '2', name: 'Glycerol', chemical_formula: 'C3H8O3' },
    { id: '3', name: 'Trehalose', chemical_formula: 'C12H22O11' },
    { id: '4', name: 'Propylene Glycol', chemical_formula: 'C3H8O2' },
    { id: '5', name: 'Ethylene Glycol', chemical_formula: 'C2H6O2' }
  ];
  
  const handleSubmit = (e) => {
    e.preventDefault();
    // This would typically save the protocol data to the backend
    console.log({
      name,
      description,
      cellTypes,
      freezingRate,
      storageTemp,
      thawingMethod,
      selectedCryoprotectants,
      steps,
      isTemplate,
      notes
    });
    
    // For now, just redirect back to protocols list
    window.location.href = '/protocols';
  };
  
  const toggleCryoprotectant = (id) => {
    if (selectedCryoprotectants.includes(id)) {
      setSelectedCryoprotectants(selectedCryoprotectants.filter(item => item !== id));
    } else {
      setSelectedCryoprotectants([...selectedCryoprotectants, id]);
    }
  };
  
  const addCellType = () => {
    if (newCellType.trim() !== '' && !cellTypes.includes(newCellType.trim())) {
      setCellTypes([...cellTypes, newCellType.trim()]);
      setNewCellType('');
    }
  };
  
  const removeCellType = (index) => {
    setCellTypes(cellTypes.filter((_, i) => i !== index));
  };
  
  const handleStepChange = (index, field, value) => {
    const updatedSteps = [...steps];
    updatedSteps[index] = { ...updatedSteps[index], [field]: value };
    setSteps(updatedSteps);
  };
  
  const addStep = () => {
    setSteps([...steps, { 
      order: steps.length + 1, 
      title: '', 
      description: '', 
      duration: '',
      isCritical: false
    }]);
  };
  
  const removeStep = (index) => {
    if (steps.length > 1) {
      const updatedSteps = steps.filter((_, i) => i !== index)
        .map((step, i) => ({ ...step, order: i + 1 }));
      setSteps(updatedSteps);
    }
  };
  
  const moveStep = (index, direction) => {
    if (
      (direction === 'up' && index === 0) || 
      (direction === 'down' && index === steps.length - 1)
    ) {
      return;
    }
    
    const updatedSteps = [...steps];
    const targetIndex = direction === 'up' ? index - 1 : index + 1;
    
    // Swap the steps
    [updatedSteps[index], updatedSteps[targetIndex]] = [updatedSteps[targetIndex], updatedSteps[index]];
    
    // Fix the order numbers
    updatedSteps.forEach((step, i) => {
      step.order = i + 1;
    });
    
    setSteps(updatedSteps);
  };

  return (
    <>
      <Head>
        <title>Create New Protocol | CryoProtect</title>
        <meta name="description" content="Create a new standardized protocol for cryopreservation experiments" />
      </Head>
    
      <div className="container mx-auto px-4 py-8">
        <div className="mb-6">
          <Link href="/protocols">
            <a className="text-primary hover:underline flex items-center">
              <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" className="mr-2">
                <path d="M19 12H5"/>
                <path d="M12 19l-7-7 7-7"/>
              </svg>
              Back to Protocols
            </a>
          </Link>
        </div>
        
        <div className="mb-8">
          <h1 className="text-3xl font-bold mb-2">Create New Protocol</h1>
          <p className="text-muted-foreground">
            Design a standardized protocol for reproducible cryopreservation procedures
          </p>
        </div>
        
        <form onSubmit={handleSubmit} className="max-w-4xl">
          <div className="bg-card rounded-lg border p-6 shadow-sm mb-6">
            <h2 className="text-xl font-semibold mb-4">Basic Information</h2>
            <div className="space-y-4">
              <div>
                <label htmlFor="name" className="block text-sm font-medium mb-1">
                  Protocol Name <span className="text-destructive">*</span>
                </label>
                <input
                  id="name"
                  type="text"
                  className="w-full px-3 py-2 border rounded-md"
                  placeholder="Enter a descriptive name for your protocol"
                  value={name}
                  onChange={(e) => setName(e.target.value)}
                  required
                />
              </div>
              
              <div>
                <label htmlFor="description" className="block text-sm font-medium mb-1">
                  Description <span className="text-destructive">*</span>
                </label>
                <textarea
                  id="description"
                  rows={3}
                  className="w-full px-3 py-2 border rounded-md"
                  placeholder="Describe the purpose and application of this protocol"
                  value={description}
                  onChange={(e) => setDescription(e.target.value)}
                  required
                />
              </div>
              
              <div className="flex items-center">
                <input
                  id="isTemplate"
                  type="checkbox"
                  className="h-4 w-4 mr-2"
                  checked={isTemplate}
                  onChange={(e) => setIsTemplate(e.target.checked)}
                />
                <label htmlFor="isTemplate" className="text-sm font-medium">
                  Make this protocol available as a template
                </label>
              </div>
            </div>
          </div>
          
          <div className="bg-card rounded-lg border p-6 shadow-sm mb-6">
            <h2 className="text-xl font-semibold mb-4">Compatible Cell Types</h2>
            <div className="space-y-4">
              <div className="flex flex-wrap gap-2">
                {cellTypes.map((type, index) => (
                  <div key={index} className="bg-muted px-3 py-1 rounded-full flex items-center">
                    <span className="text-sm">{type}</span>
                    <button 
                      type="button"
                      className="ml-2 text-muted-foreground hover:text-foreground"
                      onClick={() => removeCellType(index)}
                    >
                      <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
                        <path d="M18 6L6 18"></path>
                        <path d="M6 6l12 12"></path>
                      </svg>
                    </button>
                  </div>
                ))}
              </div>
              
              <div className="flex">
                <input
                  type="text"
                  className="w-full px-3 py-2 border rounded-l-md"
                  placeholder="Add a cell type..."
                  value={newCellType}
                  onChange={(e) => setNewCellType(e.target.value)}
                  onKeyPress={(e) => e.key === 'Enter' && (e.preventDefault(), addCellType())}
                />
                <button
                  type="button"
                  className="bg-primary text-primary-foreground px-4 rounded-r-md"
                  onClick={addCellType}
                >
                  Add
                </button>
              </div>
            </div>
          </div>
          
          <div className="bg-card rounded-lg border p-6 shadow-sm mb-6">
            <h2 className="text-xl font-semibold mb-4">Cryopreservation Parameters</h2>
            <div className="space-y-4">
              <div>
                <label className="block text-sm font-medium mb-3">
                  Cryoprotectants <span className="text-destructive">*</span>
                </label>
                <div className="grid grid-cols-1 md:grid-cols-2 gap-3">
                  {availableCryoprotectants.map(cp => (
                    <div 
                      key={cp.id} 
                      className={`flex items-center p-3 border rounded-md cursor-pointer ${
                        selectedCryoprotectants.includes(cp.id) ? 'bg-primary/10 border-primary' : ''
                      }`}
                      onClick={() => toggleCryoprotectant(cp.id)}
                    >
                      <input
                        type="checkbox"
                        className="h-4 w-4 mr-2"
                        checked={selectedCryoprotectants.includes(cp.id)}
                        onChange={() => {}}
                      />
                      <div>
                        <div className="font-medium">{cp.name}</div>
                        <div className="text-sm text-muted-foreground">{cp.chemical_formula}</div>
                      </div>
                    </div>
                  ))}
                </div>
                <div className="mt-2 flex justify-end">
                  <Link href="/molecules">
                    <a className="text-primary text-sm hover:underline">
                      Browse more molecules
                    </a>
                  </Link>
                </div>
              </div>
              
              <div>
                <label htmlFor="freezingRate" className="block text-sm font-medium mb-1">
                  Freezing Rate <span className="text-destructive">*</span>
                </label>
                <select
                  id="freezingRate"
                  className="w-full px-3 py-2 border rounded-md"
                  value={freezingRate}
                  onChange={(e) => setFreezingRate(e.target.value)}
                  required
                >
                  <option value="">Select freezing rate</option>
                  <option value="uncontrolled">Uncontrolled freezing</option>
                  <option value="-1C/min">-1°C/min (standard slow freezing)</option>
                  <option value="-2C/min">-2°C/min</option>
                  <option value="-5C/min">-5°C/min</option>
                  <option value="-10C/min">-10°C/min</option>
                  <option value="custom">Custom rate</option>
                </select>
              </div>
              
              <div>
                <label htmlFor="storageTemp" className="block text-sm font-medium mb-1">
                  Storage Temperature <span className="text-destructive">*</span>
                </label>
                <select
                  id="storageTemp"
                  className="w-full px-3 py-2 border rounded-md"
                  value={storageTemp}
                  onChange={(e) => setStorageTemp(e.target.value)}
                  required
                >
                  <option value="">Select storage temperature</option>
                  <option value="-80C">-80°C (ultra-low freezer)</option>
                  <option value="-130C">-130°C (vapor phase)</option>
                  <option value="-196C">-196°C (liquid nitrogen)</option>
                  <option value="custom">Custom temperature</option>
                </select>
              </div>
              
              <div>
                <label htmlFor="thawingMethod" className="block text-sm font-medium mb-1">
                  Thawing Method <span className="text-destructive">*</span>
                </label>
                <select
                  id="thawingMethod"
                  className="w-full px-3 py-2 border rounded-md"
                  value={thawingMethod}
                  onChange={(e) => setThawingMethod(e.target.value)}
                  required
                >
                  <option value="">Select thawing method</option>
                  <option value="37C water bath">Rapid thawing in 37°C water bath</option>
                  <option value="room temp">Room temperature thawing</option>
                  <option value="controlled">Controlled rate thawing</option>
                  <option value="custom">Custom method</option>
                </select>
              </div>
            </div>
          </div>
          
          <div className="bg-card rounded-lg border p-6 shadow-sm mb-6">
            <div className="flex justify-between items-center mb-4">
              <h2 className="text-xl font-semibold">Protocol Steps</h2>
              <button
                type="button"
                className="inline-flex items-center text-sm text-primary hover:underline"
                onClick={addStep}
              >
                <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" className="mr-1">
                  <path d="M12 5v14"></path>
                  <path d="M5 12h14"></path>
                </svg>
                Add Step
              </button>
            </div>
            
            <div className="space-y-6">
              {steps.map((step, index) => (
                <div key={index} className="rounded-md border p-4">
                  <div className="flex justify-between items-center mb-3">
                    <div className="font-medium flex items-center">
                      <span className="h-6 w-6 rounded-full bg-primary text-white text-sm flex items-center justify-center mr-2">
                        {step.order}
                      </span>
                      Step {step.order}
                    </div>
                    <div className="flex space-x-2">
                      <button
                        type="button"
                        className="text-muted-foreground hover:text-foreground"
                        onClick={() => moveStep(index, 'up')}
                        disabled={index === 0}
                      >
                        <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
                          <path d="M12 19V5"></path>
                          <path d="M5 12l7-7 7 7"></path>
                        </svg>
                      </button>
                      <button
                        type="button"
                        className="text-muted-foreground hover:text-foreground"
                        onClick={() => moveStep(index, 'down')}
                        disabled={index === steps.length - 1}
                      >
                        <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
                          <path d="M12 5v14"></path>
                          <path d="M19 12l-7 7-7-7"></path>
                        </svg>
                      </button>
                      <button
                        type="button"
                        className="text-destructive/70 hover:text-destructive"
                        onClick={() => removeStep(index)}
                      >
                        <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
                          <path d="M3 6h18"></path>
                          <path d="M19 6v14a2 2 0 0 1-2 2H7a2 2 0 0 1-2-2V6m3 0V4a2 2 0 0 1 2-2h4a2 2 0 0 1 2 2v2"></path>
                        </svg>
                      </button>
                    </div>
                  </div>
                  
                  <div className="space-y-3">
                    <div>
                      <label htmlFor={`step-${index}-title`} className="block text-sm font-medium mb-1">
                        Title <span className="text-destructive">*</span>
                      </label>
                      <input
                        id={`step-${index}-title`}
                        type="text"
                        className="w-full px-3 py-2 border rounded-md"
                        placeholder="Enter step title"
                        value={step.title}
                        onChange={(e) => handleStepChange(index, 'title', e.target.value)}
                        required
                      />
                    </div>
                    
                    <div>
                      <label htmlFor={`step-${index}-description`} className="block text-sm font-medium mb-1">
                        Description <span className="text-destructive">*</span>
                      </label>
                      <textarea
                        id={`step-${index}-description`}
                        rows={2}
                        className="w-full px-3 py-2 border rounded-md"
                        placeholder="Describe what to do in this step"
                        value={step.description}
                        onChange={(e) => handleStepChange(index, 'description', e.target.value)}
                        required
                      />
                    </div>
                    
                    <div className="flex space-x-4">
                      <div className="flex-1">
                        <label htmlFor={`step-${index}-duration`} className="block text-sm font-medium mb-1">
                          Duration
                        </label>
                        <input
                          id={`step-${index}-duration`}
                          type="text"
                          className="w-full px-3 py-2 border rounded-md"
                          placeholder="e.g., 10 minutes"
                          value={step.duration}
                          onChange={(e) => handleStepChange(index, 'duration', e.target.value)}
                        />
                      </div>
                      
                      <div className="flex items-end mb-1 pl-2">
                        <input
                          id={`step-${index}-critical`}
                          type="checkbox"
                          className="h-4 w-4 mr-2"
                          checked={step.isCritical}
                          onChange={(e) => handleStepChange(index, 'isCritical', e.target.checked)}
                        />
                        <label htmlFor={`step-${index}-critical`} className="text-sm font-medium">
                          Critical Step
                        </label>
                      </div>
                    </div>
                  </div>
                </div>
              ))}
            </div>
          </div>
          
          <div className="bg-card rounded-lg border p-6 shadow-sm mb-6">
            <h2 className="text-xl font-semibold mb-4">Additional Notes</h2>
            <div>
              <label htmlFor="notes" className="block text-sm font-medium mb-1">
                Protocol Notes
              </label>
              <textarea
                id="notes"
                rows={4}
                className="w-full px-3 py-2 border rounded-md"
                placeholder="Add any additional notes, tips, or warnings for this protocol"
                value={notes}
                onChange={(e) => setNotes(e.target.value)}
              />
            </div>
          </div>
          
          <div className="flex gap-3 justify-end">
            <Link href="/protocols">
              <a className="inline-flex items-center justify-center rounded-md bg-secondary px-4 py-2 text-sm font-medium text-secondary-foreground shadow transition-colors hover:bg-secondary/80">
                Cancel
              </a>
            </Link>
            <button 
              type="submit"
              className="inline-flex items-center justify-center rounded-md bg-primary px-4 py-2 text-sm font-medium text-primary-foreground shadow transition-colors hover:bg-primary/90"
            >
              Create Protocol
            </button>
          </div>
        </form>
      </div>
    </>
  );
}
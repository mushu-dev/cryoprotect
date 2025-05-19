import React from 'react';
import Link from 'next/link';
import Head from 'next/head';

export default function CreateExperimentPage() {
  // State for form fields
  const [title, setTitle] = React.useState('');
  const [description, setDescription] = React.useState('');
  const [protocolId, setProtocolId] = React.useState('');
  const [cellType, setCellType] = React.useState('');
  const [freezingRate, setFreezingRate] = React.useState('');
  const [storageTemp, setStorageTemp] = React.useState('');
  const [thawingMethod, setThawingMethod] = React.useState('');
  const [selectedCryoprotectants, setSelectedCryoprotectants] = React.useState([]);
  
  // Mock data
  const protocols = [
    { id: '1', name: 'Standard Cell Freezing', version: '1.2.0' },
    { id: '2', name: 'Vitrification Protocol', version: '2.0.1' },
    { id: '3', name: 'Plant Tissue Preservation', version: '1.5.0' }
  ];
  
  const availableCryoprotectants = [
    { id: '1', name: 'DMSO', chemical_formula: 'C2H6OS' },
    { id: '2', name: 'Glycerol', chemical_formula: 'C3H8O3' },
    { id: '3', name: 'Trehalose', chemical_formula: 'C12H22O11' },
    { id: '4', name: 'Propylene Glycol', chemical_formula: 'C3H8O2' },
    { id: '5', name: 'Ethylene Glycol', chemical_formula: 'C2H6O2' }
  ];
  
  const handleSubmit = (e) => {
    e.preventDefault();
    // This would typically save the experiment data to the backend
    console.log({
      title,
      description,
      protocolId,
      cellType,
      freezingRate,
      storageTemp,
      thawingMethod,
      selectedCryoprotectants
    });
    
    // For now, just redirect back to experiments list
    window.location.href = '/experiments';
  };
  
  const toggleCryoprotectant = (id) => {
    if (selectedCryoprotectants.includes(id)) {
      setSelectedCryoprotectants(selectedCryoprotectants.filter(item => item !== id));
    } else {
      setSelectedCryoprotectants([...selectedCryoprotectants, id]);
    }
  };

  return (
    <>
      <Head>
        <title>Create New Experiment | CryoProtect</title>
        <meta name="description" content="Design a new cryopreservation experiment with detailed protocols and parameters" />
      </Head>
    
      <div className="container mx-auto px-4 py-8">
        <div className="mb-6">
          <Link href="/experiments">
            <a className="text-primary hover:underline flex items-center">
              <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" className="mr-2">
                <path d="M19 12H5"/>
                <path d="M12 19l-7-7 7-7"/>
              </svg>
              Back to Experiments
            </a>
          </Link>
        </div>
        
        <div className="mb-8">
          <h1 className="text-3xl font-bold mb-2">Create New Experiment</h1>
          <p className="text-muted-foreground">
            Design a new cryopreservation experiment with detailed parameters
          </p>
        </div>
        
        <form onSubmit={handleSubmit} className="max-w-4xl">
          <div className="bg-card rounded-lg border p-6 shadow-sm mb-6">
            <h2 className="text-xl font-semibold mb-4">Basic Information</h2>
            <div className="space-y-4">
              <div>
                <label htmlFor="title" className="block text-sm font-medium mb-1">
                  Experiment Title <span className="text-destructive">*</span>
                </label>
                <input
                  id="title"
                  type="text"
                  className="w-full px-3 py-2 border rounded-md"
                  placeholder="Enter a descriptive title for your experiment"
                  value={title}
                  onChange={(e) => setTitle(e.target.value)}
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
                  placeholder="Describe the purpose and goals of this experiment"
                  value={description}
                  onChange={(e) => setDescription(e.target.value)}
                  required
                />
              </div>
              
              <div>
                <label htmlFor="protocol" className="block text-sm font-medium mb-1">
                  Protocol <span className="text-destructive">*</span>
                </label>
                <select
                  id="protocol"
                  className="w-full px-3 py-2 border rounded-md"
                  value={protocolId}
                  onChange={(e) => setProtocolId(e.target.value)}
                  required
                >
                  <option value="">Select a protocol</option>
                  {protocols.map(protocol => (
                    <option key={protocol.id} value={protocol.id}>
                      {protocol.name} (v{protocol.version})
                    </option>
                  ))}
                  <option value="custom">Create New Protocol</option>
                </select>
              </div>
              
              <div>
                <label htmlFor="cellType" className="block text-sm font-medium mb-1">
                  Cell Type <span className="text-destructive">*</span>
                </label>
                <input
                  id="cellType"
                  type="text"
                  className="w-full px-3 py-2 border rounded-md"
                  placeholder="e.g., Human dermal fibroblasts"
                  value={cellType}
                  onChange={(e) => setCellType(e.target.value)}
                  required
                />
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
          
          <div className="flex gap-3 justify-end">
            <Link href="/experiments">
              <a className="inline-flex items-center justify-center rounded-md bg-secondary px-4 py-2 text-sm font-medium text-secondary-foreground shadow transition-colors hover:bg-secondary/80">
                Cancel
              </a>
            </Link>
            <button 
              type="submit"
              className="inline-flex items-center justify-center rounded-md bg-primary px-4 py-2 text-sm font-medium text-primary-foreground shadow transition-colors hover:bg-primary/90"
            >
              Create Experiment
            </button>
          </div>
        </form>
      </div>
    </>
  );
}
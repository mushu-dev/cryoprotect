import React, { useState } from 'react';
import Link from 'next/link';
import Head from 'next/head';
import { useRouter } from 'next/router';

export default function CreateExperimentPage() {
  const router = useRouter();
  
  // Form state
  const [formData, setFormData] = useState({
    name: '',
    description: '',
    date: new Date().toISOString().split('T')[0], // Today's date
    protocol: '',
    cellType: '',
    temperature: '',
    concentration: '',
    freezingRate: '',
    storageTemp: '',
    thawingMethod: '',
    notes: '',
    cryoprotectants: []
  });
  
  // Validation state
  const [errors, setErrors] = useState({});
  const [isSubmitting, setIsSubmitting] = useState(false);
  const [submitSuccess, setSubmitSuccess] = useState(false);
  
  // Mock data
  const protocols = [
    { id: '1', name: 'Standard Cell Freezing', version: '2.1.0' },
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
  
  // Handle input changes
  const handleChange = (e) => {
    const { name, value } = e.target;
    setFormData({
      ...formData,
      [name]: value
    });
    
    // Clear error for this field if any
    if (errors[name]) {
      setErrors({
        ...errors,
        [name]: null
      });
    }
  };
  
  // Toggle cryoprotectant selection
  const toggleCryoprotectant = (id) => {
    const newCryoprotectants = formData.cryoprotectants.includes(id)
      ? formData.cryoprotectants.filter(item => item !== id)
      : [...formData.cryoprotectants, id];
    
    setFormData({
      ...formData,
      cryoprotectants: newCryoprotectants
    });
  };
  
  // Validate form
  const validateForm = () => {
    const newErrors = {};
    
    // Required fields
    if (!formData.name) newErrors.name = 'Experiment name is required';
    if (formData.name && formData.name.length < 3) newErrors.name = 'Name must be at least 3 characters';
    
    if (!formData.description) newErrors.description = 'Description is required';
    if (!formData.cellType) newErrors.cellType = 'Cell type is required';
    if (!formData.protocol) newErrors.protocol = 'Please select a protocol';
    
    // Validate numeric fields
    if (formData.temperature && isNaN(Number(formData.temperature))) {
      newErrors.temperature = 'Temperature must be a number';
    }
    
    if (formData.concentration && isNaN(Number(formData.concentration))) {
      newErrors.concentration = 'Concentration must be a number';
    }
    
    // Check if at least one cryoprotectant is selected
    if (formData.cryoprotectants.length === 0) {
      newErrors.cryoprotectants = 'Select at least one cryoprotectant';
    }
    
    return newErrors;
  };
  
  // Form submission
  const handleSubmit = async (e) => {
    e.preventDefault();
    
    // Validate form
    const formErrors = validateForm();
    
    if (Object.keys(formErrors).length > 0) {
      setErrors(formErrors);
      return;
    }
    
    setIsSubmitting(true);
    
    try {
      // This would typically be an API call
      console.log('Submitting experiment data:', formData);
      
      // Simulate API call
      await new Promise(resolve => setTimeout(resolve, 1000));
      
      setSubmitSuccess(true);
      
      // Redirect to experiments page after successful submission
      setTimeout(() => {
        router.push('/experiments');
      }, 1500);
      
    } catch (error) {
      console.error('Error submitting experiment:', error);
      setErrors({ form: 'Failed to create experiment. Please try again.' });
    } finally {
      setIsSubmitting(false);
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
          <div className="mt-2 p-3 bg-blue-50 text-blue-800 rounded-md text-sm">
            <div className="flex items-start">
              <svg className="h-5 w-5 text-blue-600 mr-2 mt-0.5" viewBox="0 0 20 20" fill="currentColor">
                <path fillRule="evenodd" d="M18 10a8 8 0 11-16 0 8 8 0 0116 0zm-7-4a1 1 0 11-2 0 1 1 0 012 0zM9 9a1 1 0 000 2v3a1 1 0 001 1h1a1 1 0 100-2v-3a1 1 0 00-1-1H9z" clipRule="evenodd" />
              </svg>
              <div>
                <strong>Enhanced Data Analysis Available</strong>
                <p className="mt-1">
                  After creating your experiment, you'll have access to new visualization tools, 
                  filtering options, and the ability to compare results across multiple experiments.
                </p>
              </div>
            </div>
          </div>
        </div>
        
        {submitSuccess && (
          <div className="mb-6 p-4 bg-green-50 border border-green-200 rounded-md" data-testid="success-message">
            <div className="flex items-center">
              <svg className="h-5 w-5 text-green-500 mr-2" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 20 20" fill="currentColor">
                <path fillRule="evenodd" d="M10 18a8 8 0 100-16 8 8 0 000 16zm3.707-9.293a1 1 0 00-1.414-1.414L9 10.586 7.707 9.293a1 1 0 00-1.414 1.414l2 2a1 1 0 001.414 0l4-4z" clipRule="evenodd" />
              </svg>
              <span className="text-green-800 font-medium">Experiment created successfully!</span>
            </div>
          </div>
        )}
        
        {errors.form && (
          <div className="mb-6 p-4 bg-red-50 border border-red-200 rounded-md" data-testid="form-error">
            <div className="flex items-center">
              <svg className="h-5 w-5 text-red-500 mr-2" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 20 20" fill="currentColor">
                <path fillRule="evenodd" d="M10 18a8 8 0 100-16 8 8 0 000 16zM8.707 7.293a1 1 0 00-1.414 1.414L8.586 10l-1.293 1.293a1 1 0 101.414 1.414L10 11.414l1.293 1.293a1 1 0 001.414-1.414L11.414 10l1.293-1.293a1 1 0 00-1.414-1.414L10 8.586 8.707 7.293z" clipRule="evenodd" />
              </svg>
              <span className="text-red-800">{errors.form}</span>
            </div>
          </div>
        )}
        
        <form onSubmit={handleSubmit} className="max-w-4xl">
          <div className="bg-card rounded-lg border p-6 shadow-sm mb-6">
            <h2 className="text-xl font-semibold mb-4">Basic Information</h2>
            <div className="space-y-4">
              <div>
                <label htmlFor="name" className="block text-sm font-medium mb-1">
                  Experiment Name <span className="text-red-500">*</span>
                </label>
                <input
                  id="name"
                  name="name"
                  type="text"
                  className={`w-full px-3 py-2 border rounded-md ${errors.name ? 'border-red-500' : ''}`}
                  placeholder="Enter a descriptive name for your experiment"
                  value={formData.name}
                  onChange={handleChange}
                />
                {errors.name && (
                  <p className="mt-1 text-sm text-red-600" data-testid="form-error">{errors.name}</p>
                )}
              </div>
              
              <div>
                <label htmlFor="description" className="block text-sm font-medium mb-1">
                  Description <span className="text-red-500">*</span>
                </label>
                <textarea
                  id="description"
                  name="description"
                  rows={3}
                  className={`w-full px-3 py-2 border rounded-md ${errors.description ? 'border-red-500' : ''}`}
                  placeholder="Describe the purpose and goals of this experiment"
                  value={formData.description}
                  onChange={handleChange}
                />
                {errors.description && (
                  <p className="mt-1 text-sm text-red-600" data-testid="form-error">{errors.description}</p>
                )}
              </div>
              
              <div>
                <label htmlFor="date" className="block text-sm font-medium mb-1">
                  Date
                </label>
                <input
                  id="date"
                  name="date"
                  type="date"
                  className="w-full px-3 py-2 border rounded-md"
                  value={formData.date}
                  onChange={handleChange}
                />
              </div>
              
              <div>
                <label htmlFor="protocol" className="block text-sm font-medium mb-1">
                  Protocol <span className="text-red-500">*</span>
                </label>
                <select
                  id="protocol"
                  name="protocol"
                  className={`w-full px-3 py-2 border rounded-md ${errors.protocol ? 'border-red-500' : ''}`}
                  value={formData.protocol}
                  onChange={handleChange}
                >
                  <option value="">Select a protocol</option>
                  {protocols.map(protocol => (
                    <option key={protocol.id} value={protocol.id}>
                      {protocol.name} (v{protocol.version})
                    </option>
                  ))}
                </select>
                {errors.protocol && (
                  <p className="mt-1 text-sm text-red-600" data-testid="form-error">{errors.protocol}</p>
                )}
              </div>
              
              <div>
                <label htmlFor="cellType" className="block text-sm font-medium mb-1">
                  Cell Type <span className="text-red-500">*</span>
                </label>
                <input
                  id="cellType"
                  name="cellType"
                  type="text"
                  className={`w-full px-3 py-2 border rounded-md ${errors.cellType ? 'border-red-500' : ''}`}
                  placeholder="e.g., Human dermal fibroblasts"
                  value={formData.cellType}
                  onChange={handleChange}
                />
                {errors.cellType && (
                  <p className="mt-1 text-sm text-red-600" data-testid="form-error">{errors.cellType}</p>
                )}
              </div>
            </div>
          </div>
          
          <div className="bg-card rounded-lg border p-6 shadow-sm mb-6">
            <h2 className="text-xl font-semibold mb-4">Cryopreservation Parameters</h2>
            <div className="space-y-4">
              <div>
                <label className="block text-sm font-medium mb-3">
                  Cryoprotectants <span className="text-red-500">*</span>
                </label>
                <div className="grid grid-cols-1 md:grid-cols-2 gap-3">
                  {availableCryoprotectants.map(cp => (
                    <div 
                      key={cp.id} 
                      className={`flex items-center p-3 border rounded-md cursor-pointer ${
                        formData.cryoprotectants.includes(cp.id) ? 'bg-primary/10 border-primary' : ''
                      } ${errors.cryoprotectants ? 'border-red-500' : ''}`}
                      onClick={() => toggleCryoprotectant(cp.id)}
                    >
                      <input
                        type="checkbox"
                        className="h-4 w-4 mr-2"
                        checked={formData.cryoprotectants.includes(cp.id)}
                        onChange={() => {}}
                      />
                      <div>
                        <div className="font-medium">{cp.name}</div>
                        <div className="text-sm text-muted-foreground">{cp.chemical_formula}</div>
                      </div>
                    </div>
                  ))}
                </div>
                {errors.cryoprotectants && (
                  <p className="mt-1 text-sm text-red-600" data-testid="form-error">{errors.cryoprotectants}</p>
                )}
              </div>
              
              <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                <div>
                  <label htmlFor="temperature" className="block text-sm font-medium mb-1">
                    Temperature (°C)
                  </label>
                  <input
                    id="temperature"
                    name="temperature"
                    type="number"
                    className={`w-full px-3 py-2 border rounded-md ${errors.temperature ? 'border-red-500' : ''}`}
                    placeholder="e.g., 25"
                    value={formData.temperature}
                    onChange={handleChange}
                  />
                  {errors.temperature && (
                    <p className="mt-1 text-sm text-red-600" data-testid="form-error">{errors.temperature}</p>
                  )}
                </div>
                
                <div>
                  <label htmlFor="concentration" className="block text-sm font-medium mb-1">
                    Concentration (%)
                  </label>
                  <input
                    id="concentration"
                    name="concentration"
                    type="number"
                    className={`w-full px-3 py-2 border rounded-md ${errors.concentration ? 'border-red-500' : ''}`}
                    placeholder="e.g., 10"
                    value={formData.concentration}
                    onChange={handleChange}
                  />
                  {errors.concentration && (
                    <p className="mt-1 text-sm text-red-600" data-testid="form-error">{errors.concentration}</p>
                  )}
                </div>
              </div>
              
              <div>
                <label htmlFor="freezingRate" className="block text-sm font-medium mb-1">
                  Freezing Rate
                </label>
                <select
                  id="freezingRate"
                  name="freezingRate"
                  className="w-full px-3 py-2 border rounded-md"
                  value={formData.freezingRate}
                  onChange={handleChange}
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
                  Storage Temperature
                </label>
                <select
                  id="storageTemp"
                  name="storageTemp"
                  className="w-full px-3 py-2 border rounded-md"
                  value={formData.storageTemp}
                  onChange={handleChange}
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
                  Thawing Method
                </label>
                <select
                  id="thawingMethod"
                  name="thawingMethod"
                  className="w-full px-3 py-2 border rounded-md"
                  value={formData.thawingMethod}
                  onChange={handleChange}
                >
                  <option value="">Select thawing method</option>
                  <option value="37C water bath">Rapid thawing in 37°C water bath</option>
                  <option value="room temp">Room temperature thawing</option>
                  <option value="controlled">Controlled rate thawing</option>
                  <option value="custom">Custom method</option>
                </select>
              </div>
              
              <div>
                <label htmlFor="notes" className="block text-sm font-medium mb-1">
                  Notes
                </label>
                <textarea
                  id="notes"
                  name="notes"
                  rows={3}
                  className="w-full px-3 py-2 border rounded-md"
                  placeholder="Add any additional notes or observations"
                  value={formData.notes}
                  onChange={handleChange}
                />
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
              disabled={isSubmitting}
            >
              {isSubmitting ? 'Creating...' : 'Create Experiment'}
            </button>
          </div>
        </form>
      </div>
    </>
  );
}
/**
 * CryoProtect Analyzer Web Interface
 * API Module - Handles communication with the Flask API
 */

const API = (function() {
  // Base URL for API requests
  const API_BASE_URL = '/api/v1';
  
  // Default request options
  const defaultOptions = {
    headers: {
      'Content-Type': 'application/json'
    }
  };
  
  /**
   * Make an API request
   * @param {string} endpoint - API endpoint
   * @param {Object} options - Request options
   * @returns {Promise} Promise that resolves with the response data
   */
  async function request(endpoint, options = {}) {
    // Merge default options with provided options
    const requestOptions = {
      ...defaultOptions,
      ...options,
      headers: {
        ...defaultOptions.headers,
        ...options.headers
      }
    };
    
    // Add authentication token if available
    const token = Auth.getToken();
    if (token) {
      requestOptions.headers['Authorization'] = `Bearer ${token}`;
    }
    
    try {
      const response = await fetch(`${API_BASE_URL}${endpoint}`, requestOptions);
      
      // Parse response JSON
      const data = await response.json();
      
      // Check if response is successful
      if (!response.ok) {
        throw new Error(data.message || `API request failed with status ${response.status}`);
      }
      
      return data;
    } catch (error) {
      console.error('API request error:', error);
      throw error;
    }
  }
  
  /**
   * Get all molecules
   * @returns {Promise} Promise that resolves with the molecules data
   */
  function getMolecules() {
    return request('/molecules');
  }
  
  /**
   * Get a specific molecule
   * @param {string} moleculeId - The molecule ID
   * @returns {Promise} Promise that resolves with the molecule data
   */
  function getMolecule(moleculeId) {
    return request(`/molecules/${moleculeId}`);
  }
  
  /**
   * Import a molecule from PubChem
   * @param {number} cid - PubChem Compound ID
   * @returns {Promise} Promise that resolves with the imported molecule data
   */
  function importMolecule(cid) {
    return request('/molecules', {
      method: 'POST',
      body: JSON.stringify({ cid: parseInt(cid) })
    });
  }
  
  /**
   * Get all mixtures
   * @returns {Promise} Promise that resolves with the mixtures data
   */
  function getMixtures() {
    return request('/mixtures');
  }
  
  /**
   * Get a specific mixture
   * @param {string} mixtureId - The mixture ID
   * @returns {Promise} Promise that resolves with the mixture data
   */
  function getMixture(mixtureId) {
    return request(`/mixtures/${mixtureId}`);
  }
  
  /**
   * Create a new mixture
   * @param {Object} mixtureData - The mixture data
   * @returns {Promise} Promise that resolves with the created mixture data
   */
  function createMixture(mixtureData) {
    return request('/mixtures', {
      method: 'POST',
      body: JSON.stringify(mixtureData)
    });
  }
  
  /**
   * Get predictions for a mixture
   * @param {string} mixtureId - The mixture ID
   * @returns {Promise} Promise that resolves with the predictions data
   */
  function getPredictions(mixtureId) {
    return request(`/mixtures/${mixtureId}/predictions`);
  }
  
  /**
   * Create a prediction for a mixture
   * @param {string} mixtureId - The mixture ID
   * @param {Object} predictionData - The prediction data
   * @returns {Promise} Promise that resolves with the created prediction data
   */
  function createPrediction(mixtureId, predictionData) {
    return request(`/mixtures/${mixtureId}/predictions`, {
      method: 'POST',
      body: JSON.stringify(predictionData)
    });
  }
  
  /**
   * Get experiments for a mixture
   * @param {string} mixtureId - The mixture ID
   * @returns {Promise} Promise that resolves with the experiments data
   */
  function getExperiments(mixtureId) {
    return request(`/mixtures/${mixtureId}/experiments`);
  }
  
  /**
   * Create an experiment for a mixture
   * @param {string} mixtureId - The mixture ID
   * @param {Object} experimentData - The experiment data
   * @returns {Promise} Promise that resolves with the created experiment data
   */
  function createExperiment(mixtureId, experimentData) {
    return request(`/mixtures/${mixtureId}/experiments`, {
      method: 'POST',
      body: JSON.stringify(experimentData)
    });
  }
  
  /**
   * Get a comparison between prediction and experiment for a mixture
   * @param {string} mixtureId - The mixture ID
   * @param {string} propertyName - The property name to compare
   * @returns {Promise} Promise that resolves with the comparison data
   */
  function getComparison(mixtureId, propertyName) {
    return request(`/mixtures/${mixtureId}/comparisons?property_name=${encodeURIComponent(propertyName)}`);
  }
  
  /**
   * Get all property types
   * @returns {Promise} Promise that resolves with the property types data
   */
  function getPropertyTypes() {
    // This is a mock function since the API doesn't have a direct endpoint for property types
    // In a real implementation, this would call an API endpoint
    return Promise.resolve([
      { name: 'Molecular Weight', units: 'g/mol', data_type: 'numeric' },
      { name: 'LogP', data_type: 'numeric' },
      { name: 'TPSA', units: 'Å²', data_type: 'numeric' },
      { name: 'H-Bond Donors', data_type: 'numeric' },
      { name: 'H-Bond Acceptors', data_type: 'numeric' },
      { name: 'Toxicity', data_type: 'text' },
      { name: 'Stability', data_type: 'text' },
      { name: 'Environmental Safety', data_type: 'text' },
      { name: 'Freezing Point', units: '°C', data_type: 'numeric' },
      { name: 'Viscosity', units: 'cP', data_type: 'numeric' },
      { name: 'Cell Viability', units: '%', data_type: 'numeric' },
      { name: 'Ice Crystal Formation', data_type: 'text' }
    ]);
  }
  
  /**
   * Get all calculation methods
   * @returns {Promise} Promise that resolves with the calculation methods data
   */
  function getCalculationMethods() {
    // This is a mock function since the API doesn't have a direct endpoint for calculation methods
    // In a real implementation, this would call an API endpoint
    return Promise.resolve([
      { name: 'PubChem Properties', description: 'Properties fetched directly from PubChem' },
      { name: 'CryoProtect Scoring', description: 'Scoring algorithm based on weighted properties' },
      { name: 'Experimental Validation', description: 'Results from laboratory experiments' },
      { name: 'Molecular Dynamics', description: 'Simulation of molecular interactions' },
      { name: 'QSAR Model', description: 'Quantitative Structure-Activity Relationship model' }
    ]);
  }
  
  // Public API
  return {
    getMolecules,
    getMolecule,
    importMolecule,
    getMixtures,
    getMixture,
    createMixture,
    getPredictions,
    createPrediction,
    getExperiments,
    createExperiment,
    getComparison,
    getPropertyTypes,
    getCalculationMethods
  };
})();
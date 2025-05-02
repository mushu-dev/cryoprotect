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
  /**
   * Make an API request with enhanced error handling
   * @param {string} endpoint - API endpoint
   * @param {Object} options - Request options
   * @param {Object} retryOptions - Retry options
   * @returns {Promise} Promise that resolves with the response data
   */
  async function request(endpoint, options = {}, retryOptions = { retries: 0, delay: 1000 }) {
    // Merge default options with provided options
    const requestOptions = {
      ...defaultOptions,
      ...options,
      headers: {
        ...defaultOptions.headers,
        ...options.headers
      }
    };
    
    // Add authentication token only for modifying requests (not for GET)
    const method = (requestOptions.method || 'GET').toUpperCase();
    const token = Auth.getToken();
    if (token && method !== 'GET') {
      requestOptions.headers['Authorization'] = `Bearer ${token}`;
    }
    
    // Add request ID for tracking
    const requestId = Date.now().toString(36) + Math.random().toString(36).substr(2);
    requestOptions.headers['X-Request-ID'] = requestId;
    
    try {
      // Attempt the fetch request
      const response = await fetch(`${API_BASE_URL}${endpoint}`, requestOptions);
      
      // Try to parse response as JSON
      let data;
      const contentType = response.headers.get('content-type');
      if (contentType && contentType.includes('application/json')) {
        data = await response.json();
      } else {
        // Handle non-JSON responses
        const text = await response.text();
        data = { message: text };
      }
      
      // Check if response is successful
      if (!response.ok) {
        // Create an enhanced error object with additional context
        const error = createApiError(response.status, data, endpoint);
        
        // Log the error with request details
        console.error(`API Error [${requestId}]:`, {
          endpoint,
          status: response.status,
          method: requestOptions.method || 'GET',
          error: error.message,
          details: error.details
        });
        
        throw error;
      }
      
      return data;
    } catch (error) {
      // Check if it's a network error and retries are available
      if (error.name === 'TypeError' && retryOptions.retries > 0) {
        console.warn(`Network error, retrying... (${retryOptions.retries} attempts left)`);
        
        // Wait before retrying
        await new Promise(resolve => setTimeout(resolve, retryOptions.delay));
        
        // Retry with one less retry attempt
        return request(endpoint, options, {
          retries: retryOptions.retries - 1,
          delay: retryOptions.delay * 1.5 // Exponential backoff
        });
      }
      
      // If it's already an ApiError, just rethrow it
      if (error.isApiError) {
        throw error;
      }
      
      // For other errors, create a generic API error
      const apiError = new Error(error.message || 'An unexpected error occurred');
      apiError.isApiError = true;
      apiError.originalError = error;
      apiError.errorType = 'NETWORK_ERROR';
      apiError.userMessage = 'Unable to connect to the server. Please check your internet connection and try again.';
      
      console.error('API request error:', apiError);
      throw apiError;
    }
  }
  
  /**
   * Create an API error with enhanced context
   * @param {number} status - HTTP status code
   * @param {Object} data - Response data
   * @param {string} endpoint - API endpoint
   * @returns {Error} Enhanced error object
   */
  function createApiError(status, data, endpoint) {
    const error = new Error(data.message || `API request failed with status ${status}`);
    error.isApiError = true;
    error.status = status;
    error.endpoint = endpoint;
    error.details = data.details || {};
    
    // Categorize errors for better handling
    if (status === 401) {
      error.errorType = 'AUTHENTICATION_ERROR';
      error.userMessage = 'Your session has expired. Please log in again.';
    } else if (status === 403) {
      error.errorType = 'AUTHORIZATION_ERROR';
      error.userMessage = 'You do not have permission to perform this action.';
    } else if (status === 404) {
      error.errorType = 'NOT_FOUND_ERROR';
      error.userMessage = 'The requested resource was not found.';
    } else if (status === 422) {
      error.errorType = 'VALIDATION_ERROR';
      error.userMessage = 'There was a problem with the data you submitted.';
      
      // Format validation errors for better display
      if (data.details && data.details.errors) {
        error.validationErrors = data.details.errors;
        error.userMessage = formatValidationErrors(data.details.errors);
      }
    } else if (status >= 500) {
      error.errorType = 'SERVER_ERROR';
      error.userMessage = 'The server encountered an error. Please try again later.';
    } else {
      error.errorType = 'UNKNOWN_ERROR';
      error.userMessage = 'An unexpected error occurred. Please try again.';
    }
    
    return error;
  }
  
  /**
   * Format validation errors into a user-friendly message
   * @param {Object} errors - Validation errors object
   * @returns {string} Formatted error message
   */
  function formatValidationErrors(errors) {
    if (!errors || typeof errors !== 'object') {
      return 'Invalid data submitted.';
    }
    
    const errorMessages = [];
    
    for (const field in errors) {
      if (Object.prototype.hasOwnProperty.call(errors, field)) {
        const fieldErrors = errors[field];
        const fieldName = formatFieldName(field);
        
        if (Array.isArray(fieldErrors)) {
          errorMessages.push(`${fieldName}: ${fieldErrors.join(', ')}`);
        } else {
          errorMessages.push(`${fieldName}: ${fieldErrors}`);
        }
      }
    }
    
    return errorMessages.length > 0
      ? `Please correct the following: ${errorMessages.join('; ')}`
      : 'There was a problem with the data you submitted.';
  }
  
  /**
   * Format a field name for display
   * @param {string} field - Field name
   * @returns {string} Formatted field name
   */
  function formatFieldName(field) {
    return field
      .replace(/_/g, ' ')
      .replace(/([A-Z])/g, ' $1')
      .replace(/^./, str => str.toUpperCase())
      .trim();
  }
  
  /**
   * Get molecules with pagination support
   * @param {Object} options - Pagination and filtering options
   * @param {number} options.limit - Maximum number of items to return (default: 20)
   * @param {number} options.offset - Number of items to skip (default: 0)
   * @param {string} options.search - Search term for filtering molecules
   * @param {string} options.sort - Field to sort by
   * @param {string} options.order - Sort order ('asc' or 'desc')
   * @returns {Promise} Promise that resolves with the molecules data and pagination metadata
   */
  function getMolecules(options = {}) {
    // Build query parameters
    const params = new URLSearchParams();
    
    // Add pagination parameters
    if (options.limit !== undefined) {
      params.append('limit', options.limit);
    }
    
    if (options.offset !== undefined) {
      params.append('offset', options.offset);
    }
    
    // Add search parameter if provided
    if (options.search) {
      params.append('search', options.search);
    }
    
    // Add sort parameters if provided
    if (options.sort) {
      params.append('sort', options.sort);
    }
    
    if (options.order) {
      params.append('order', options.order);
    }
    
    // Build the URL with query parameters
    const url = params.toString() ? `/molecules?${params.toString()}` : '/molecules';
    
    return request(url);
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
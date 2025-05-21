/**
 * API utility functions for fetching data from the CryoProtect backend
 * Supports both REST API and Convex
 * Enhanced with circuit breaker pattern, retry mechanisms, and fallbacks
 */

// Convex client import disabled for minimal deployment
// import { convex } from '../src/convex/client';

// Base API URLs
const API_ENDPOINTS = {
  HEROKU_API: 'https://cryoprotect-8030e4025428.herokuapp.com/api',
  RDKIT_API: 'https://cryoprotect-rdkit.fly.dev',
  CONVEX: 'https://primary-meerkat-478.convex.cloud'
};

// Disable Convex for minimal deployment
const USE_CONVEX = false;

// Fallback mock data in case the API is not available
const MOCK_DATA = {
  molecules: [
    {
      id: 1,
      name: 'Glycerol',
      formula: 'C3H8O3',
      pubchem_cid: '753',
      smiles: 'C(C(CO)O)O',
      molecular_weight: 92.09,
      is_cryoprotectant: true,
      description: 'A common cryoprotectant used in various applications.'
    },
    {
      id: 2,
      name: 'Dimethyl Sulfoxide (DMSO)',
      formula: 'C2H6OS',
      pubchem_cid: '679',
      smiles: 'CS(=O)C',
      molecular_weight: 78.13,
      is_cryoprotectant: true,
      description: 'A widely used penetrating cryoprotectant.'
    },
    {
      id: 3,
      name: 'Ethylene Glycol',
      formula: 'C2H6O2',
      pubchem_cid: '174',
      smiles: 'OCCO',
      molecular_weight: 62.07,
      is_cryoprotectant: true,
      description: 'Used in cryopreservation of embryos and tissues.'
    },
    {
      id: 4,
      name: 'Propylene Glycol',
      formula: 'C3H8O2',
      pubchem_cid: '1030',
      smiles: 'CC(O)CO',
      molecular_weight: 76.09,
      is_cryoprotectant: true,
      description: 'Used as a cryoprotectant for various biological materials.'
    }
  ],
  mixtures: [
    {
      id: 1,
      name: 'VS55',
      description: 'A cryoprotectant mixture used for organ preservation.',
      components: [
        {
          molecule: {
            id: 1,
            name: 'Glycerol',
            formula: 'C3H8O3'
          },
          concentration: 15,
          concentration_unit: '%',
          role: 'Primary cryoprotectant'
        },
        {
          molecule: {
            id: 2,
            name: 'Dimethyl Sulfoxide (DMSO)',
            formula: 'C2H6OS'
          },
          concentration: 20,
          concentration_unit: '%',
          role: 'Penetrating cryoprotectant'
        },
        {
          molecule: {
            id: 3,
            name: 'Ethylene Glycol',
            formula: 'C2H6O2'
          },
          concentration: 20,
          concentration_unit: '%',
          role: 'Penetrating cryoprotectant'
        }
      ]
    },
    {
      id: 2,
      name: 'Standard Cell Freezing Medium',
      description: 'Common mixture for cell line preservation.',
      components: [
        {
          molecule: {
            id: 2,
            name: 'Dimethyl Sulfoxide (DMSO)',
            formula: 'C2H6OS'
          },
          concentration: 10,
          concentration_unit: '%',
          role: 'Primary cryoprotectant'
        }
      ]
    }
  ],
  // Mock SVG data for molecule depictions (minimal example)
  molecular_depictions: {
    'C(C(CO)O)O': `<svg version="1.1" id="mol-glycerol" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" x="0px" y="0px" width="300px" height="200px" viewBox="0 0 300 200" enable-background="new 0 0 300 200">
      <rect fill="#FFFFFF" width="300" height="200"/>
      <g transform="translate(10,10) scale(2.8)">
        <path fill="none" stroke="#000000" stroke-width="1.5" d="M 50,34.6 L 60,40 L 70,34.6"/>
        <path fill="none" stroke="#000000" stroke-width="1.5" d="M 60,40 L 60,50"/>
        <path fill="none" stroke="#000000" stroke-width="1.5" d="M 60,50 L 70,55.4"/>
        <path fill="none" stroke="#000000" stroke-width="1.5" d="M 60,50 L 50,55.4"/>
        <text x="49" y="32" font-family="sans-serif" font-size="10" fill="#FF0000">OH</text>
        <text x="69" y="32" font-family="sans-serif" font-size="10" fill="#FF0000">OH</text>
        <text x="69" y="62" font-family="sans-serif" font-size="10" fill="#FF0000">OH</text>
      </g>
    </svg>`,
    'CS(=O)C': `<svg version="1.1" id="mol-dmso" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" x="0px" y="0px" width="300px" height="200px" viewBox="0 0 300 200" enable-background="new 0 0 300 200">
      <rect fill="#FFFFFF" width="300" height="200"/>
      <g transform="translate(10,10) scale(2.8)">
        <path fill="none" stroke="#000000" stroke-width="1.5" d="M 50,50 L 60,50"/>
        <path fill="none" stroke="#000000" stroke-width="1.5" d="M 60,50 L 70,50"/>
        <path fill="none" stroke="#000000" stroke-width="1.5" d="M 60,50 L 60,40"/>
        <text x="56" y="54" font-family="sans-serif" font-size="10" fill="#AA5500">S</text>
        <text x="59" y="36" font-family="sans-serif" font-size="10" fill="#FF0000">O</text>
      </g>
    </svg>`
  },
  // Precomputed molecular properties for fallback
  molecular_properties: {
    'C(C(CO)O)O': {
      molecular_weight: 92.09,
      logp: -1.76,
      num_atoms: 14,
      num_rings: 0,
      formula: 'C3H8O3'
    },
    'CS(=O)C': {
      molecular_weight: 78.13,
      logp: -0.62,
      num_atoms: 10,
      num_rings: 0,
      formula: 'C2H6OS'
    },
    'OCCO': {
      molecular_weight: 62.07,
      logp: -1.2,
      num_atoms: 10,
      num_rings: 0,
      formula: 'C2H6O2'
    },
    'CC(O)CO': {
      molecular_weight: 76.09,
      logp: -0.92,
      num_atoms: 13,
      num_rings: 0,
      formula: 'C3H8O2'
    }
  }
};

// Utility to check if we're running in the browser
const isBrowser = typeof window !== 'undefined';

// Circuit breaker state - persistent across requests in the same session
const circuitBreaker = {
  HEROKU_API: {
    failureCount: 0,
    lastFailureTime: null,
    state: 'CLOSED', // CLOSED, OPEN, HALF_OPEN
    lastSuccessTime: Date.now()
  },
  RDKIT_API: {
    failureCount: 0,
    lastFailureTime: null,
    state: 'CLOSED', // CLOSED, OPEN, HALF_OPEN
    lastSuccessTime: Date.now()
  }
};

// Circuit breaker settings
const CIRCUIT_BREAKER_THRESHOLD = 5; // Number of failures before opening circuit
const CIRCUIT_BREAKER_TIMEOUT = 30000; // Time to wait before trying again (half-open)
const MAX_RETRIES = 3; // Maximum number of retries for a request
const RETRY_BACKOFF_MS = 1000; // Base backoff time in milliseconds
const REQUEST_TIMEOUT_MS = 10000; // Default request timeout (10 seconds)

/**
 * Check if circuit breaker is open for the given service
 * @param {string} service - Service to check (HEROKU_API, RDKIT_API)
 * @returns {boolean} - True if circuit is open
 */
function isCircuitOpen(service) {
  const breaker = circuitBreaker[service];
  
  if (!breaker) return false;
  
  if (breaker.state === 'OPEN') {
    // Check if enough time has passed to try again
    const now = Date.now();
    if (now - breaker.lastFailureTime > CIRCUIT_BREAKER_TIMEOUT) {
      // Move to half-open state
      breaker.state = 'HALF_OPEN';
      console.info(`Circuit breaker for ${service} moving to HALF_OPEN state`);
      return false;
    }
    return true;
  }
  
  return false;
}

/**
 * Record a success for the given service
 * @param {string} service - Service to record success for
 */
function recordSuccess(service) {
  const breaker = circuitBreaker[service];
  if (!breaker) return;
  
  // Update success metrics
  breaker.lastSuccessTime = Date.now();
  
  if (breaker.state === 'HALF_OPEN') {
    // Reset circuit breaker on successful test request
    breaker.failureCount = 0;
    breaker.state = 'CLOSED';
    console.info(`Circuit breaker for ${service} reset to CLOSED state`);
  } else if (breaker.state === 'CLOSED') {
    // In closed state, reset failure count after a successful request
    // This helps recover from intermittent failures
    breaker.failureCount = 0;
  }
}

/**
 * Record a failure for the given service
 * @param {string} service - Service to record failure for
 */
function recordFailure(service) {
  const breaker = circuitBreaker[service];
  if (!breaker) return;
  
  // Increment failure count and update timestamp
  breaker.failureCount++;
  breaker.lastFailureTime = Date.now();
  
  if (breaker.state === 'CLOSED' && breaker.failureCount >= CIRCUIT_BREAKER_THRESHOLD) {
    // Trip the circuit breaker when failure threshold is reached
    breaker.state = 'OPEN';
    console.warn(`Circuit breaker for ${service} tripped to OPEN state after ${CIRCUIT_BREAKER_THRESHOLD} failures`);
  } else if (breaker.state === 'HALF_OPEN') {
    // If test request fails in half-open state, return to open state
    breaker.state = 'OPEN';
    console.warn(`Circuit breaker for ${service} returned to OPEN state after failed test request`);
  }
}

/**
 * Enhanced fetch wrapper with error handling, retries, circuit breaker, and fallback to mock data
 * 
 * Features:
 * - Circuit breaker pattern to prevent cascading failures
 * - Exponential backoff retry mechanism
 * - Request timeout handling
 * - Graceful degradation with mock data fallbacks
 * - Detailed logging for debugging
 * 
 * @param {string} endpoint - API endpoint to fetch from (without base URL)
 * @param {Object} options - Fetch options
 * @param {string} mockDataKey - Key to access mock data if API fails
 * @param {number} retryCount - Current retry count (used internally)
 * @returns {Promise<any>} - Parsed response data
 */
async function fetchWithErrorHandling(endpoint, options = {}, mockDataKey = null, retryCount = 0) {
  const service = 'HEROKU_API';
  const fullUrl = `${API_ENDPOINTS.HEROKU_API}${endpoint}`;
  
  // Create a request ID for logging
  const requestId = Math.random().toString(36).substring(2, 10);
  
  // Check circuit breaker
  if (isCircuitOpen(service)) {
    console.warn(`[${requestId}] Circuit open for ${service}, falling back to mock data`);
    return getMockData(endpoint, mockDataKey);
  }
  
  // Create AbortController for timeout
  const controller = typeof AbortController !== 'undefined' ? new AbortController() : null;
  const timeoutDuration = options.timeout || REQUEST_TIMEOUT_MS;
  const timeoutId = controller ? setTimeout(() => controller.abort(), timeoutDuration) : null;
  
  try {
    console.info(`[${requestId}] Fetching from ${fullUrl} (attempt ${retryCount + 1}/${MAX_RETRIES + 1})`);
    
    const response = await fetch(fullUrl, {
      ...options,
      headers: {
        'Content-Type': 'application/json',
        'Accept': 'application/json',
        'X-Request-ID': requestId,
        ...options.headers,
      },
      mode: 'cors',
      credentials: 'omit',
      signal: options.signal || controller?.signal
    });

    // Clear timeout if set
    if (timeoutId) clearTimeout(timeoutId);

    if (!response.ok) {
      throw new Error(`API error: ${response.status} ${response.statusText}`);
    }

    const data = await response.json();
    
    // Record success for circuit breaker
    recordSuccess(service);
    
    return data;
  } catch (error) {
    // Clear timeout if still active
    if (timeoutId) clearTimeout(timeoutId);
    
    // Handle timeout errors specifically
    const isTimeout = error.name === 'AbortError';
    const errorMsg = isTimeout 
      ? `Request timed out after ${timeoutDuration}ms` 
      : error.message;
    
    console.error(`[${requestId}] API request failed (attempt ${retryCount + 1}/${MAX_RETRIES + 1}): ${errorMsg}`);
    
    // Record failure for circuit breaker
    recordFailure(service);
    
    // Check if we should retry
    if (retryCount < MAX_RETRIES) {
      const backoffTime = RETRY_BACKOFF_MS * Math.pow(2, retryCount);
      console.info(`[${requestId}] Retrying in ${backoffTime}ms...`);
      
      return new Promise(resolve => {
        setTimeout(() => {
          resolve(fetchWithErrorHandling(endpoint, options, mockDataKey, retryCount + 1));
        }, backoffTime);
      });
    }
    
    // Fall back to mock data after all retries are exhausted
    console.warn(`[${requestId}] All retries failed, falling back to mock data`);
    return getMockData(endpoint, mockDataKey);
  }
}

/**
 * Get mock data for an endpoint
 * @param {string} endpoint - API endpoint
 * @param {string} mockDataKey - Key to access mock data
 * @returns {any} - Mock data
 */
function getMockData(endpoint, mockDataKey) {
  // If we have mock data for this key, return it
  if (mockDataKey && MOCK_DATA[mockDataKey]) {
    console.info(`Using mock data for ${endpoint} with key ${mockDataKey}`);
    
    // Special case for single item requests
    if (endpoint.includes('/')) {
      const parts = endpoint.split('/');
      const id = parts[parts.length - 1];
      
      // If the mockDataKey points to an array, find the item with matching ID
      const items = MOCK_DATA[mockDataKey];
      
      if (items && Array.isArray(items)) {
        const item = items.find(i => String(i.id) === String(id));
        if (item) {
          return item;
        }
      }
    }
    
    return MOCK_DATA[mockDataKey];
  }
  
  // Default mock data for unknown endpoints
  if (endpoint.includes('molecules')) {
    return Array.isArray(MOCK_DATA.molecules) ? MOCK_DATA.molecules[0] : null;
  } else if (endpoint.includes('mixtures')) {
    return Array.isArray(MOCK_DATA.mixtures) ? MOCK_DATA.mixtures[0] : null;
  }
  
  // No mock data available
  throw new Error(`No mock data available for ${endpoint}`);
}

/**
 * Fetch wrapper for RDKit API with circuit breaker
 * @param {string} endpoint - API endpoint to fetch from (e.g., '/calculate', '/depict')
 * @param {Object} data - Request data to send
 * @param {Object} options - Additional fetch options
 * @param {number} retryCount - Current retry count (used internally)
 * @returns {Promise<any>} - Parsed response data
 */
async function fetchRDKitAPI(endpoint, data, options = {}, retryCount = 0) {
  const service = 'RDKIT_API';
  const fullUrl = `${API_ENDPOINTS.RDKIT_API}${endpoint}`;
  
  // Create a request ID for logging
  const requestId = Math.random().toString(36).substring(2, 10);
  
  // Check circuit breaker
  if (isCircuitOpen(service)) {
    console.warn(`[${requestId}] Circuit open for ${service}`);
    throw new Error(`RDKit service unavailable (circuit open)`);
  }
  
  // Create AbortController for timeout
  const controller = typeof AbortController !== 'undefined' ? new AbortController() : null;
  const timeoutDuration = options.timeout || REQUEST_TIMEOUT_MS;
  const timeoutId = controller ? setTimeout(() => controller.abort(), timeoutDuration) : null;
  
  try {
    console.info(`[${requestId}] Fetching from ${fullUrl} (attempt ${retryCount + 1}/${MAX_RETRIES + 1})`);
    
    const response = await fetch(fullUrl, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
        'Accept': 'application/json',
        'X-Request-ID': requestId,
        ...options.headers
      },
      body: JSON.stringify(data),
      mode: 'cors',
      credentials: 'omit',
      signal: controller?.signal
    });

    // Clear timeout
    if (timeoutId) clearTimeout(timeoutId);

    if (!response.ok) {
      throw new Error(`RDKit API error: ${response.status} ${response.statusText}`);
    }

    const responseData = await response.json();
    
    // Record success
    recordSuccess(service);
    
    return responseData;
  } catch (error) {
    // Clear timeout if still active
    if (timeoutId) clearTimeout(timeoutId);
    
    // Handle timeout errors specifically
    const isTimeout = error.name === 'AbortError';
    const errorMsg = isTimeout 
      ? `Request timed out after ${timeoutDuration}ms` 
      : error.message;
    
    console.error(`[${requestId}] RDKit API request failed (attempt ${retryCount + 1}/${MAX_RETRIES + 1}): ${errorMsg}`);
    
    // Record failure
    recordFailure(service);
    
    // Check if we should retry
    if (retryCount < MAX_RETRIES) {
      const backoffTime = RETRY_BACKOFF_MS * Math.pow(2, retryCount);
      console.info(`[${requestId}] Retrying RDKit API in ${backoffTime}ms...`);
      
      return new Promise(resolve => {
        setTimeout(() => {
          resolve(fetchRDKitAPI(endpoint, data, options, retryCount + 1));
        }, backoffTime);
      });
    }
    
    throw error;
  }
}

/**
 * Get all molecules
 * @returns {Promise<Array>} List of molecules
 */
export async function getMolecules() {
  // Convex functionality disabled for minimal deployment
  // We'll always use the REST API or mock data
  
  // Fall back to REST API
  return fetchWithErrorHandling('/molecules', {}, 'molecules');
}

/**
 * Get a single molecule by ID
 * @param {string|number} id - Molecule ID
 * @returns {Promise<Object>} Molecule data
 */
export async function getMolecule(id) {
  // Convex functionality disabled for minimal deployment
  // We'll always use the REST API or mock data
  
  // Fall back to REST API
  return fetchWithErrorHandling(`/molecules/${id}`, {}, 'molecules');
}

/**
 * Get all mixtures 
 * @returns {Promise<Array>} List of mixtures
 */
export async function getMixtures() {
  // Convex functionality disabled for minimal deployment
  // We'll always use the REST API or mock data
  
  // Fall back to REST API
  return fetchWithErrorHandling('/mixtures', {}, 'mixtures');
}

/**
 * Get a single mixture by ID
 * @param {string|number} id - Mixture ID
 * @returns {Promise<Object>} Mixture data
 */
export async function getMixture(id) {
  // Convex functionality disabled for minimal deployment
  // We'll always use the REST API or mock data
  
  // Fall back to REST API
  return fetchWithErrorHandling(`/mixtures/${id}`, {}, 'mixtures');
}

/**
 * Get RDKit molecular properties with circuit breaker and fallback
 * @param {string} smiles - SMILES notation of the molecule
 * @returns {Promise<Object>} RDKit calculation results
 */
export async function getRDKitProperties(smiles) {
  if (!smiles) {
    return {
      error: 'No SMILES provided',
      molecular_weight: 0,
      logp: 0,
      num_atoms: 0,
      num_rings: 0
    };
  }
  
  try {
    // Request properties from RDKit service
    const data = await fetchRDKitAPI('/calculate', { smiles });
    return data;
  } catch (error) {
    console.error('RDKit properties calculation failed:', error);
    
    // Use fallback property calculation
    if (MOCK_DATA.molecular_properties[smiles]) {
      console.info(`Using pre-computed properties for ${smiles}`);
      return {
        ...MOCK_DATA.molecular_properties[smiles],
        source: 'fallback'
      };
    }
    
    // If no pre-computed data, use rough estimation
    return {
      molecular_weight: estimateMolecularWeight(smiles),
      logp: 0,
      num_atoms: countAtoms(smiles),
      num_rings: estimateRingCount(smiles),
      error: 'Failed to fetch properties from RDKit service, using estimated values',
      source: 'estimated'
    };
  }
}

/**
 * Generate a 2D depiction of a molecule with circuit breaker and fallback
 * @param {string} smiles - SMILES notation of the molecule
 * @param {Object} options - Optional rendering options (width, height, includeMoleculeDetails)
 * @returns {Promise<string>} SVG depiction
 */
export async function getMoleculeDepiction(smiles, options = {}) {
  if (!smiles) return '';
  
  const requestData = { 
    smiles, 
    width: options.width || 300, 
    height: options.height || 200,
    includeMoleculeDetails: options.includeMoleculeDetails || false
  };
  
  try {
    // Request depiction from RDKit service
    const data = await fetchRDKitAPI('/depict', requestData);
    return data.svg || '';
  } catch (error) {
    console.error('RDKit depiction generation failed:', error);
    
    // Try to use cached SVG depiction if available
    if (MOCK_DATA.molecular_depictions[smiles]) {
      console.info(`Using cached depiction for ${smiles}`);
      return MOCK_DATA.molecular_depictions[smiles];
    }
    
    // Otherwise return empty string to show fallback UI
    return '';
  }
}

/**
 * Very basic estimation of molecular weight from SMILES (fallback method)
 * This is a very rough approximation and should only be used when RDKit is unavailable
 * @param {string} smiles - SMILES notation
 * @returns {number} Estimated molecular weight
 */
function estimateMolecularWeight(smiles) {
  // This is an extremely simplified approach
  // A real implementation would need to parse the SMILES properly
  
  // Count some common atoms
  const carbonCount = (smiles.match(/C/g) || []).length;
  const hydrogenCount = (smiles.match(/H/g) || []).length;
  const oxygenCount = (smiles.match(/O/g) || []).length;
  const nitrogenCount = (smiles.match(/N/g) || []).length;
  const sulfurCount = (smiles.match(/S/g) || []).length;
  const phosphorusCount = (smiles.match(/P/g) || []).length;
  const fluorineCount = (smiles.match(/F/g) || []).length;
  const chlorineCount = (smiles.match(/Cl/g) || []).length;
  const bromineCount = (smiles.match(/Br/g) || []).length;
  const iodineCount = (smiles.match(/I/g) || []).length;
  
  // Approximate weights
  const C_WEIGHT = 12.01;
  const H_WEIGHT = 1.01;
  const O_WEIGHT = 16.00;
  const N_WEIGHT = 14.01;
  const S_WEIGHT = 32.07;
  const P_WEIGHT = 30.97;
  const F_WEIGHT = 19.00;
  const CL_WEIGHT = 35.45;
  const BR_WEIGHT = 79.90;
  const I_WEIGHT = 126.90;
  
  // Rough estimate of weight from atom counts
  const weight = (
    carbonCount * C_WEIGHT + 
    hydrogenCount * H_WEIGHT + 
    oxygenCount * O_WEIGHT + 
    nitrogenCount * N_WEIGHT +
    sulfurCount * S_WEIGHT +
    phosphorusCount * P_WEIGHT +
    fluorineCount * F_WEIGHT +
    chlorineCount * CL_WEIGHT +
    bromineCount * BR_WEIGHT +
    iodineCount * I_WEIGHT
  );
  
  // Add an estimate for implicit hydrogens
  // This is a very rough heuristic
  const implicitHydrogens = Math.max(0, carbonCount * 2);
  
  return Math.round((weight + implicitHydrogens * H_WEIGHT) * 100) / 100;
}

/**
 * Count atoms in SMILES (fallback method)
 * @param {string} smiles - SMILES notation
 * @returns {number} Approximate atom count
 */
function countAtoms(smiles) {
  // This is a very simplified approach
  // We'll count the uppercase letters as a rough proxy for atoms
  // We'll also check for some common two-letter atoms
  
  // Basic pattern for any atom symbol (uppercase letter possibly followed by lowercase letter)
  const atomPattern = /[A-Z][a-z]?/g;
  const matches = smiles.match(atomPattern) || [];
  
  return matches.length;
}

/**
 * Estimate the number of rings in a molecule from SMILES
 * This is a very basic heuristic for when RDKit is unavailable
 * @param {string} smiles - SMILES notation
 * @returns {number} Approximate ring count
 */
function estimateRingCount(smiles) {
  // Look for digit pairs which indicate ring closures in SMILES
  const ringDigits = smiles.match(/\d/g) || [];
  
  // Each ring should be counted once (two digits per ring)
  return Math.floor(ringDigits.length / 2);
}

/**
 * Check Convex connection status
 * @returns {Promise<boolean>} True if Convex is connected
 */
export async function checkConvexConnection() {
  // Convex functionality disabled for minimal deployment
  return false;
}

/**
 * Add a new molecule using Convex
 * @param {Object} moleculeData - Data for the new molecule
 * @returns {Promise<Object>} Result of the operation
 */
export async function addMolecule(moleculeData) {
  // Convex functionality disabled for minimal deployment
  throw new Error('Adding molecules is not available in this minimal deployment');
}

/**
 * Delete a molecule using Convex
 * @param {string} id - ID of the molecule to delete
 * @returns {Promise<Object>} Result of the operation
 */
export async function deleteMolecule(id) {
  // Convex functionality disabled for minimal deployment
  throw new Error('Deleting molecules is not available in this minimal deployment');
}

export default {
  getMolecules,
  getMolecule,
  getMixtures,
  getMixture,
  getRDKitProperties,
  getMoleculeDepiction,
  checkConvexConnection,
  addMolecule,
  deleteMolecule
};
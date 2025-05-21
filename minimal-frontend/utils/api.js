/**
 * API utility functions for fetching data from the CryoProtect backend
 * Supports both REST API and Convex
 */

import { convex } from '../src/convex/client';

// Base API URLs
const API_ENDPOINTS = {
  HEROKU_API: 'https://cryoprotect-8030e4025428.herokuapp.com/api',
  RDKIT_API: 'https://cryoprotect-rdkit.fly.dev',
  CONVEX: 'https://primary-meerkat-478.convex.cloud'
};

// Check whether to use Convex
const USE_CONVEX = process.env.NEXT_PUBLIC_USE_CONVEX === 'true';

// Fallback mock data in case the API is not available
const MOCK_DATA = {
  molecules: [
    {
      id: 1,
      name: 'Glycerol',
      formula: 'C3H8O3',
      pubchem_cid: '753',
      molecular_weight: 92.09,
      is_cryoprotectant: true,
      description: 'A common cryoprotectant used in various applications.'
    },
    {
      id: 2,
      name: 'Dimethyl Sulfoxide (DMSO)',
      formula: 'C2H6OS',
      pubchem_cid: '679',
      molecular_weight: 78.13,
      is_cryoprotectant: true,
      description: 'A widely used penetrating cryoprotectant.'
    },
    {
      id: 3,
      name: 'Ethylene Glycol',
      formula: 'C2H6O2',
      pubchem_cid: '174',
      molecular_weight: 62.07,
      is_cryoprotectant: true,
      description: 'Used in cryopreservation of embryos and tissues.'
    },
    {
      id: 4,
      name: 'Propylene Glycol',
      formula: 'C3H8O2',
      pubchem_cid: '1030',
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
  ]
};

// Utility to check if we're running in the browser
const isBrowser = typeof window !== 'undefined';

/**
 * Fetch wrapper with error handling and fallback to mock data
 * @param {string} endpoint - API endpoint to fetch from (without base URL)
 * @param {Object} options - Fetch options
 * @param {string} mockDataKey - Key to access mock data if API fails
 * @returns {Promise<any>} - Parsed response data
 */
async function fetchWithErrorHandling(endpoint, options = {}, mockDataKey = null) {
  try {
    // Apply circuit breaker here if needed

    const response = await fetch(`${API_ENDPOINTS.HEROKU_API}${endpoint}`, {
      ...options,
      headers: {
        'Content-Type': 'application/json',
        'Accept': 'application/json',
        ...options.headers,
      },
      mode: 'cors',
      credentials: 'omit'
    });

    if (!response.ok) {
      throw new Error(`API error: ${response.status} ${response.statusText}`);
    }

    return await response.json();
  } catch (error) {
    console.error('API request failed:', error);
    
    // If we have mock data for this endpoint, return it
    if (mockDataKey && MOCK_DATA[mockDataKey]) {
      console.info(`Falling back to mock data for ${endpoint}`);
      return MOCK_DATA[mockDataKey];
    }
    
    // For single item requests, try to find it in mock data
    if (endpoint.includes('/') && mockDataKey) {
      const id = endpoint.split('/').pop();
      const items = MOCK_DATA[mockDataKey];
      
      if (items && Array.isArray(items)) {
        const item = items.find(i => i.id.toString() === id.toString());
        if (item) {
          console.info(`Falling back to mock data for ${endpoint}`);
          return item;
        }
      }
    }
    
    throw error;
  }
}

/**
 * Get all molecules
 * @returns {Promise<Array>} List of molecules
 */
export async function getMolecules() {
  // If we're using Convex and in the browser, try to use it
  if (USE_CONVEX && isBrowser && convex) {
    try {
      // Use our Convex client directly to fetch molecules
      const { api } = await import('../convex/_generated/api');
      const molecules = await convex.query(api.molecules.query.getRecentMolecules, { limit: 20 });
      
      if (molecules && Array.isArray(molecules)) {
        // Adapt Convex data to match the API format
        return molecules.map(molecule => ({
          id: molecule._id,
          name: molecule.name,
          pubchem_cid: molecule.pubchemCid,
          formula: molecule.formula,
          molecular_weight: molecule.molecularWeight,
          is_cryoprotectant: molecule.isCryoprotectant,
          description: molecule.description
        }));
      }
    } catch (error) {
      console.error('Convex query failed:', error);
    }
  }
  
  // Fall back to REST API
  return fetchWithErrorHandling('/molecules', {}, 'molecules');
}

/**
 * Get a single molecule by ID
 * @param {string|number} id - Molecule ID
 * @returns {Promise<Object>} Molecule data
 */
export async function getMolecule(id) {
  // If we're using Convex and in the browser, try to use it
  if (USE_CONVEX && isBrowser && convex) {
    try {
      // Use our Convex client directly to fetch a single molecule
      const { api } = await import('../convex/_generated/api');
      const molecule = await convex.query(api.molecules.query.getMolecule, { id });
      
      if (molecule) {
        // Adapt Convex data to match the API format
        return {
          id: molecule._id,
          name: molecule.name,
          pubchem_cid: molecule.pubchemCid,
          formula: molecule.formula,
          molecular_weight: molecule.molecularWeight,
          is_cryoprotectant: molecule.isCryoprotectant,
          description: molecule.description
        };
      }
    } catch (error) {
      console.error('Convex query failed:', error);
    }
  }
  
  // Fall back to REST API
  return fetchWithErrorHandling(`/molecules/${id}`, {}, 'molecules');
}

/**
 * Get all mixtures 
 * @returns {Promise<Array>} List of mixtures
 */
export async function getMixtures() {
  // If we're using Convex and in the browser, try to use it
  if (USE_CONVEX && isBrowser && convex) {
    try {
      // Use our Convex client directly to fetch mixtures
      const { api } = await import('../convex/_generated/api');
      const mixtures = await convex.query(api.mixtures.mixtures.getRecentMixtures, { limit: 20 });
      
      if (mixtures && Array.isArray(mixtures)) {
        // Adapt Convex data to match the API format
        return mixtures.map(mixture => ({
          id: mixture._id,
          name: mixture.name,
          description: mixture.description,
          freezing_point: mixture.freezingPoint,
          components: mixture.components || []
        }));
      }
    } catch (error) {
      console.error('Convex query failed:', error);
    }
  }
  
  // Fall back to REST API
  return fetchWithErrorHandling('/mixtures', {}, 'mixtures');
}

/**
 * Get a single mixture by ID
 * @param {string|number} id - Mixture ID
 * @returns {Promise<Object>} Mixture data
 */
export async function getMixture(id) {
  // If we're using Convex and in the browser, try to use it
  if (USE_CONVEX && isBrowser && convex) {
    try {
      // Use our Convex client directly to fetch a single mixture
      const { api } = await import('../convex/_generated/api');
      const mixture = await convex.query(api.mixtures.mixtures.getMixture, { id });
      
      if (mixture) {
        // Adapt Convex data to match the API format
        return {
          id: mixture._id,
          name: mixture.name,
          description: mixture.description,
          freezing_point: mixture.freezingPoint,
          components: mixture.components || []
        };
      }
    } catch (error) {
      console.error('Convex query failed:', error);
    }
  }
  
  // Fall back to REST API
  return fetchWithErrorHandling(`/mixtures/${id}`, {}, 'mixtures');
}

/**
 * Get RDKit molecular properties
 * @param {string} smiles - SMILES notation of the molecule
 * @returns {Promise<Object>} RDKit calculation results
 */
export async function getRDKitProperties(smiles) {
  try {
    const response = await fetch(`${API_ENDPOINTS.RDKIT_API}/calculate`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
        'Accept': 'application/json',
      },
      body: JSON.stringify({ smiles }),
      mode: 'cors',
      credentials: 'omit'
    });

    if (!response.ok) {
      throw new Error(`RDKit API error: ${response.status} ${response.statusText}`);
    }

    return await response.json();
  } catch (error) {
    console.error('RDKit API request failed:', error);
    
    // Return minimal mock data for RDKit
    return {
      molecular_weight: 0,
      logp: 0,
      num_atoms: 0,
      num_rings: 0,
      error: 'Failed to fetch properties from RDKit service'
    };
  }
}

/**
 * Generate a 2D depiction of a molecule
 * @param {string} smiles - SMILES notation of the molecule
 * @returns {Promise<string>} SVG depiction
 */
export async function getMoleculeDepiction(smiles) {
  try {
    const response = await fetch(`${API_ENDPOINTS.RDKIT_API}/depict`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
        'Accept': 'application/json',
      },
      body: JSON.stringify({ smiles, width: 300, height: 200 }),
      mode: 'cors',
      credentials: 'omit'
    });

    if (!response.ok) {
      throw new Error(`RDKit API error: ${response.status} ${response.statusText}`);
    }

    const data = await response.json();
    return data.svg || '';
  } catch (error) {
    console.error('RDKit API request failed:', error);
    return '';
  }
}

/**
 * Check Convex connection status
 * @returns {Promise<boolean>} True if Convex is connected
 */
export async function checkConvexConnection() {
  if (!USE_CONVEX || !isBrowser || !convex) {
    return false;
  }
  
  try {
    return await convex.health();
  } catch (error) {
    console.error('Convex health check failed:', error);
    return false;
  }
}

/**
 * Add a new molecule using Convex
 * @param {Object} moleculeData - Data for the new molecule
 * @returns {Promise<Object>} Result of the operation
 */
export async function addMolecule(moleculeData) {
  if (!USE_CONVEX || !isBrowser || !convex) {
    throw new Error('Convex is not available');
  }
  
  try {
    const { api } = await import('../convex/_generated/api');
    const result = await convex.mutation(api.molecules.mutations.addMolecule, moleculeData);
    return result;
  } catch (error) {
    console.error('Failed to add molecule:', error);
    throw error;
  }
}

/**
 * Delete a molecule using Convex
 * @param {string} id - ID of the molecule to delete
 * @returns {Promise<Object>} Result of the operation
 */
export async function deleteMolecule(id) {
  if (!USE_CONVEX || !isBrowser || !convex) {
    throw new Error('Convex is not available');
  }
  
  try {
    const { api } = await import('../convex/_generated/api');
    const result = await convex.mutation(api.molecules.mutations.deleteMolecule, { id });
    return result;
  } catch (error) {
    console.error('Failed to delete molecule:', error);
    throw error;
  }
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
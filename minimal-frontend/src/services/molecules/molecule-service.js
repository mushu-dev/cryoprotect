/**
 * Molecule Service for CryoProtect
 * 
 * Provides methods for accessing molecule data from both the API and Convex.
 * Automatically selects the appropriate data source based on the NEXT_PUBLIC_USE_CONVEX
 * environment variable.
 */
import axios from 'axios';
import { isEnabled } from '../../../shared/convex/client';

// Base URL for API requests
const BASE_URL = process.env.NEXT_PUBLIC_API_URL || '/api';

/**
 * Molecule Service with both API and Convex support
 */
export const moleculeService = {
  /**
   * Get a list of molecules with pagination
   * 
   * @param {Object} params - Query parameters
   * @param {number} params.page - Page number for pagination
   * @param {number} params.limit - Number of items per page
   * @param {string} params.search - Search query
   * @param {string} params.sort_by - Field to sort by
   * @param {string} params.sort_order - Sort order (asc or desc)
   * @param {boolean} params.is_cryoprotectant - Filter by cryoprotectant status
   * @returns {Promise<Object>} - Molecules response with pagination info
   */
  async getMolecules(params = {}) {
    try {
      // Use Convex if enabled
      if (isEnabled() && typeof window !== 'undefined' && window.convexMolecules) {
        // Call Convex function to get molecules
        const result = await window.convexMolecules.getAll({
          limit: params.limit || 10,
          search: params.search || '',
          sortBy: params.sort_by,
          sortOrder: params.sort_order,
          isCryoprotectant: params.is_cryoprotectant
        });
        
        return {
          molecules: result.molecules || [],
          total: result.total || 0,
          page: params.page || 1,
          limit: params.limit || 10
        };
      }
      
      // Fall back to API
      const response = await axios.get(`${BASE_URL}/v1/molecules`, { 
        params,
        headers: {
          'Content-Type': 'application/json',
        }
      });
      return response.data;
    } catch (error) {
      console.error('Error fetching molecules:', error);
      throw error;
    }
  },

  /**
   * Get a single molecule by ID
   * 
   * @param {string} id - Molecule ID
   * @returns {Promise<Object>} - Molecule data
   */
  async getMolecule(id) {
    try {
      // Use Convex if enabled
      if (isEnabled() && typeof window !== 'undefined' && window.convexMolecules) {
        const result = await window.convexMolecules.getById(id);
        return result || null;
      }
      
      // Fall back to API
      const response = await axios.get(`${BASE_URL}/v1/molecules/${id}`, {
        headers: {
          'Content-Type': 'application/json',
        }
      });
      return response.data;
    } catch (error) {
      console.error(`Error fetching molecule ${id}:`, error);
      throw error;
    }
  },

  /**
   * Get properties for a molecule
   * 
   * @param {string} id - Molecule ID
   * @returns {Promise<Array>} - Array of molecule properties
   */
  async getMoleculeProperties(id) {
    try {
      // Use Convex if enabled
      if (isEnabled() && typeof window !== 'undefined' && window.convexMolecules) {
        const result = await window.convexMolecules.getProperties(id);
        return result.properties || [];
      }
      
      // Fall back to API
      const response = await axios.get(`${BASE_URL}/v1/molecules/${id}/properties`, {
        headers: {
          'Content-Type': 'application/json',
        }
      });
      return response.data.properties || [];
    } catch (error) {
      console.error(`Error fetching properties for molecule ${id}:`, error);
      throw error;
    }
  },

  /**
   * Search molecules by name or SMILES
   * 
   * @param {string} query - Search query
   * @param {number} limit - Maximum number of results
   * @returns {Promise<Array>} - Array of molecules matching the search
   */
  async searchMolecules(query, limit = 10) {
    try {
      // Use Convex if enabled
      if (isEnabled() && typeof window !== 'undefined' && window.convexMolecules) {
        const result = await window.convexMolecules.search(query, limit);
        return result.results || [];
      }
      
      // Fall back to API
      const response = await axios.get(`${BASE_URL}/v1/molecules/search`, {
        params: { query, limit },
        headers: {
          'Content-Type': 'application/json',
        }
      });
      return response.data.results || [];
    } catch (error) {
      console.error('Error searching molecules:', error);
      throw error;
    }
  },

  /**
   * Import a molecule from PubChem by CID
   * 
   * @param {string} cid - PubChem CID
   * @returns {Promise<Object>} - Imported molecule data
   */
  async importFromPubChem(cid) {
    try {
      // Use Convex if enabled
      if (isEnabled() && typeof window !== 'undefined' && window.convexMolecules) {
        const result = await window.convexMolecules.importFromPubChem(cid);
        return result || null;
      }
      
      // Fall back to API
      const response = await axios.post(
        `${BASE_URL}/v1/molecules/import/pubchem`,
        { cid },
        {
          headers: {
            'Content-Type': 'application/json',
          }
        }
      );
      return response.data;
    } catch (error) {
      console.error(`Error importing molecule from PubChem (CID: ${cid}):`, error);
      throw error;
    }
  },
  
  /**
   * Initialize Convex functions for molecules
   * 
   * @param {Object} convexFunctions - Object containing Convex query functions
   */
  initConvex(convexFunctions) {
    if (typeof window !== 'undefined') {
      window.convexMolecules = convexFunctions;
    }
  }
};
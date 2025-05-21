/**
 * Mixture Service for CryoProtect
 * 
 * Provides methods for accessing mixture data from both the API and Convex.
 * Automatically selects the appropriate data source based on the NEXT_PUBLIC_USE_CONVEX
 * environment variable.
 */
import axios from 'axios';
import { isEnabled } from '../../../shared/convex/client';

// Base URL for API requests
const BASE_URL = process.env.NEXT_PUBLIC_API_URL || '/api';

/**
 * Mixture Service with both API and Convex support
 */
export const mixtureService = {
  /**
   * Get a list of mixtures with pagination
   * 
   * @param {Object} params - Query parameters
   * @param {number} params.page - Page number for pagination
   * @param {number} params.limit - Number of items per page
   * @param {string} params.search - Search query
   * @param {string} params.sort_by - Field to sort by
   * @param {string} params.sort_order - Sort order (asc or desc)
   * @param {boolean} params.is_cryoprotectant_mixture - Filter by cryoprotectant status
   * @param {string} params.user_id - Filter by user ID
   * @returns {Promise<Object>} - Mixtures response with pagination info
   */
  async getMixtures(params = {}) {
    try {
      // Use Convex if enabled
      if (isEnabled() && typeof window !== 'undefined' && window.convexMixtures) {
        // Call Convex function to get mixtures
        const result = await window.convexMixtures.getAll({
          limit: params.limit || 10,
          search: params.search || '',
          sortBy: params.sort_by,
          sortOrder: params.sort_order,
          isCryoprotectant: params.is_cryoprotectant_mixture,
          userId: params.user_id
        });
        
        return {
          mixtures: result.mixtures || [],
          total: result.total || 0,
          page: params.page || 1,
          limit: params.limit || 10
        };
      }
      
      // Fall back to API
      const response = await axios.get(`${BASE_URL}/v1/mixtures`, { 
        params,
        headers: {
          'Content-Type': 'application/json',
        }
      });
      return response.data;
    } catch (error) {
      console.error('Error fetching mixtures:', error);
      throw error;
    }
  },

  /**
   * Get a single mixture by ID
   * 
   * @param {string} id - Mixture ID
   * @returns {Promise<Object>} - Mixture data
   */
  async getMixture(id) {
    try {
      // Use Convex if enabled
      if (isEnabled() && typeof window !== 'undefined' && window.convexMixtures) {
        const result = await window.convexMixtures.getById(id);
        return result || null;
      }
      
      // Fall back to API
      const response = await axios.get(`${BASE_URL}/v1/mixtures/${id}`, {
        headers: {
          'Content-Type': 'application/json',
        }
      });
      return response.data;
    } catch (error) {
      console.error(`Error fetching mixture ${id}:`, error);
      throw error;
    }
  },

  /**
   * Create a new mixture
   * 
   * @param {Object} data - Mixture data
   * @param {string} data.name - Mixture name
   * @param {string} data.description - Mixture description
   * @param {Array} data.components - Array of mixture components
   * @returns {Promise<Object>} - Created mixture data
   */
  async createMixture(data) {
    try {
      // Use Convex if enabled
      if (isEnabled() && typeof window !== 'undefined' && window.convexMixtures) {
        const result = await window.convexMixtures.create(data);
        return result || null;
      }
      
      // Fall back to API
      const response = await axios.post(
        `${BASE_URL}/v1/mixtures`,
        data,
        {
          headers: {
            'Content-Type': 'application/json',
          }
        }
      );
      return response.data;
    } catch (error) {
      console.error('Error creating mixture:', error);
      throw error;
    }
  },

  /**
   * Update a mixture
   * 
   * @param {string} id - Mixture ID
   * @param {Object} data - Data to update
   * @returns {Promise<Object>} - Updated mixture data
   */
  async updateMixture(id, data) {
    try {
      // Use Convex if enabled
      if (isEnabled() && typeof window !== 'undefined' && window.convexMixtures) {
        const result = await window.convexMixtures.update(id, data);
        return result || null;
      }
      
      // Fall back to API
      const response = await axios.patch(
        `${BASE_URL}/v1/mixtures/${id}`,
        data,
        {
          headers: {
            'Content-Type': 'application/json',
          }
        }
      );
      return response.data;
    } catch (error) {
      console.error(`Error updating mixture ${id}:`, error);
      throw error;
    }
  },

  /**
   * Delete a mixture
   * 
   * @param {string} id - Mixture ID
   * @returns {Promise<void>}
   */
  async deleteMixture(id) {
    try {
      // Use Convex if enabled
      if (isEnabled() && typeof window !== 'undefined' && window.convexMixtures) {
        await window.convexMixtures.delete(id);
        return;
      }
      
      // Fall back to API
      await axios.delete(`${BASE_URL}/v1/mixtures/${id}`, {
        headers: {
          'Content-Type': 'application/json',
        }
      });
    } catch (error) {
      console.error(`Error deleting mixture ${id}:`, error);
      throw error;
    }
  },

  /**
   * Add a component to a mixture
   * 
   * @param {string} mixtureId - Mixture ID
   * @param {Object} data - Component data
   * @returns {Promise<Object>} - Created component data
   */
  async addComponent(mixtureId, data) {
    try {
      // Use Convex if enabled
      if (isEnabled() && typeof window !== 'undefined' && window.convexMixtures) {
        const result = await window.convexMixtures.addComponent(mixtureId, data);
        return result || null;
      }
      
      // Fall back to API
      const response = await axios.post(
        `${BASE_URL}/v1/mixtures/${mixtureId}/components`,
        data,
        {
          headers: {
            'Content-Type': 'application/json',
          }
        }
      );
      return response.data;
    } catch (error) {
      console.error(`Error adding component to mixture ${mixtureId}:`, error);
      throw error;
    }
  },

  /**
   * Update a mixture component
   * 
   * @param {string} mixtureId - Mixture ID
   * @param {string} componentId - Component ID
   * @param {Object} data - Data to update
   * @returns {Promise<Object>} - Updated component data
   */
  async updateComponent(mixtureId, componentId, data) {
    try {
      // Use Convex if enabled
      if (isEnabled() && typeof window !== 'undefined' && window.convexMixtures) {
        const result = await window.convexMixtures.updateComponent(mixtureId, componentId, data);
        return result || null;
      }
      
      // Fall back to API
      const response = await axios.patch(
        `${BASE_URL}/v1/mixtures/${mixtureId}/components/${componentId}`,
        data,
        {
          headers: {
            'Content-Type': 'application/json',
          }
        }
      );
      return response.data;
    } catch (error) {
      console.error(`Error updating component ${componentId} in mixture ${mixtureId}:`, error);
      throw error;
    }
  },

  /**
   * Remove a component from a mixture
   * 
   * @param {string} mixtureId - Mixture ID
   * @param {string} componentId - Component ID
   * @returns {Promise<void>}
   */
  async removeComponent(mixtureId, componentId) {
    try {
      // Use Convex if enabled
      if (isEnabled() && typeof window !== 'undefined' && window.convexMixtures) {
        await window.convexMixtures.removeComponent(mixtureId, componentId);
        return;
      }
      
      // Fall back to API
      await axios.delete(
        `${BASE_URL}/v1/mixtures/${mixtureId}/components/${componentId}`,
        {
          headers: {
            'Content-Type': 'application/json',
          }
        }
      );
    } catch (error) {
      console.error(`Error removing component ${componentId} from mixture ${mixtureId}:`, error);
      throw error;
    }
  },

  /**
   * Get cryoprotection score for a mixture
   * 
   * @param {string} mixtureId - Mixture ID
   * @returns {Promise<Object>} - Score data
   */
  async getCryoprotectionScore(mixtureId) {
    try {
      // Use Convex if enabled
      if (isEnabled() && typeof window !== 'undefined' && window.convexMixtures) {
        const result = await window.convexMixtures.getCryoprotectionScore(mixtureId);
        return result || { score: 0 };
      }
      
      // Fall back to API
      const response = await axios.get(
        `${BASE_URL}/v1/mixtures/${mixtureId}/cryoprotection-score`,
        {
          headers: {
            'Content-Type': 'application/json',
          }
        }
      );
      return response.data;
    } catch (error) {
      console.error(`Error fetching cryoprotection score for mixture ${mixtureId}:`, error);
      throw error;
    }
  },

  /**
   * Search mixtures by name
   * 
   * @param {string} query - Search query
   * @param {number} limit - Maximum number of results
   * @returns {Promise<Array>} - Array of mixtures matching the search
   */
  async searchMixtures(query, limit = 10) {
    try {
      // Use Convex if enabled
      if (isEnabled() && typeof window !== 'undefined' && window.convexMixtures) {
        const result = await window.convexMixtures.search(query, limit);
        return result.results || [];
      }
      
      // Fall back to API
      const response = await axios.get(`${BASE_URL}/v1/mixtures/search`, {
        params: { query, limit },
        headers: {
          'Content-Type': 'application/json',
        }
      });
      return response.data.results || [];
    } catch (error) {
      console.error('Error searching mixtures:', error);
      throw error;
    }
  },
  
  /**
   * Initialize Convex functions for mixtures
   * 
   * @param {Object} convexFunctions - Object containing Convex query functions
   */
  initConvex(convexFunctions) {
    if (typeof window !== 'undefined') {
      window.convexMixtures = convexFunctions;
    }
  }
};
import axios from 'axios';
import { Molecule } from '@/features/molecules/services/molecule-service';

export interface MixtureComponent {
  id: string;
  mixture_id: string;
  molecule_id: string;
  molecule?: Molecule;
  concentration: number;
  concentration_unit: string;
  created_at?: string;
  updated_at?: string;
}

export interface Mixture {
  id: string;
  name: string;
  description?: string;
  is_cryoprotectant_mixture?: boolean;
  user_id?: string;
  created_at?: string;
  updated_at?: string;
  components?: MixtureComponent[];
  total_concentration?: number;
  cryoprotection_score?: number;
}

export interface MixturesResponse {
  mixtures: Mixture[];
  total: number;
  page: number;
  limit: number;
}

export interface MixtureParams {
  page?: number;
  limit?: number;
  search?: string;
  sort_by?: string;
  sort_order?: 'asc' | 'desc';
  is_cryoprotectant_mixture?: boolean;
  user_id?: string;
}

export interface CreateMixtureComponentData {
  molecule_id: string;
  concentration: number;
  concentration_unit: string;
}

export interface CreateMixtureData {
  name: string;
  description?: string;
  components: CreateMixtureComponentData[];
}

export interface UpdateMixtureData {
  name?: string;
  description?: string;
  is_cryoprotectant_mixture?: boolean;
}

const BASE_URL = process.env.NEXT_PUBLIC_API_URL || '/api';

export const mixtureService = {
  /**
   * Get a list of mixtures with pagination
   */
  async getMixtures(params: MixtureParams = {}): Promise<MixturesResponse> {
    try {
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
   */
  async getMixture(id: string): Promise<Mixture> {
    try {
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
   */
  async createMixture(data: CreateMixtureData): Promise<Mixture> {
    try {
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
   */
  async updateMixture(id: string, data: UpdateMixtureData): Promise<Mixture> {
    try {
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
   */
  async deleteMixture(id: string): Promise<void> {
    try {
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
   */
  async addComponent(mixtureId: string, data: CreateMixtureComponentData): Promise<MixtureComponent> {
    try {
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
   */
  async updateComponent(
    mixtureId: string, 
    componentId: string, 
    data: { concentration?: number; concentration_unit?: string }
  ): Promise<MixtureComponent> {
    try {
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
   */
  async removeComponent(mixtureId: string, componentId: string): Promise<void> {
    try {
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
   */
  async getCryoprotectionScore(mixtureId: string): Promise<{ score: number; details?: any }> {
    try {
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
   */
  async searchMixtures(query: string, limit: number = 10): Promise<Mixture[]> {
    try {
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
  }
};
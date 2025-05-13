import axios from 'axios';

export interface Molecule {
  id: string;
  name: string;
  smiles: string;
  pubchem_cid?: string;
  formula?: string;
  molecular_weight?: number;
  is_cryoprotectant?: boolean;
  created_at?: string;
  updated_at?: string;
}

export interface MoleculeProperty {
  id: string;
  molecule_id: string;
  property_name: string;
  property_value: string | number;
  property_unit?: string;
  calculation_method?: string;
  created_at?: string;
  updated_at?: string;
}

export interface MoleculesResponse {
  molecules: Molecule[];
  total: number;
  page: number;
  limit: number;
}

export interface MoleculeParams {
  page?: number;
  limit?: number;
  search?: string;
  sort_by?: string;
  sort_order?: 'asc' | 'desc';
  is_cryoprotectant?: boolean;
}

export interface MoleculeError {
  status: number;
  message: string;
}

const BASE_URL = process.env.NEXT_PUBLIC_API_URL || '/api';

export const moleculeService = {
  /**
   * Get a list of molecules with pagination
   */
  async getMolecules(params: MoleculeParams = {}): Promise<MoleculesResponse> {
    try {
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
   */
  async getMolecule(id: string): Promise<Molecule> {
    try {
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
   */
  async getMoleculeProperties(id: string): Promise<MoleculeProperty[]> {
    try {
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
   */
  async searchMolecules(query: string, limit: number = 10): Promise<Molecule[]> {
    try {
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
   */
  async importFromPubChem(cid: string): Promise<Molecule> {
    try {
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
  }
};
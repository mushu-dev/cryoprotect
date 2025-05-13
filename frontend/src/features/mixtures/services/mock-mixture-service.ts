import { Mixture, MixtureParams, MixturesResponse } from './mixture-service';

// Mock data for mixtures
const mockMixtures: Mixture[] = [
  {
    id: 'mix-1',
    name: 'Glycerol-DMSO Mixture',
    description: 'Standard cryoprotective mixture with glycerol and DMSO',
    is_cryoprotectant_mixture: true,
    cryoprotection_score: 8.5,
    total_concentration: 30,
    created_at: '2025-05-01T12:00:00Z',
    updated_at: '2025-05-01T12:00:00Z',
    components: [
      {
        id: 'comp-1',
        mixture_id: 'mix-1',
        molecule: {
          id: 'mol-1',
          name: 'Glycerol',
          smiles: 'C(C(CO)O)O',
          formula: 'C3H8O3',
          is_cryoprotectant: true,
        },
        concentration: 15,
        concentration_unit: '%',
      },
      {
        id: 'comp-2',
        mixture_id: 'mix-1',
        molecule: {
          id: 'mol-2',
          name: 'Dimethyl sulfoxide',
          smiles: 'CS(=O)C',
          formula: 'C2H6OS',
          is_cryoprotectant: true,
        },
        concentration: 15,
        concentration_unit: '%',
      },
    ],
  },
  {
    id: 'mix-2',
    name: 'Ethylene Glycol Solution',
    description: 'Simple ethylene glycol solution for cryopreservation',
    is_cryoprotectant_mixture: true,
    cryoprotection_score: 7.2,
    total_concentration: 20,
    created_at: '2025-05-01T12:00:00Z',
    updated_at: '2025-05-01T12:00:00Z',
    components: [
      {
        id: 'comp-3',
        mixture_id: 'mix-2',
        molecule: {
          id: 'mol-3',
          name: 'Ethylene glycol',
          smiles: 'C(CO)O',
          formula: 'C2H6O2',
          is_cryoprotectant: true,
        },
        concentration: 20,
        concentration_unit: '%',
      },
    ],
  },
  {
    id: 'mix-3',
    name: 'Trehalose Buffer',
    description: 'Buffered trehalose solution for protein preservation',
    is_cryoprotectant_mixture: true,
    cryoprotection_score: 6.8,
    total_concentration: 10,
    created_at: '2025-05-01T12:00:00Z',
    updated_at: '2025-05-01T12:00:00Z',
    components: [
      {
        id: 'comp-4',
        mixture_id: 'mix-3',
        molecule: {
          id: 'mol-5',
          name: 'Trehalose',
          smiles: 'C(C1C(C(C(C(O1)OC2C(C(C(C(O2)CO)O)O)O)O)O)O)O',
          formula: 'C12H22O11',
          is_cryoprotectant: true,
        },
        concentration: 10,
        concentration_unit: '%',
      },
    ],
  },
];

export const mockMixtureService = {
  /**
   * Get a list of mixtures with pagination
   */
  async getMixtures(params: MixtureParams = {}): Promise<MixturesResponse> {
    // Simulate API delay
    await new Promise((resolve) => setTimeout(resolve, 500));
    
    const page = params.page || 1;
    const limit = params.limit || 10;
    const startIndex = (page - 1) * limit;
    const endIndex = startIndex + limit;
    
    let filteredMixtures = [...mockMixtures];
    
    // Filter by search term
    if (params.search) {
      const searchTerm = params.search.toLowerCase();
      filteredMixtures = filteredMixtures.filter(
        (mix) => mix.name.toLowerCase().includes(searchTerm) || 
                 mix.description?.toLowerCase().includes(searchTerm)
      );
    }
    
    // Filter by cryoprotectant status
    if (params.is_cryoprotectant_mixture !== undefined) {
      filteredMixtures = filteredMixtures.filter(
        (mix) => mix.is_cryoprotectant_mixture === params.is_cryoprotectant_mixture
      );
    }
    
    // Sort
    if (params.sort_by) {
      const sortField = params.sort_by as keyof Mixture;
      const sortOrder = params.sort_order === 'desc' ? -1 : 1;
      
      filteredMixtures.sort((a, b) => {
        const aValue = a[sortField];
        const bValue = b[sortField];
        
        if (aValue === undefined) return sortOrder === 1 ? 1 : -1;
        if (bValue === undefined) return sortOrder === 1 ? -1 : 1;
        
        if (typeof aValue === 'string' && typeof bValue === 'string') {
          return sortOrder * aValue.localeCompare(bValue);
        }
        
        return sortOrder * ((aValue as number) - (bValue as number));
      });
    }
    
    // Paginate
    const paginatedMixtures = filteredMixtures.slice(startIndex, endIndex);
    
    return {
      mixtures: paginatedMixtures,
      total: filteredMixtures.length,
      page: page,
      limit: limit,
    };
  },

  /**
   * Get a single mixture by ID
   */
  async getMixture(id: string): Promise<Mixture> {
    // Simulate API delay
    await new Promise((resolve) => setTimeout(resolve, 300));
    
    const mixture = mockMixtures.find((mix) => mix.id === id);
    
    if (!mixture) {
      throw new Error(`Mixture with ID ${id} not found`);
    }
    
    return mixture;
  },

  /**
   * Get cryoprotection score for a mixture
   */
  async getCryoprotectionScore(id: string): Promise<{ score: number }> {
    // Simulate API delay
    await new Promise((resolve) => setTimeout(resolve, 300));
    
    const mixture = mockMixtures.find((mix) => mix.id === id);
    
    if (!mixture) {
      throw new Error(`Mixture with ID ${id} not found`);
    }
    
    return { score: mixture.cryoprotection_score || 0 };
  },

  /**
   * Create a new mixture
   */
  async createMixture(data: Partial<Mixture>): Promise<Mixture> {
    // Simulate API delay
    await new Promise((resolve) => setTimeout(resolve, 500));
    
    const newMixture: Mixture = {
      id: `mix-${mockMixtures.length + 1}`,
      name: data.name || 'New Mixture',
      description: data.description || '',
      is_cryoprotectant_mixture: data.is_cryoprotectant_mixture || false,
      total_concentration: data.total_concentration || 0,
      cryoprotection_score: data.cryoprotection_score || 0,
      created_at: new Date().toISOString(),
      updated_at: new Date().toISOString(),
      components: data.components || [],
    };
    
    return newMixture;
  },

  /**
   * Add component to a mixture
   */
  async addComponent(mixtureId: string, componentData: any): Promise<any> {
    // Simulate API delay
    await new Promise((resolve) => setTimeout(resolve, 300));
    
    const mixture = mockMixtures.find((mix) => mix.id === mixtureId);
    
    if (!mixture) {
      throw new Error(`Mixture with ID ${mixtureId} not found`);
    }
    
    const newComponent = {
      id: `comp-${Math.random().toString(36).substr(2, 9)}`,
      mixture_id: mixtureId,
      molecule: componentData.molecule,
      concentration: componentData.concentration || 0,
      concentration_unit: componentData.concentration_unit || '%',
    };
    
    return newComponent;
  }
};
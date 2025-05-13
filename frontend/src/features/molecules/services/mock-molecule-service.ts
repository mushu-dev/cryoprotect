import { Molecule, MoleculeParams, MoleculeProperty, MoleculesResponse } from './molecule-service';

// Mock data for molecules
const mockMolecules: Molecule[] = [
  {
    id: 'mol-1',
    name: 'Glycerol',
    smiles: 'C(C(CO)O)O',
    pubchem_cid: '753',
    formula: 'C3H8O3',
    molecular_weight: 92.09,
    is_cryoprotectant: true,
    created_at: '2025-05-01T12:00:00Z',
    updated_at: '2025-05-01T12:00:00Z',
  },
  {
    id: 'mol-2',
    name: 'Dimethyl sulfoxide',
    smiles: 'CS(=O)C',
    pubchem_cid: '679',
    formula: 'C2H6OS',
    molecular_weight: 78.13,
    is_cryoprotectant: true,
    created_at: '2025-05-01T12:00:00Z',
    updated_at: '2025-05-01T12:00:00Z',
  },
  {
    id: 'mol-3',
    name: 'Ethylene glycol',
    smiles: 'C(CO)O',
    pubchem_cid: '174',
    formula: 'C2H6O2',
    molecular_weight: 62.07,
    is_cryoprotectant: true,
    created_at: '2025-05-01T12:00:00Z',
    updated_at: '2025-05-01T12:00:00Z',
  },
  {
    id: 'mol-4',
    name: 'Propylene glycol',
    smiles: 'CC(CO)O',
    pubchem_cid: '1030',
    formula: 'C3H8O2',
    molecular_weight: 76.09,
    is_cryoprotectant: true,
    created_at: '2025-05-01T12:00:00Z',
    updated_at: '2025-05-01T12:00:00Z',
  },
  {
    id: 'mol-5',
    name: 'Trehalose',
    smiles: 'C(C1C(C(C(C(O1)OC2C(C(C(C(O2)CO)O)O)O)O)O)O)O',
    pubchem_cid: '7427',
    formula: 'C12H22O11',
    molecular_weight: 342.30,
    is_cryoprotectant: true,
    created_at: '2025-05-01T12:00:00Z',
    updated_at: '2025-05-01T12:00:00Z',
  },
];

// Mock data for properties
const mockProperties: Record<string, MoleculeProperty[]> = {
  'mol-1': [
    {
      id: 'prop-1-1',
      molecule_id: 'mol-1',
      property_name: 'LogP',
      property_value: -1.76,
      property_unit: '',
      calculation_method: 'Computed',
    },
    {
      id: 'prop-1-2',
      molecule_id: 'mol-1',
      property_name: 'Hydrogen Bond Donors',
      property_value: 3,
      calculation_method: 'Computed',
    },
    {
      id: 'prop-1-3',
      molecule_id: 'mol-1',
      property_name: 'Hydrogen Bond Acceptors',
      property_value: 3,
      calculation_method: 'Computed',
    },
  ],
  'mol-2': [
    {
      id: 'prop-2-1',
      molecule_id: 'mol-2',
      property_name: 'LogP',
      property_value: -1.35,
      property_unit: '',
      calculation_method: 'Computed',
    },
    {
      id: 'prop-2-2',
      molecule_id: 'mol-2',
      property_name: 'Hydrogen Bond Donors',
      property_value: 0,
      calculation_method: 'Computed',
    },
  ],
};

export const mockMoleculeService = {
  /**
   * Get a list of molecules with pagination
   */
  async getMolecules(params: MoleculeParams = {}): Promise<MoleculesResponse> {
    // Simulate API delay
    await new Promise((resolve) => setTimeout(resolve, 500));
    
    const page = params.page || 1;
    const limit = params.limit || 10;
    const startIndex = (page - 1) * limit;
    const endIndex = startIndex + limit;
    
    let filteredMolecules = [...mockMolecules];
    
    // Filter by search term
    if (params.search) {
      const searchTerm = params.search.toLowerCase();
      filteredMolecules = filteredMolecules.filter(
        (mol) => mol.name.toLowerCase().includes(searchTerm) || mol.formula?.toLowerCase().includes(searchTerm)
      );
    }
    
    // Filter by cryoprotectant status
    if (params.is_cryoprotectant !== undefined) {
      filteredMolecules = filteredMolecules.filter(
        (mol) => mol.is_cryoprotectant === params.is_cryoprotectant
      );
    }
    
    // Sort
    if (params.sort_by) {
      const sortField = params.sort_by as keyof Molecule;
      const sortOrder = params.sort_order === 'desc' ? -1 : 1;
      
      filteredMolecules.sort((a, b) => {
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
    const paginatedMolecules = filteredMolecules.slice(startIndex, endIndex);
    
    return {
      molecules: paginatedMolecules,
      total: filteredMolecules.length,
      page: page,
      limit: limit,
    };
  },

  /**
   * Get a single molecule by ID
   */
  async getMolecule(id: string): Promise<Molecule> {
    // Simulate API delay
    await new Promise((resolve) => setTimeout(resolve, 300));
    
    const molecule = mockMolecules.find((mol) => mol.id === id);
    
    if (!molecule) {
      throw new Error(`Molecule with ID ${id} not found`);
    }
    
    return molecule;
  },

  /**
   * Get properties for a molecule
   */
  async getMoleculeProperties(id: string): Promise<MoleculeProperty[]> {
    // Simulate API delay
    await new Promise((resolve) => setTimeout(resolve, 300));
    
    // Check if molecule exists
    if (!mockMolecules.find((mol) => mol.id === id)) {
      throw new Error(`Molecule with ID ${id} not found`);
    }
    
    // Return properties or empty array if none exist
    return mockProperties[id] || [];
  },

  /**
   * Search molecules by name or SMILES
   */
  async searchMolecules(query: string, limit: number = 10): Promise<Molecule[]> {
    // Simulate API delay
    await new Promise((resolve) => setTimeout(resolve, 300));
    
    if (!query) {
      return [];
    }
    
    const searchTerm = query.toLowerCase();
    const results = mockMolecules.filter(
      (mol) => 
        mol.name.toLowerCase().includes(searchTerm) || 
        mol.smiles.toLowerCase().includes(searchTerm) ||
        mol.formula?.toLowerCase().includes(searchTerm) ||
        mol.pubchem_cid?.includes(searchTerm)
    );
    
    return results.slice(0, limit);
  },

  /**
   * Import a molecule from PubChem by CID
   */
  async importFromPubChem(cid: string): Promise<Molecule> {
    // Simulate API delay
    await new Promise((resolve) => setTimeout(resolve, 1000));
    
    // Check if molecule already exists
    const existingMolecule = mockMolecules.find((mol) => mol.pubchem_cid === cid);
    
    if (existingMolecule) {
      return existingMolecule;
    }
    
    // Simulate creating a new molecule
    const newMolecule: Molecule = {
      id: `mol-${mockMolecules.length + 1}`,
      name: `Imported Molecule ${cid}`,
      smiles: 'CC(=O)OC1=CC=CC=C1C(=O)O',
      pubchem_cid: cid,
      formula: 'C9H8O4',
      molecular_weight: 180.16,
      is_cryoprotectant: false,
      created_at: new Date().toISOString(),
      updated_at: new Date().toISOString(),
    };
    
    // In a real implementation, we would add this to the list
    // mockMolecules.push(newMolecule);
    
    return newMolecule;
  }
};
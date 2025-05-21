const axios = require('axios');
const MockAdapter = require('axios-mock-adapter');

// Create a new instance of axios mock adapter
const mock = new MockAdapter(axios);

// Sample test data
const sampleMolecules = [
  {
    id: '1',
    name: 'Glycerol',
    smiles: 'C(C(CO)O)O',
    pubchem_cid: '753',
    formula: 'C3H8O3',
    molecular_weight: 92.09,
    is_cryoprotectant: true
  },
  {
    id: '2',
    name: 'Dimethyl sulfoxide',
    smiles: 'CS(=O)C',
    pubchem_cid: '679',
    formula: 'C2H6OS',
    molecular_weight: 78.13,
    is_cryoprotectant: true
  }
];

// Define a simple molecule service
const moleculeService = {
  async getMolecules() {
    try {
      const response = await axios.get('/api/v1/molecules');
      return response.data;
    } catch (error) {
      throw error;
    }
  },
  
  async getMolecule(id) {
    try {
      const response = await axios.get(`/api/v1/molecules/${id}`);
      return response.data;
    } catch (error) {
      throw error;
    }
  }
};

// Define a simple mixture service
const mixtureService = {
  async getMixtures() {
    try {
      const response = await axios.get('/api/v1/mixtures');
      return response.data;
    } catch (error) {
      throw error;
    }
  },
  
  async getMixture(id) {
    try {
      const response = await axios.get(`/api/v1/mixtures/${id}`);
      return response.data;
    } catch (error) {
      throw error;
    }
  }
};

describe('Molecule API Tests', () => {
  beforeEach(() => {
    mock.reset();
  });

  test('getMolecules should return a list of molecules', async () => {
    const paginatedResponse = {
      molecules: sampleMolecules,
      total: 2,
      page: 1,
      limit: 10
    };
    
    mock.onGet('/api/v1/molecules').reply(200, paginatedResponse);
    
    const result = await moleculeService.getMolecules();
    
    expect(result).toEqual(paginatedResponse);
    expect(result.molecules.length).toBe(2);
    expect(result.total).toBe(2);
  });

  test('getMolecule should return a single molecule by ID', async () => {
    mock.onGet('/api/v1/molecules/1').reply(200, sampleMolecules[0]);
    
    const result = await moleculeService.getMolecule('1');
    
    expect(result).toEqual(sampleMolecules[0]);
    expect(result.name).toBe('Glycerol');
  });

  test('getMolecules should handle error', async () => {
    mock.onGet('/api/v1/molecules').reply(500, { message: 'Server error' });
    
    await expect(moleculeService.getMolecules()).rejects.toThrow();
  });
});

describe('Mixture API Tests', () => {
  beforeEach(() => {
    mock.reset();
  });

  test('getMixtures should return a list of mixtures', async () => {
    const paginatedResponse = {
      mixtures: [
        {
          id: '1',
          name: 'Glycerol-DMSO Mixture',
          components: [
            { molecule_id: '1', concentration: 10 },
            { molecule_id: '2', concentration: 5 }
          ]
        }
      ],
      total: 1,
      page: 1,
      limit: 10
    };
    
    mock.onGet('/api/v1/mixtures').reply(200, paginatedResponse);
    
    const result = await mixtureService.getMixtures();
    
    expect(result).toEqual(paginatedResponse);
    expect(result.mixtures.length).toBe(1);
    expect(result.total).toBe(1);
  });

  test('getMixture should return a single mixture by ID', async () => {
    const mixture = {
      id: '1',
      name: 'Glycerol-DMSO Mixture',
      components: [
        { molecule_id: '1', concentration: 10 },
        { molecule_id: '2', concentration: 5 }
      ]
    };
    
    mock.onGet('/api/v1/mixtures/1').reply(200, mixture);
    
    const result = await mixtureService.getMixture('1');
    
    expect(result).toEqual(mixture);
    expect(result.name).toBe('Glycerol-DMSO Mixture');
  });

  test('getMixtures should handle error', async () => {
    mock.onGet('/api/v1/mixtures').reply(500, { message: 'Server error' });
    
    await expect(mixtureService.getMixtures()).rejects.toThrow();
  });
});
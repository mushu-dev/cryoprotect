const axios = require('axios');
const MockAdapter = require('axios-mock-adapter');
const { moleculeService } = require('../../src/features/molecules/services/molecule-service');
const { mixtureService } = require('../../src/features/mixtures/services/mixture-service');

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

const sampleMoleculeProperties = [
  {
    id: 'prop1',
    molecule_id: '1',
    property_name: 'LogP',
    property_value: -1.76,
    property_unit: null,
    calculation_method: 'RDKit'
  },
  {
    id: 'prop2',
    molecule_id: '1',
    property_name: 'Molecular Weight',
    property_value: 92.09,
    property_unit: 'g/mol',
    calculation_method: 'RDKit'
  }
];

const sampleMixtures = [
  {
    id: '1',
    name: 'Glycerol-DMSO Mixture',
    description: 'A common cryoprotectant mixture',
    is_cryoprotectant_mixture: true,
    components: [
      {
        id: 'comp1',
        mixture_id: '1',
        molecule_id: '1',
        molecule: sampleMolecules[0],
        concentration: 10,
        concentration_unit: '%w/v'
      },
      {
        id: 'comp2',
        mixture_id: '1',
        molecule_id: '2',
        molecule: sampleMolecules[1],
        concentration: 5,
        concentration_unit: '%w/v'
      }
    ],
    cryoprotection_score: 85
  }
];

describe('Molecule API Tests', () => {
  beforeEach(() => {
    mock.reset();
  });

  test('getMolecules should return a list of molecules with pagination', async () => {
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

  test('getMoleculeProperties should return properties for a molecule', async () => {
    mock.onGet('/api/v1/molecules/1/properties').reply(200, { properties: sampleMoleculeProperties });
    
    const result = await moleculeService.getMoleculeProperties('1');
    
    expect(result).toEqual(sampleMoleculeProperties);
    expect(result.length).toBe(2);
    expect(result[0].property_name).toBe('LogP');
  });

  test('searchMolecules should return molecules matching query', async () => {
    mock.onGet('/api/v1/molecules/search').reply(200, { results: [sampleMolecules[0]] });
    
    const result = await moleculeService.searchMolecules('Glyc');
    
    expect(result).toEqual([sampleMolecules[0]]);
    expect(result.length).toBe(1);
    expect(result[0].name).toBe('Glycerol');
  });

  test('importFromPubChem should import a molecule by CID', async () => {
    mock.onPost('/api/v1/molecules/import/pubchem').reply(200, sampleMolecules[0]);
    
    const result = await moleculeService.importFromPubChem('753');
    
    expect(result).toEqual(sampleMolecules[0]);
    expect(result.pubchem_cid).toBe('753');
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

  test('getMixtures should return a list of mixtures with pagination', async () => {
    const paginatedResponse = {
      mixtures: sampleMixtures,
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
    mock.onGet('/api/v1/mixtures/1').reply(200, sampleMixtures[0]);
    
    const result = await mixtureService.getMixture('1');
    
    expect(result).toEqual(sampleMixtures[0]);
    expect(result.name).toBe('Glycerol-DMSO Mixture');
    expect(result.components.length).toBe(2);
  });

  test('createMixture should create a new mixture', async () => {
    const createData = {
      name: 'New Test Mixture',
      description: 'A test mixture',
      components: [
        {
          molecule_id: '1',
          concentration: 15,
          concentration_unit: '%w/v'
        }
      ]
    };
    
    const newMixture = {
      id: '3',
      name: 'New Test Mixture',
      description: 'A test mixture',
      is_cryoprotectant_mixture: false,
      components: [
        {
          id: 'newcomp1',
          mixture_id: '3',
          molecule_id: '1',
          concentration: 15,
          concentration_unit: '%w/v'
        }
      ]
    };
    
    mock.onPost('/api/v1/mixtures').reply(201, newMixture);
    
    const result = await mixtureService.createMixture(createData);
    
    expect(result).toEqual(newMixture);
    expect(result.id).toBe('3');
    expect(result.name).toBe('New Test Mixture');
  });

  test('updateMixture should update a mixture', async () => {
    const updateData = {
      name: 'Updated Mixture Name',
      is_cryoprotectant_mixture: true
    };
    
    const updatedMixture = {
      ...sampleMixtures[0],
      name: 'Updated Mixture Name'
    };
    
    mock.onPatch('/api/v1/mixtures/1').reply(200, updatedMixture);
    
    const result = await mixtureService.updateMixture('1', updateData);
    
    expect(result).toEqual(updatedMixture);
    expect(result.name).toBe('Updated Mixture Name');
  });

  test('deleteMixture should delete a mixture', async () => {
    mock.onDelete('/api/v1/mixtures/1').reply(204);
    
    await expect(mixtureService.deleteMixture('1')).resolves.not.toThrow();
  });

  test('addComponent should add a component to a mixture', async () => {
    const componentData = {
      molecule_id: '2',
      concentration: 8,
      concentration_unit: '%w/v'
    };
    
    const newComponent = {
      id: 'newcomp2',
      mixture_id: '1',
      molecule_id: '2',
      concentration: 8,
      concentration_unit: '%w/v'
    };
    
    mock.onPost('/api/v1/mixtures/1/components').reply(201, newComponent);
    
    const result = await mixtureService.addComponent('1', componentData);
    
    expect(result).toEqual(newComponent);
    expect(result.concentration).toBe(8);
  });

  test('updateComponent should update a mixture component', async () => {
    const updateData = {
      concentration: 12,
      concentration_unit: '%v/v'
    };
    
    const updatedComponent = {
      id: 'comp1',
      mixture_id: '1',
      molecule_id: '1',
      concentration: 12,
      concentration_unit: '%v/v'
    };
    
    mock.onPatch('/api/v1/mixtures/1/components/comp1').reply(200, updatedComponent);
    
    const result = await mixtureService.updateComponent('1', 'comp1', updateData);
    
    expect(result).toEqual(updatedComponent);
    expect(result.concentration).toBe(12);
    expect(result.concentration_unit).toBe('%v/v');
  });

  test('removeComponent should remove a component from a mixture', async () => {
    mock.onDelete('/api/v1/mixtures/1/components/comp1').reply(204);
    
    await expect(mixtureService.removeComponent('1', 'comp1')).resolves.not.toThrow();
  });

  test('getCryoprotectionScore should return score for a mixture', async () => {
    const scoreData = {
      score: 85,
      details: {
        hydrogen_bonding_score: 92,
        viscosity_score: 78,
        toxicity_score: 65
      }
    };
    
    mock.onGet('/api/v1/mixtures/1/cryoprotection-score').reply(200, scoreData);
    
    const result = await mixtureService.getCryoprotectionScore('1');
    
    expect(result).toEqual(scoreData);
    expect(result.score).toBe(85);
    expect(result.details.hydrogen_bonding_score).toBe(92);
  });

  test('searchMixtures should return mixtures matching query', async () => {
    mock.onGet('/api/v1/mixtures/search').reply(200, { results: [sampleMixtures[0]] });
    
    const result = await mixtureService.searchMixtures('Glycerol');
    
    expect(result).toEqual([sampleMixtures[0]]);
    expect(result.length).toBe(1);
    expect(result[0].name).toBe('Glycerol-DMSO Mixture');
  });

  test('getMixtures should handle error', async () => {
    mock.onGet('/api/v1/mixtures').reply(500, { message: 'Server error' });
    
    await expect(mixtureService.getMixtures()).rejects.toThrow();
  });
});

describe('Integration Tests for Molecule and Mixture APIs', () => {
  beforeEach(() => {
    mock.reset();
  });

  test('Creating a mixture with molecules from API', async () => {
    // First fetch molecules
    mock.onGet('/api/v1/molecules').reply(200, {
      molecules: sampleMolecules,
      total: 2,
      page: 1,
      limit: 10
    });
    
    // Then create a mixture with those molecules
    const createData = {
      name: 'Integration Test Mixture',
      description: 'Created in integration test',
      components: [
        {
          molecule_id: '1',
          concentration: 15,
          concentration_unit: '%w/v'
        },
        {
          molecule_id: '2',
          concentration: 7.5,
          concentration_unit: '%w/v'
        }
      ]
    };
    
    const newMixture = {
      id: '4',
      name: 'Integration Test Mixture',
      description: 'Created in integration test',
      is_cryoprotectant_mixture: false,
      components: [
        {
          id: 'intcomp1',
          mixture_id: '4',
          molecule_id: '1',
          concentration: 15,
          concentration_unit: '%w/v'
        },
        {
          id: 'intcomp2',
          mixture_id: '4',
          molecule_id: '2',
          concentration: 7.5,
          concentration_unit: '%w/v'
        }
      ]
    };
    
    mock.onPost('/api/v1/mixtures').reply(201, newMixture);
    
    // Get the molecules
    const moleculesResult = await moleculeService.getMolecules();
    
    // Create the mixture with those molecules
    const mixtureResult = await mixtureService.createMixture(createData);
    
    // Verify the results
    expect(moleculesResult.molecules.length).toBe(2);
    expect(mixtureResult.id).toBe('4');
    expect(mixtureResult.components.length).toBe(2);
    expect(mixtureResult.components[0].molecule_id).toBe('1');
    expect(mixtureResult.components[1].molecule_id).toBe('2');
  });

  test('Getting molecules and their properties', async () => {
    // Set up mocks
    mock.onGet('/api/v1/molecules').reply(200, {
      molecules: sampleMolecules,
      total: 2,
      page: 1,
      limit: 10
    });
    
    mock.onGet('/api/v1/molecules/1/properties').reply(200, { 
      properties: sampleMoleculeProperties 
    });
    
    // Get molecules and then properties for the first molecule
    const moleculesResult = await moleculeService.getMolecules();
    const firstMoleculeId = moleculesResult.molecules[0].id;
    const propertiesResult = await moleculeService.getMoleculeProperties(firstMoleculeId);
    
    // Verify the results
    expect(moleculesResult.molecules.length).toBe(2);
    expect(firstMoleculeId).toBe('1');
    expect(propertiesResult.length).toBe(2);
    expect(propertiesResult[0].property_name).toBe('LogP');
    expect(propertiesResult[1].property_name).toBe('Molecular Weight');
  });

  test('Creating a mixture and getting its cryoprotection score', async () => {
    // Set up mocks
    const createData = {
      name: 'Score Test Mixture',
      description: 'Created for score testing',
      components: [
        {
          molecule_id: '1',
          concentration: 20,
          concentration_unit: '%w/v'
        }
      ]
    };
    
    const newMixture = {
      id: '5',
      name: 'Score Test Mixture',
      description: 'Created for score testing',
      is_cryoprotectant_mixture: false,
      components: [
        {
          id: 'scorecomp1',
          mixture_id: '5',
          molecule_id: '1',
          concentration: 20,
          concentration_unit: '%w/v'
        }
      ]
    };
    
    const scoreData = {
      score: 78,
      details: {
        hydrogen_bonding_score: 85,
        viscosity_score: 70,
        toxicity_score: 60
      }
    };
    
    mock.onPost('/api/v1/mixtures').reply(201, newMixture);
    mock.onGet('/api/v1/mixtures/5/cryoprotection-score').reply(200, scoreData);
    
    // Create the mixture
    const mixtureResult = await mixtureService.createMixture(createData);
    
    // Get the cryoprotection score
    const scoreResult = await mixtureService.getCryoprotectionScore(mixtureResult.id);
    
    // Verify the results
    expect(mixtureResult.id).toBe('5');
    expect(scoreResult.score).toBe(78);
    expect(scoreResult.details.hydrogen_bonding_score).toBe(85);
  });
});
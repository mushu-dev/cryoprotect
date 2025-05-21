/**
 * Tests for molecule functions
 */

import { createMockContext } from './setup';
import { createMolecule, createMoleculesBatch } from '../../convex/molecules/create';
import { getMolecule, searchMolecules } from '../../convex/molecules/query';
import { updateMolecule, deprecateMolecule } from '../../convex/molecules/update';
import { deleteMolecule } from '../../convex/molecules/delete';
import { isValidMoleculeName, isValidPubchemCid, isValidInchiKey } from '../../convex/molecules/validation';

describe('Molecule Validation', () => {
  test('isValidMoleculeName validates correctly', () => {
    expect(isValidMoleculeName('Water')).toBe(true);
    expect(isValidMoleculeName('')).toBe(false);
    expect(isValidMoleculeName('  ')).toBe(false);
  });

  test('isValidPubchemCid validates correctly', () => {
    expect(isValidPubchemCid('962')).toBe(true);
    expect(isValidPubchemCid('abc')).toBe(false);
    expect(isValidPubchemCid('')).toBe(false);
  });

  test('isValidInchiKey validates correctly', () => {
    expect(isValidInchiKey('XLYOFNOQVPJJNP-UHFFFAOYSA-N')).toBe(true);
    expect(isValidInchiKey('invalid')).toBe(false);
    expect(isValidInchiKey('')).toBe(false);
  });
});

describe('Molecule CRUD Functions', () => {
  // Setup mock context and data
  const ctx = createMockContext();
  const mockMoleculeId = '01234567890123456789012345678901';
  const mockMolecule = {
    _id: { id: mockMoleculeId, tableName: 'molecules' },
    _creationTime: Date.now(),
    name: 'Water',
    pubchemCid: '962',
    canonicalSmiles: 'O',
    inchiKey: 'XLYOFNOQVPJJNP-UHFFFAOYSA-N',
    formula: 'H2O',
    status: 'active',
    consolidated: false,
    createdAt: Date.now(),
    updatedAt: Date.now(),
  };

  beforeEach(() => {
    // Reset mock functions
    jest.resetAllMocks();
    
    // Setup default mocks
    ctx.db.get.mockResolvedValue(mockMolecule);
    ctx.db.insert.mockResolvedValue({ id: mockMoleculeId, tableName: 'molecules' });
    ctx.db.patch.mockResolvedValue({ ...mockMolecule, name: 'Updated Water' });
    
    // Mock query results
    const mockQueryResult = {
      page: [mockMolecule],
      continueCursor: null,
    };
    
    ctx.db.query().filter().collect.mockResolvedValue([mockMolecule]);
    ctx.db.query().withIndex().eq().collect.mockResolvedValue([mockMolecule]);
    ctx.db.query().collect.mockResolvedValue([mockMolecule]);
  });

  describe('createMolecule', () => {
    test('creates a molecule with valid input', async () => {
      const input = {
        name: 'Water',
        pubchemCid: '962',
        canonicalSmiles: 'O',
        inchiKey: 'XLYOFNOQVPJJNP-UHFFFAOYSA-N',
        formula: 'H2O',
      };
      
      // Mock that no duplicates are found
      ctx.db.query().withIndex().eq().collect.mockResolvedValue([]);
      
      const result = await createMolecule.handler(ctx as any, { input });
      
      expect(ctx.db.insert).toHaveBeenCalledWith('molecules', expect.objectContaining({
        name: 'Water',
        pubchemCid: '962',
        status: 'active',
      }));
      expect(result).toEqual(mockMolecule);
    });

    test('throws error on duplicate molecule', async () => {
      const input = {
        name: 'Water',
        pubchemCid: '962',
      };
      
      // Mock that a duplicate is found
      ctx.db.query().withIndex().eq().collect.mockResolvedValue([mockMolecule]);
      
      await expect(createMolecule.handler(ctx as any, { input }))
        .rejects.toThrow('Duplicate molecule found');
    });
  });

  describe('getMolecule', () => {
    test('retrieves a molecule by ID', async () => {
      const result = await getMolecule.handler(ctx as any, { 
        id: { id: mockMoleculeId, tableName: 'molecules' },
        includeProperties: false,
      });
      
      expect(ctx.db.get).toHaveBeenCalledWith({ id: mockMoleculeId, tableName: 'molecules' });
      expect(result).toEqual(mockMolecule);
    });

    test('returns null for non-existent molecule', async () => {
      ctx.db.get.mockResolvedValue(null);
      
      const result = await getMolecule.handler(ctx as any, { 
        id: { id: 'nonexistent', tableName: 'molecules' },
        includeProperties: false,
      });
      
      expect(result).toBeNull();
    });
  });

  describe('updateMolecule', () => {
    test('updates a molecule with valid input', async () => {
      const input = {
        name: 'Updated Water',
      };
      
      const result = await updateMolecule.handler(ctx as any, { 
        id: { id: mockMoleculeId, tableName: 'molecules' },
        input,
      });
      
      expect(ctx.db.patch).toHaveBeenCalledWith(
        { id: mockMoleculeId, tableName: 'molecules' },
        expect.objectContaining({
          name: 'Updated Water',
        })
      );
      expect(result.name).toBe('Updated Water');
    });

    test('throws error when updating consolidated molecule', async () => {
      // Mock a consolidated molecule
      ctx.db.get.mockResolvedValue({
        ...mockMolecule,
        consolidated: true,
      });
      
      await expect(updateMolecule.handler(ctx as any, { 
        id: { id: mockMoleculeId, tableName: 'molecules' },
        input: { name: 'Should Fail' },
      })).rejects.toThrow('Cannot update a consolidated molecule');
    });
  });

  describe('deleteMolecule', () => {
    test('deletes a molecule with no references', async () => {
      // Mock no references to this molecule
      ctx.db.query().withIndex().first.mockResolvedValue(null);
      
      const result = await deleteMolecule.handler(ctx as any, { 
        id: { id: mockMoleculeId, tableName: 'molecules' },
      });
      
      expect(ctx.db.delete).toHaveBeenCalledWith({ id: mockMoleculeId, tableName: 'molecules' });
      expect(result.success).toBe(true);
    });

    test('throws error when deleting molecule with references', async () => {
      // Mock that there are properties referencing this molecule
      ctx.db.query().withIndex().first.mockResolvedValue({ _id: 'property1' });
      
      await expect(deleteMolecule.handler(ctx as any, { 
        id: { id: mockMoleculeId, tableName: 'molecules' },
      })).rejects.toThrow('Cannot delete molecule with properties');
    });

    test('deletes a molecule with references when force=true', async () => {
      // Mock that there are properties referencing this molecule
      ctx.db.query().withIndex().collect.mockResolvedValue([{ _id: 'property1' }]);
      
      const result = await deleteMolecule.handler(ctx as any, { 
        id: { id: mockMoleculeId, tableName: 'molecules' },
        force: true,
      });
      
      // Should delete all references and the molecule
      expect(ctx.db.delete).toHaveBeenCalled();
      expect(result.success).toBe(true);
    });
  });
});
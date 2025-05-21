/**
 * Molecules hook
 * Provides React Query hooks for molecule data access
 */
import { useQuery, useMutation, useQueryClient } from '@tanstack/react-query';
import { getMoleculeService, markApiConnectionFailed } from '../services/service-factory';

export const MOLECULES_QUERY_KEY = 'molecules';

/**
 * Hook to get a list of molecules with optional filtering
 * @param {Object} params - Query parameters
 * @returns {Object} React Query result
 */
export function useMolecules(params = {}) {
  const moleculeService = getMoleculeService();
  
  return useQuery({
    queryKey: [MOLECULES_QUERY_KEY, params],
    queryFn: async () => {
      try {
        return await moleculeService.getMolecules(params);
      } catch (error) {
        // If the API call fails, mark the connection as failed for future mock data use
        markApiConnectionFailed();
        throw error;
      }
    },
    staleTime: 5 * 60 * 1000, // 5 minutes
  });
}

/**
 * Hook to get a single molecule by ID
 * @param {string} id - Molecule ID
 * @returns {Object} React Query result
 */
export function useMolecule(id) {
  const moleculeService = getMoleculeService();
  
  return useQuery({
    queryKey: [MOLECULES_QUERY_KEY, id],
    queryFn: async () => {
      try {
        return await moleculeService.getMolecule(id);
      } catch (error) {
        markApiConnectionFailed();
        throw error;
      }
    },
    staleTime: 5 * 60 * 1000, // 5 minutes
    enabled: !!id,
  });
}

/**
 * Hook to get properties for a molecule
 * @param {string} id - Molecule ID
 * @returns {Object} React Query result
 */
export function useMoleculeProperties(id) {
  const moleculeService = getMoleculeService();
  
  return useQuery({
    queryKey: [MOLECULES_QUERY_KEY, id, 'properties'],
    queryFn: async () => {
      try {
        return await moleculeService.getMoleculeProperties(id);
      } catch (error) {
        markApiConnectionFailed();
        throw error;
      }
    },
    staleTime: 5 * 60 * 1000, // 5 minutes
    enabled: !!id,
  });
}

/**
 * Hook to search molecules by query string
 * @param {string} query - Search query
 * @param {number} limit - Maximum number of results
 * @returns {Object} React Query result
 */
export function useSearchMolecules(query, limit = 10) {
  const moleculeService = getMoleculeService();
  
  return useQuery({
    queryKey: [MOLECULES_QUERY_KEY, 'search', query, limit],
    queryFn: async () => {
      try {
        return await moleculeService.searchMolecules(query, limit);
      } catch (error) {
        markApiConnectionFailed();
        throw error;
      }
    },
    staleTime: 5 * 60 * 1000, // 5 minutes
    enabled: query.length > 2, // Only search when query has 3+ characters
  });
}

/**
 * Hook to import a molecule from PubChem
 * @returns {Object} React Query mutation
 */
export function useImportFromPubChem() {
  const queryClient = useQueryClient();
  const moleculeService = getMoleculeService();
  
  return useMutation({
    mutationFn: async (cid) => {
      try {
        return await moleculeService.importFromPubChem(cid);
      } catch (error) {
        markApiConnectionFailed();
        throw error;
      }
    },
    onSuccess: (data) => {
      // Invalidate molecules query to refetch list after import
      queryClient.invalidateQueries({ queryKey: [MOLECULES_QUERY_KEY] });
      // Update cache with the new molecule
      queryClient.setQueryData([MOLECULES_QUERY_KEY, data.id], data);
    },
  });
}
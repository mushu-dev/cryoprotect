import { useQuery, useMutation, useQueryClient } from '@tanstack/react-query';
import { getMoleculeService, markApiConnectionFailed } from '@/features/common/api/service-factory';
import { 
  type Molecule, 
  type MoleculeParams,
  type MoleculeProperty 
} from '../services/molecule-service';

export const MOLECULES_QUERY_KEY = 'molecules';

export function useMolecules(params: MoleculeParams = {}) {
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

export function useMolecule(id: string) {
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

export function useMoleculeProperties(id: string) {
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

export function useSearchMolecules(query: string, limit: number = 10) {
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

export function useImportFromPubChem() {
  const queryClient = useQueryClient();
  const moleculeService = getMoleculeService();
  
  return useMutation({
    mutationFn: async (cid: string) => {
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
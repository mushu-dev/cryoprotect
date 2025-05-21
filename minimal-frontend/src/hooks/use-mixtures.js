/**
 * Mixtures hook
 * Provides React Query hooks for mixture data access
 */
import { useQuery, useMutation, useQueryClient } from '@tanstack/react-query';
import { getMixtureService, markApiConnectionFailed } from '../services/service-factory';

export const MIXTURES_QUERY_KEY = 'mixtures';

/**
 * Hook to get a list of mixtures with optional filtering
 * @param {Object} params - Query parameters
 * @returns {Object} React Query result
 */
export function useMixtures(params = {}) {
  const mixtureService = getMixtureService();
  
  return useQuery({
    queryKey: [MIXTURES_QUERY_KEY, params],
    queryFn: async () => {
      try {
        return await mixtureService.getMixtures(params);
      } catch (error) {
        markApiConnectionFailed();
        throw error;
      }
    },
    staleTime: 5 * 60 * 1000, // 5 minutes
  });
}

/**
 * Hook to get a single mixture by ID
 * @param {string} id - Mixture ID
 * @returns {Object} React Query result
 */
export function useMixture(id) {
  const mixtureService = getMixtureService();
  
  return useQuery({
    queryKey: [MIXTURES_QUERY_KEY, id],
    queryFn: async () => {
      try {
        return await mixtureService.getMixture(id);
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
 * Hook to create a new mixture
 * @returns {Object} React Query mutation
 */
export function useCreateMixture() {
  const queryClient = useQueryClient();
  const mixtureService = getMixtureService();
  
  return useMutation({
    mutationFn: async (data) => {
      try {
        return await mixtureService.createMixture(data);
      } catch (error) {
        markApiConnectionFailed();
        throw error;
      }
    },
    onSuccess: (data) => {
      // Invalidate mixtures query to refetch list after creation
      queryClient.invalidateQueries({ queryKey: [MIXTURES_QUERY_KEY] });
      // Update cache with the new mixture
      queryClient.setQueryData([MIXTURES_QUERY_KEY, data.id], data);
    },
  });
}

/**
 * Hook to update an existing mixture
 * @returns {Object} React Query mutation
 */
export function useUpdateMixture() {
  const queryClient = useQueryClient();
  const mixtureService = getMixtureService();
  
  return useMutation({
    mutationFn: ({ id, data }) => mixtureService.updateMixture(id, data),
    onSuccess: (data) => {
      // Invalidate mixtures query to refetch list after update
      queryClient.invalidateQueries({ queryKey: [MIXTURES_QUERY_KEY] });
      // Update cache with the updated mixture
      queryClient.setQueryData([MIXTURES_QUERY_KEY, data.id], data);
    },
  });
}

/**
 * Hook to delete a mixture
 * @returns {Object} React Query mutation
 */
export function useDeleteMixture() {
  const queryClient = useQueryClient();
  const mixtureService = getMixtureService();
  
  return useMutation({
    mutationFn: (id) => mixtureService.deleteMixture(id),
    onSuccess: (_, id) => {
      // Invalidate mixtures query to refetch list after deletion
      queryClient.invalidateQueries({ queryKey: [MIXTURES_QUERY_KEY] });
      // Remove deleted mixture from cache
      queryClient.removeQueries({ queryKey: [MIXTURES_QUERY_KEY, id] });
    },
  });
}

/**
 * Hook to add a component to a mixture
 * @returns {Object} React Query mutation
 */
export function useAddComponent() {
  const queryClient = useQueryClient();
  const mixtureService = getMixtureService();
  
  return useMutation({
    mutationFn: ({ mixtureId, data }) => mixtureService.addComponent(mixtureId, data),
    onSuccess: (_, { mixtureId }) => {
      // Invalidate specific mixture query to refetch after update
      queryClient.invalidateQueries({ queryKey: [MIXTURES_QUERY_KEY, mixtureId] });
    },
  });
}

/**
 * Hook to update a component in a mixture
 * @returns {Object} React Query mutation
 */
export function useUpdateComponent() {
  const queryClient = useQueryClient();
  const mixtureService = getMixtureService();
  
  return useMutation({
    mutationFn: ({ mixtureId, componentId, data }) => 
      mixtureService.updateComponent(mixtureId, componentId, data),
    onSuccess: (_, { mixtureId }) => {
      // Invalidate specific mixture query to refetch after update
      queryClient.invalidateQueries({ queryKey: [MIXTURES_QUERY_KEY, mixtureId] });
    },
  });
}

/**
 * Hook to remove a component from a mixture
 * @returns {Object} React Query mutation
 */
export function useRemoveComponent() {
  const queryClient = useQueryClient();
  const mixtureService = getMixtureService();
  
  return useMutation({
    mutationFn: ({ mixtureId, componentId }) => 
      mixtureService.removeComponent(mixtureId, componentId),
    onSuccess: (_, { mixtureId }) => {
      // Invalidate specific mixture query to refetch after update
      queryClient.invalidateQueries({ queryKey: [MIXTURES_QUERY_KEY, mixtureId] });
    },
  });
}

/**
 * Hook to get the cryoprotection score for a mixture
 * @param {string} mixtureId - Mixture ID
 * @returns {Object} React Query result
 */
export function useCryoprotectionScore(mixtureId) {
  const mixtureService = getMixtureService();
  
  return useQuery({
    queryKey: [MIXTURES_QUERY_KEY, mixtureId, 'cryoprotection-score'],
    queryFn: async () => {
      try {
        return await mixtureService.getCryoprotectionScore(mixtureId);
      } catch (error) {
        markApiConnectionFailed();
        throw error;
      }
    },
    staleTime: 5 * 60 * 1000, // 5 minutes
    enabled: !!mixtureId,
  });
}

/**
 * Hook to search mixtures by query string
 * @param {string} query - Search query
 * @param {number} limit - Maximum number of results
 * @returns {Object} React Query result
 */
export function useSearchMixtures(query, limit = 10) {
  const mixtureService = getMixtureService();
  
  return useQuery({
    queryKey: [MIXTURES_QUERY_KEY, 'search', query, limit],
    queryFn: async () => {
      try {
        return await mixtureService.searchMixtures(query, limit);
      } catch (error) {
        markApiConnectionFailed();
        throw error;
      }
    },
    staleTime: 5 * 60 * 1000, // 5 minutes
    enabled: query.length > 2, // Only search when query has 3+ characters
  });
}
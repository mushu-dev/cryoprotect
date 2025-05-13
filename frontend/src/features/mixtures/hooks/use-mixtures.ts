import { useQuery, useMutation, useQueryClient } from '@tanstack/react-query';
import { 
  mixtureService,
  type Mixture,
  type MixtureParams,
  type CreateMixtureData,
  type UpdateMixtureData,
  type CreateMixtureComponentData,
  type MixtureComponent
} from '../services/mixture-service';

export const MIXTURES_QUERY_KEY = 'mixtures';

export function useMixtures(params: MixtureParams = {}) {
  return useQuery({
    queryKey: [MIXTURES_QUERY_KEY, params],
    queryFn: () => mixtureService.getMixtures(params),
    staleTime: 5 * 60 * 1000, // 5 minutes
  });
}

export function useMixture(id: string) {
  return useQuery({
    queryKey: [MIXTURES_QUERY_KEY, id],
    queryFn: () => mixtureService.getMixture(id),
    staleTime: 5 * 60 * 1000, // 5 minutes
    enabled: !!id,
  });
}

export function useCreateMixture() {
  const queryClient = useQueryClient();
  
  return useMutation({
    mutationFn: (data: CreateMixtureData) => mixtureService.createMixture(data),
    onSuccess: (data) => {
      // Invalidate mixtures query to refetch list after creation
      queryClient.invalidateQueries({ queryKey: [MIXTURES_QUERY_KEY] });
      // Update cache with the new mixture
      queryClient.setQueryData([MIXTURES_QUERY_KEY, data.id], data);
    },
  });
}

export function useUpdateMixture() {
  const queryClient = useQueryClient();
  
  return useMutation({
    mutationFn: ({ id, data }: { id: string; data: UpdateMixtureData }) => 
      mixtureService.updateMixture(id, data),
    onSuccess: (data) => {
      // Invalidate mixtures query to refetch list after update
      queryClient.invalidateQueries({ queryKey: [MIXTURES_QUERY_KEY] });
      // Update cache with the updated mixture
      queryClient.setQueryData([MIXTURES_QUERY_KEY, data.id], data);
    },
  });
}

export function useDeleteMixture() {
  const queryClient = useQueryClient();
  
  return useMutation({
    mutationFn: (id: string) => mixtureService.deleteMixture(id),
    onSuccess: (_, id) => {
      // Invalidate mixtures query to refetch list after deletion
      queryClient.invalidateQueries({ queryKey: [MIXTURES_QUERY_KEY] });
      // Remove deleted mixture from cache
      queryClient.removeQueries({ queryKey: [MIXTURES_QUERY_KEY, id] });
    },
  });
}

export function useAddComponent() {
  const queryClient = useQueryClient();
  
  return useMutation({
    mutationFn: ({ mixtureId, data }: { mixtureId: string; data: CreateMixtureComponentData }) => 
      mixtureService.addComponent(mixtureId, data),
    onSuccess: (_, { mixtureId }) => {
      // Invalidate specific mixture query to refetch after update
      queryClient.invalidateQueries({ queryKey: [MIXTURES_QUERY_KEY, mixtureId] });
    },
  });
}

export function useUpdateComponent() {
  const queryClient = useQueryClient();
  
  return useMutation({
    mutationFn: ({ 
      mixtureId, 
      componentId, 
      data 
    }: { 
      mixtureId: string; 
      componentId: string; 
      data: { concentration?: number; concentration_unit?: string } 
    }) => mixtureService.updateComponent(mixtureId, componentId, data),
    onSuccess: (_, { mixtureId }) => {
      // Invalidate specific mixture query to refetch after update
      queryClient.invalidateQueries({ queryKey: [MIXTURES_QUERY_KEY, mixtureId] });
    },
  });
}

export function useRemoveComponent() {
  const queryClient = useQueryClient();
  
  return useMutation({
    mutationFn: ({ mixtureId, componentId }: { mixtureId: string; componentId: string }) => 
      mixtureService.removeComponent(mixtureId, componentId),
    onSuccess: (_, { mixtureId }) => {
      // Invalidate specific mixture query to refetch after update
      queryClient.invalidateQueries({ queryKey: [MIXTURES_QUERY_KEY, mixtureId] });
    },
  });
}

export function useCryoprotectionScore(mixtureId: string) {
  return useQuery({
    queryKey: [MIXTURES_QUERY_KEY, mixtureId, 'cryoprotection-score'],
    queryFn: () => mixtureService.getCryoprotectionScore(mixtureId),
    staleTime: 5 * 60 * 1000, // 5 minutes
    enabled: !!mixtureId,
  });
}

export function useSearchMixtures(query: string, limit: number = 10) {
  return useQuery({
    queryKey: [MIXTURES_QUERY_KEY, 'search', query, limit],
    queryFn: () => mixtureService.searchMixtures(query, limit),
    staleTime: 5 * 60 * 1000, // 5 minutes
    enabled: query.length > 2, // Only search when query has 3+ characters
  });
}
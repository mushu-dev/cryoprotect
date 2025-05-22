/**
 * Enhanced Molecules Hooks
 * Provides React Query hooks for enhanced molecular data with real-time updates
 */
import React from 'react';
import { useQuery, useMutation, useQueryClient } from '@tanstack/react-query';
import { enhancedConvexMoleculeService } from '../services/enhanced-convex-service';

export const ENHANCED_MOLECULES_QUERY_KEY = 'enhanced_molecules';

/**
 * Hook to get a list of molecules with enhanced properties and scores
 * @param {Object} params - Query parameters
 * @returns {Object} React Query result with enhanced data
 */
export function useEnhancedMolecules(params = {}) {
  return useQuery({
    queryKey: [ENHANCED_MOLECULES_QUERY_KEY, params],
    queryFn: () => enhancedConvexMoleculeService.getMolecules(params),
    staleTime: 2 * 60 * 1000, // 2 minutes (shorter for real-time feel)
    cacheTime: 10 * 60 * 1000, // 10 minutes
    keepPreviousData: true, // Keep previous data while loading new
    retry: 2,
    retryDelay: attemptIndex => Math.min(1000 * 2 ** attemptIndex, 30000),
  });
}

/**
 * Hook to get a single molecule with all enhanced data
 * @param {string} id - Molecule ID
 * @returns {Object} React Query result
 */
export function useEnhancedMolecule(id) {
  return useQuery({
    queryKey: [ENHANCED_MOLECULES_QUERY_KEY, id],
    queryFn: () => enhancedConvexMoleculeService.getMolecule(id),
    staleTime: 5 * 60 * 1000, // 5 minutes
    cacheTime: 15 * 60 * 1000, // 15 minutes
    enabled: !!id,
    retry: 2,
  });
}

/**
 * Hook to get molecule properties with enhanced data
 * @param {string} id - Molecule ID
 * @returns {Object} React Query result
 */
export function useEnhancedMoleculeProperties(id) {
  return useQuery({
    queryKey: [ENHANCED_MOLECULES_QUERY_KEY, id, 'properties'],
    queryFn: () => enhancedConvexMoleculeService.getMoleculeProperties(id),
    staleTime: 5 * 60 * 1000, // 5 minutes
    cacheTime: 15 * 60 * 1000, // 15 minutes
    enabled: !!id,
    retry: 2,
  });
}

/**
 * Hook to search molecules with enhanced scoring
 * @param {string} query - Search query
 * @param {number} limit - Maximum number of results
 * @returns {Object} React Query result
 */
export function useEnhancedSearchMolecules(query, limit = 10) {
  return useQuery({
    queryKey: [ENHANCED_MOLECULES_QUERY_KEY, 'search', query, limit],
    queryFn: () => enhancedConvexMoleculeService.searchMolecules(query, limit),
    staleTime: 3 * 60 * 1000, // 3 minutes
    cacheTime: 10 * 60 * 1000, // 10 minutes
    enabled: query.length > 2, // Only search when query has 3+ characters
    retry: 1,
  });
}

/**
 * Hook to get cryoprotectant ranking
 * @param {string} category - Category filter (optional)
 * @param {number} limit - Maximum number of results
 * @returns {Object} React Query result
 */
export function useCryoprotectantRanking(category = null, limit = 10) {
  return useQuery({
    queryKey: [ENHANCED_MOLECULES_QUERY_KEY, 'ranking', category, limit],
    queryFn: () => enhancedConvexMoleculeService.getCryoprotectantRanking(category, limit),
    staleTime: 10 * 60 * 1000, // 10 minutes
    cacheTime: 30 * 60 * 1000, // 30 minutes
    retry: 2,
  });
}

/**
 * Hook to get database quality metrics
 * @returns {Object} React Query result
 */
export function useQualityMetrics() {
  return useQuery({
    queryKey: [ENHANCED_MOLECULES_QUERY_KEY, 'quality_metrics'],
    queryFn: () => enhancedConvexMoleculeService.getQualityMetrics(),
    staleTime: 15 * 60 * 1000, // 15 minutes
    cacheTime: 60 * 60 * 1000, // 1 hour
    retry: 1,
  });
}

/**
 * Hook to get client performance metrics
 * @returns {Object} Current performance metrics
 */
export function useClientMetrics() {
  return useQuery({
    queryKey: [ENHANCED_MOLECULES_QUERY_KEY, 'client_metrics'],
    queryFn: () => enhancedConvexMoleculeService.getClientMetrics(),
    staleTime: 30 * 1000, // 30 seconds
    cacheTime: 2 * 60 * 1000, // 2 minutes
    refetchInterval: 30 * 1000, // Refresh every 30 seconds for live metrics
    retry: false, // Don't retry metrics queries
  });
}

/**
 * Hook to invalidate and refresh all molecular data
 * @returns {Function} Function to trigger refresh
 */
export function useRefreshMolecules() {
  const queryClient = useQueryClient();
  
  return () => {
    // Clear service cache
    enhancedConvexMoleculeService.clearCache();
    
    // Invalidate all related queries
    queryClient.invalidateQueries({ 
      queryKey: [ENHANCED_MOLECULES_QUERY_KEY] 
    });
    
    console.log('🔄 Enhanced molecular data refreshed');
  };
}

/**
 * Custom hook for real-time molecule updates
 * @param {Array} moleculeIds - Array of molecule IDs to watch
 * @returns {Object} Real-time update status
 */
export function useRealTimeMolecules(moleculeIds = []) {
  const queryClient = useQueryClient();
  
  return useQuery({
    queryKey: [ENHANCED_MOLECULES_QUERY_KEY, 'realtime', moleculeIds],
    queryFn: async () => {
      // Get current data for all molecules
      const molecules = await Promise.all(
        moleculeIds.map(id => 
          enhancedConvexMoleculeService.getMolecule(id).catch(() => null)
        )
      );
      
      return {
        molecules: molecules.filter(mol => mol !== null),
        lastUpdate: new Date().toISOString(),
        count: molecules.filter(mol => mol !== null).length
      };
    },
    staleTime: 0, // Always fresh for real-time
    cacheTime: 1 * 60 * 1000, // 1 minute cache
    refetchInterval: 30 * 1000, // Refresh every 30 seconds
    enabled: moleculeIds.length > 0,
    retry: 1,
    onSuccess: (data) => {
      // Update individual molecule caches with fresh data
      data.molecules.forEach(molecule => {
        queryClient.setQueryData(
          [ENHANCED_MOLECULES_QUERY_KEY, molecule.id], 
          molecule
        );
      });
    }
  });
}

/**
 * Hook for batch operations on molecules
 * @returns {Object} Mutation functions for batch operations
 */
export function useBatchMoleculeOperations() {
  const queryClient = useQueryClient();
  
  const refreshBatch = useMutation({
    mutationFn: async (moleculeIds) => {
      // Refresh multiple molecules at once
      const results = await Promise.allSettled(
        moleculeIds.map(id => enhancedConvexMoleculeService.getMolecule(id))
      );
      
      return {
        successful: results.filter(r => r.status === 'fulfilled').length,
        failed: results.filter(r => r.status === 'rejected').length,
        total: moleculeIds.length
      };
    },
    onSuccess: (data, moleculeIds) => {
      // Invalidate queries for updated molecules
      moleculeIds.forEach(id => {
        queryClient.invalidateQueries({ 
          queryKey: [ENHANCED_MOLECULES_QUERY_KEY, id] 
        });
      });
      
      console.log(`🔄 Batch refresh: ${data.successful}/${data.total} successful`);
    }
  });
  
  return {
    refreshBatch
  };
}

/**
 * Hook for advanced filtering and sorting
 * @param {Object} molecules - Array of molecules to filter/sort
 * @param {Object} filters - Filter configuration
 * @returns {Array} Filtered and sorted molecules
 */
export function useAdvancedMoleculeFiltering(molecules = [], filters = {}) {
  return React.useMemo(() => {
    if (!molecules || molecules.length === 0) return [];
    
    let filtered = [...molecules];
    
    // Apply molecular weight filter
    if (filters.molecularWeightRange) {
      const [min, max] = filters.molecularWeightRange;
      filtered = filtered.filter(mol => 
        mol.molecularWeight >= min && mol.molecularWeight <= max
      );
    }
    
    // Apply LogP filter
    if (filters.logPRange) {
      const [min, max] = filters.logPRange;
      filtered = filtered.filter(mol => 
        mol.logP >= min && mol.logP <= max
      );
    }
    
    // Apply cryoprotectant score filter
    if (filters.minCryoprotectantScore) {
      filtered = filtered.filter(mol => 
        mol.cryoprotectantScore >= filters.minCryoprotectantScore
      );
    }
    
    // Apply category filter
    if (filters.category) {
      filtered = filtered.filter(mol => 
        mol.cryoprotectantCategory === filters.category
      );
    }
    
    // Apply advanced sorting
    if (filters.sortBy) {
      filtered.sort((a, b) => {
        const aValue = a[filters.sortBy];
        const bValue = b[filters.sortBy];
        
        if (aValue == null && bValue == null) return 0;
        if (aValue == null) return 1;
        if (bValue == null) return -1;
        
        const direction = filters.sortOrder === 'desc' ? -1 : 1;
        
        if (typeof aValue === 'string') {
          return direction * aValue.localeCompare(bValue);
        } else {
          return direction * (aValue - bValue);
        }
      });
    }
    
    return filtered;
  }, [molecules, filters]);
}
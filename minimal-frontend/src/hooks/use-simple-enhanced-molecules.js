/**
 * Simple Enhanced Molecules Hooks
 * Provides hooks for enhanced molecular data without React Query dependency
 */
import { useState, useEffect, useCallback } from 'react';
import { enhancedConvexMoleculeService } from '../services/enhanced-convex-service';

/**
 * Simple hook to get a list of molecules with enhanced properties
 */
export function useSimpleEnhancedMolecules(params = {}) {
  const [data, setData] = useState(null);
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState(null);

  const fetchMolecules = useCallback(async () => {
    setIsLoading(true);
    setError(null);
    
    try {
      const result = await enhancedConvexMoleculeService.getMolecules(params);
      setData(result);
    } catch (err) {
      setError(err);
      console.error('Error fetching enhanced molecules:', err);
    } finally {
      setIsLoading(false);
    }
  }, [params]);

  useEffect(() => {
    fetchMolecules();
  }, [fetchMolecules]);

  return {
    data,
    isLoading,
    error,
    refetch: fetchMolecules
  };
}

/**
 * Simple hook to get a single molecule with enhanced data
 */
export function useSimpleEnhancedMolecule(id) {
  const [data, setData] = useState(null);
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState(null);

  const fetchMolecule = useCallback(async () => {
    if (!id) return;
    
    setIsLoading(true);
    setError(null);
    
    try {
      const result = await enhancedConvexMoleculeService.getMolecule(id);
      setData(result);
    } catch (err) {
      setError(err);
      console.error(`Error fetching enhanced molecule ${id}:`, err);
    } finally {
      setIsLoading(false);
    }
  }, [id]);

  useEffect(() => {
    fetchMolecule();
  }, [fetchMolecule]);

  return {
    data,
    isLoading,
    error,
    refetch: fetchMolecule
  };
}

/**
 * Simple hook to get quality metrics
 */
export function useSimpleQualityMetrics() {
  const [data, setData] = useState(null);
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState(null);

  const fetchMetrics = useCallback(async () => {
    setIsLoading(true);
    setError(null);
    
    try {
      const result = await enhancedConvexMoleculeService.getQualityMetrics();
      setData(result);
    } catch (err) {
      setError(err);
      console.error('Error fetching quality metrics:', err);
    } finally {
      setIsLoading(false);
    }
  }, []);

  useEffect(() => {
    fetchMetrics();
  }, [fetchMetrics]);

  return {
    data,
    isLoading,
    error,
    refetch: fetchMetrics
  };
}

/**
 * Simple hook to get client performance metrics
 */
export function useSimpleClientMetrics() {
  const [data, setData] = useState(null);

  useEffect(() => {
    const updateMetrics = () => {
      const metrics = enhancedConvexMoleculeService.getClientMetrics();
      setData(metrics);
    };

    // Update immediately
    updateMetrics();

    // Update every 30 seconds
    const interval = setInterval(updateMetrics, 30000);

    return () => clearInterval(interval);
  }, []);

  return { data };
}

/**
 * Simple hook to refresh all molecular data
 */
export function useSimpleRefreshMolecules() {
  return useCallback(() => {
    enhancedConvexMoleculeService.clearCache();
    console.log('🔄 Enhanced molecular data cache cleared');
  }, []);
}
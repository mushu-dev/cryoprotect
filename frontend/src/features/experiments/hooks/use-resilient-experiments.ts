/**
 * React hooks for working with experiments using the resilient API service
 */

import { useState, useEffect, useCallback } from 'react';
import { UUID } from 'crypto';
import { 
  createResilienceEnhancedExperimentService,
  ResilienceEnhancedExperimentService 
} from '../services/enhanced-experiment-service';
import { MockExperimentService } from '../services/experiment-service';
import { MOCK_EXPERIMENTS } from '../data/mock-experiments';

// Create a mock service with initial data
const createMockService = () => {
  const mockService = new MockExperimentService();
  
  // Pre-populate with mock data
  MOCK_EXPERIMENTS.forEach(experiment => {
    // @ts-ignore - Handle type incompatibility for demo
    mockService['experiments'].push(experiment);
  });
  
  return mockService;
};

// Create the service with fallback to mock data
const experimentService = createResilienceEnhancedExperimentService(createMockService());

/**
 * Hook for listing experiments with filtering
 */
export const useExperimentsList = (params?: {
  page?: number;
  perPage?: number;
  status?: string;
  experimentType?: string;
  researcher?: string;
  tissueTypeId?: UUID;
  dateFrom?: string;
  dateTo?: string;
  tags?: string[];
  search?: string;
  sortBy?: string;
  sortOrder?: 'asc' | 'desc';
}) => {
  const [experiments, setExperiments] = useState<any[]>([]);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<Error | null>(null);
  const [pagination, setPagination] = useState({
    total: 0,
    page: params?.page || 1,
    perPage: params?.perPage || 10,
    totalPages: 0
  });

  const fetchExperiments = useCallback(async () => {
    setLoading(true);
    setError(null);
    
    try {
      // Convert params to API format
      const apiParams = params ? {
        page: params.page,
        per_page: params.perPage,
        status: params.status,
        experiment_type: params.experimentType,
        researcher: params.researcher,
        tissue_type_id: params.tissueTypeId,
        date_from: params.dateFrom,
        date_to: params.dateTo,
        tags: params.tags,
        search: params.search,
        sort_by: params.sortBy,
        sort_order: params.sortOrder
      } : undefined;
      
      const response = await experimentService.getExperiments(apiParams);
      
      setExperiments(response.data);
      setPagination({
        total: response.total,
        page: response.page,
        perPage: response.per_page,
        totalPages: response.total_pages
      });
    } catch (err) {
      setError(err instanceof Error ? err : new Error('Failed to fetch experiments'));
      console.error('Error fetching experiments:', err);
    } finally {
      setLoading(false);
    }
  }, [params]);

  useEffect(() => {
    fetchExperiments();
  }, [fetchExperiments]);

  const refetch = () => {
    fetchExperiments();
  };

  return { experiments, loading, error, pagination, refetch };
};

/**
 * Hook for getting a single experiment by ID
 */
export const useExperiment = (id?: UUID) => {
  const [experiment, setExperiment] = useState<any | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<Error | null>(null);

  const fetchExperiment = useCallback(async () => {
    if (!id) return;
    
    setLoading(true);
    setError(null);
    
    try {
      const data = await experimentService.getExperiment(id);
      setExperiment(data);
    } catch (err) {
      setError(err instanceof Error ? err : new Error(`Failed to fetch experiment ${id}`));
      console.error(`Error fetching experiment ${id}:`, err);
    } finally {
      setLoading(false);
    }
  }, [id]);

  useEffect(() => {
    fetchExperiment();
  }, [fetchExperiment]);

  const refetch = () => {
    fetchExperiment();
  };

  return { experiment, loading, error, refetch };
};

/**
 * Hook for creating a new experiment
 */
export const useCreateExperiment = () => {
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<Error | null>(null);
  const [createdExperiment, setCreatedExperiment] = useState<any | null>(null);

  const createExperiment = async (experimentData: any) => {
    setLoading(true);
    setError(null);
    
    try {
      const data = await experimentService.createExperiment(experimentData);
      setCreatedExperiment(data);
      return data;
    } catch (err) {
      setError(err instanceof Error ? err : new Error('Failed to create experiment'));
      console.error('Error creating experiment:', err);
      throw err;
    } finally {
      setLoading(false);
    }
  };

  return { createExperiment, loading, error, createdExperiment };
};

/**
 * Hook for updating an experiment
 */
export const useUpdateExperiment = (id?: UUID) => {
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<Error | null>(null);
  const [updatedExperiment, setUpdatedExperiment] = useState<any | null>(null);

  const updateExperiment = async (experimentData: any) => {
    if (!id) {
      const error = new Error('Experiment ID is required');
      setError(error);
      throw error;
    }
    
    setLoading(true);
    setError(null);
    
    try {
      const data = await experimentService.updateExperiment(id, experimentData);
      setUpdatedExperiment(data);
      return data;
    } catch (err) {
      setError(err instanceof Error ? err : new Error(`Failed to update experiment ${id}`));
      console.error(`Error updating experiment ${id}:`, err);
      throw err;
    } finally {
      setLoading(false);
    }
  };

  return { updateExperiment, loading, error, updatedExperiment };
};

/**
 * Hook for deleting an experiment
 */
export const useDeleteExperiment = () => {
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<Error | null>(null);
  const [deleted, setDeleted] = useState(false);

  const deleteExperiment = async (id: UUID) => {
    setLoading(true);
    setError(null);
    setDeleted(false);
    
    try {
      await experimentService.deleteExperiment(id);
      setDeleted(true);
    } catch (err) {
      setError(err instanceof Error ? err : new Error(`Failed to delete experiment ${id}`));
      console.error(`Error deleting experiment ${id}:`, err);
      throw err;
    } finally {
      setLoading(false);
    }
  };

  return { deleteExperiment, loading, error, deleted };
};

/**
 * Hook for getting experiment results
 */
export const useExperimentResults = (experimentId?: UUID) => {
  const [results, setResults] = useState<any[]>([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<Error | null>(null);

  const fetchResults = useCallback(async () => {
    if (!experimentId) return;
    
    setLoading(true);
    setError(null);
    
    try {
      const data = await experimentService.getExperimentResults(experimentId);
      setResults(data);
    } catch (err) {
      setError(err instanceof Error ? err : new Error(`Failed to fetch results for experiment ${experimentId}`));
      console.error(`Error fetching results for experiment ${experimentId}:`, err);
    } finally {
      setLoading(false);
    }
  }, [experimentId]);

  useEffect(() => {
    fetchResults();
  }, [fetchResults]);

  const refetch = () => {
    fetchResults();
  };

  return { results, loading, error, refetch };
};

/**
 * Hook for adding a result to an experiment
 */
export const useAddExperimentResult = (experimentId?: UUID) => {
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<Error | null>(null);
  const [addedResult, setAddedResult] = useState<any | null>(null);

  const addResult = async (resultData: any) => {
    if (!experimentId) {
      const error = new Error('Experiment ID is required');
      setError(error);
      throw error;
    }
    
    setLoading(true);
    setError(null);
    
    try {
      const data = await experimentService.addExperimentResult({
        ...resultData,
        experiment_id: experimentId
      });
      setAddedResult(data);
      return data;
    } catch (err) {
      setError(err instanceof Error ? err : new Error('Failed to add experiment result'));
      console.error('Error adding experiment result:', err);
      throw err;
    } finally {
      setLoading(false);
    }
  };

  return { addResult, loading, error, addedResult };
};

/**
 * Hook for analyzing experiments
 */
export const useExperimentAnalysis = () => {
  const [analysisResult, setAnalysisResult] = useState<any | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<Error | null>(null);

  const analyzeExperiments = async (
    experimentIds: UUID[],
    parameters?: {
      compareWith?: UUID[];
      analysisTypes?: string[];
      timeRange?: [string, string];
    }
  ) => {
    setLoading(true);
    setError(null);
    
    try {
      // Convert parameters to API format
      const apiParams = parameters ? {
        compare_with: parameters.compareWith,
        analysis_type: parameters.analysisTypes,
        time_range: parameters.timeRange
      } : undefined;
      
      const data = await experimentService.analyzeExperiments(experimentIds, apiParams);
      setAnalysisResult(data);
      return data;
    } catch (err) {
      setError(err instanceof Error ? err : new Error('Failed to analyze experiments'));
      console.error('Error analyzing experiments:', err);
      throw err;
    } finally {
      setLoading(false);
    }
  };

  return { analyzeExperiments, loading, error, analysisResult };
};

/**
 * Hook for experiment comparison functionality
 */
export const useExperimentComparison = () => {
  const [comparisonData, setComparisonData] = useState<{
    experiments: any[];
    metrics: string[];
    data: Record<string, any>;
  } | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<Error | null>(null);

  const compareExperiments = async (experimentIds: UUID[]) => {
    setLoading(true);
    setError(null);
    
    try {
      // Fetch all experiments in parallel
      const experiments = await Promise.all(
        experimentIds.map(id => experimentService.getExperiment(id))
      );
      
      // Fetch results for each experiment
      const experimentsWithResults = await Promise.all(
        experiments.map(async experiment => {
          const results = await experimentService.getExperimentResults(experiment.id);
          return {
            ...experiment,
            results
          };
        })
      );
      
      // Extract common metrics from all experiments
      const metricsSet = new Set<string>();
      experimentsWithResults.forEach(experiment => {
        experiment.results.forEach((result: any) => {
          if (result.viability_percentage !== undefined) metricsSet.add('viability_percentage');
          if (result.recovery_rate !== undefined) metricsSet.add('recovery_rate');
          if (result.functionality_score !== undefined) metricsSet.add('functionality_score');
          
          // Add any custom metrics from result_details
          if (result.result_details) {
            Object.keys(result.result_details).forEach(key => metricsSet.add(key));
          }
        });
      });
      
      const metrics = Array.from(metricsSet);
      
      // Organize comparison data by metric
      const data: Record<string, any> = {};
      
      metrics.forEach(metric => {
        data[metric] = {
          labels: experimentsWithResults.map(exp => exp.name),
          datasets: [
            {
              label: metric,
              data: experimentsWithResults.map(exp => {
                // Calculate average value for this metric across all results
                const values = exp.results
                  .map((result: any) => {
                    if (metric === 'viability_percentage') return result.viability_percentage;
                    if (metric === 'recovery_rate') return result.recovery_rate;
                    if (metric === 'functionality_score') return result.functionality_score;
                    return result.result_details?.[metric];
                  })
                  .filter((value: any) => value !== undefined);
                
                if (values.length === 0) return null;
                
                // Calculate average
                const sum = values.reduce((a: number, b: number) => a + b, 0);
                return sum / values.length;
              }),
            }
          ]
        };
      });
      
      setComparisonData({
        experiments: experimentsWithResults,
        metrics,
        data
      });
      
      return {
        experiments: experimentsWithResults,
        metrics,
        data
      };
    } catch (err) {
      setError(err instanceof Error ? err : new Error('Failed to compare experiments'));
      console.error('Error comparing experiments:', err);
      throw err;
    } finally {
      setLoading(false);
    }
  };

  return { compareExperiments, loading, error, comparisonData };
};
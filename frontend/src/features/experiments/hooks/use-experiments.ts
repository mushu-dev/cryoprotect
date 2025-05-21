'use client'

import { useState, useEffect, useCallback } from 'react'
import { 
  Experiment, 
  ExperimentService, 
  ExperimentListParams, 
  ExperimentServiceImpl,
  MockExperimentService,
  ExperimentResult
} from '../services/experiment-service'
import { UUID } from 'crypto'

const useMockData = process.env.NEXT_PUBLIC_USE_MOCK_DATA === 'true'

export function useExperiments(initialParams?: ExperimentListParams) {
  const [experiments, setExperiments] = useState<Experiment[]>([])
  const [loading, setLoading] = useState(true)
  const [error, setError] = useState<Error | null>(null)
  const [totalCount, setTotalCount] = useState(0)
  const [page, setPage] = useState(initialParams?.page || 1)
  const [perPage, setPerPage] = useState(initialParams?.per_page || 10)
  const [params, setParams] = useState<ExperimentListParams>(initialParams || {})
  
  const experimentService: ExperimentService = useMockData 
    ? new MockExperimentService()
    : new ExperimentServiceImpl()
  
  const fetchExperiments = useCallback(async () => {
    try {
      setLoading(true)
      setError(null)
      
      const queryParams: ExperimentListParams = {
        ...params,
        page,
        per_page: perPage
      }
      
      const response = await experimentService.getExperiments(queryParams)
      
      setExperiments(response.data)
      setTotalCount(response.total)
    } catch (err) {
      setError(err instanceof Error ? err : new Error('Failed to fetch experiments'))
    } finally {
      setLoading(false)
    }
  }, [experimentService, params, page, perPage])
  
  useEffect(() => {
    fetchExperiments()
  }, [fetchExperiments])
  
  const updateParams = useCallback((newParams: Partial<ExperimentListParams>) => {
    setParams(prevParams => ({
      ...prevParams,
      ...newParams
    }))
    
    // Reset to first page when changing filters
    if (!('page' in newParams)) {
      setPage(1)
    }
  }, [])
  
  const changePage = useCallback((newPage: number) => {
    setPage(newPage)
  }, [])
  
  const changePerPage = useCallback((newPerPage: number) => {
    setPerPage(newPerPage)
    setPage(1) // Reset to first page when changing items per page
  }, [])
  
  const refreshExperiments = useCallback(() => {
    fetchExperiments()
  }, [fetchExperiments])
  
  return {
    experiments,
    loading,
    error,
    totalCount,
    page,
    perPage,
    totalPages: Math.ceil(totalCount / perPage),
    updateParams,
    changePage,
    changePerPage,
    refreshExperiments
  }
}

export function useExperiment(id: UUID | null) {
  const [experiment, setExperiment] = useState<Experiment | null>(null)
  const [loading, setLoading] = useState(false)
  const [error, setError] = useState<Error | null>(null)
  
  const experimentService: ExperimentService = useMockData 
    ? new MockExperimentService()
    : new ExperimentServiceImpl()
  
  const fetchExperiment = useCallback(async () => {
    if (!id) return
    
    try {
      setLoading(true)
      setError(null)
      
      const response = await experimentService.getExperiment(id)
      setExperiment(response)
    } catch (err) {
      setError(err instanceof Error ? err : new Error('Failed to fetch experiment'))
    } finally {
      setLoading(false)
    }
  }, [experimentService, id])
  
  useEffect(() => {
    fetchExperiment()
  }, [fetchExperiment])
  
  const updateExperiment = useCallback(async (updates: Partial<Experiment>) => {
    if (!id || !experiment) return null
    
    try {
      setLoading(true)
      setError(null)
      
      const updatedExperiment = await experimentService.updateExperiment(id, updates)
      setExperiment(updatedExperiment)
      return updatedExperiment
    } catch (err) {
      setError(err instanceof Error ? err : new Error('Failed to update experiment'))
      return null
    } finally {
      setLoading(false)
    }
  }, [experimentService, id, experiment])
  
  const fetchResults = useCallback(async () => {
    if (!id) return []
    
    try {
      setLoading(true)
      setError(null)
      
      const results = await experimentService.getExperimentResults(id)
      
      // Update the experiment with the results
      if (experiment) {
        setExperiment({
          ...experiment,
          results
        })
      }
      
      return results
    } catch (err) {
      setError(err instanceof Error ? err : new Error('Failed to fetch experiment results'))
      return []
    } finally {
      setLoading(false)
    }
  }, [experimentService, id, experiment])
  
  const addResult = useCallback(async (result: Omit<ExperimentResult, 'id' | 'provenance'>) => {
    if (!id) return null
    
    try {
      setLoading(true)
      setError(null)
      
      const newResult = await experimentService.addExperimentResult(result)
      
      // Update the experiment with the new result
      if (experiment) {
        const updatedResults = experiment.results ? [...experiment.results, newResult] : [newResult]
        setExperiment({
          ...experiment,
          results: updatedResults
        })
      }
      
      return newResult
    } catch (err) {
      setError(err instanceof Error ? err : new Error('Failed to add experiment result'))
      return null
    } finally {
      setLoading(false)
    }
  }, [experimentService, id, experiment])
  
  const refreshExperiment = useCallback(() => {
    fetchExperiment()
  }, [fetchExperiment])
  
  return {
    experiment,
    loading,
    error,
    updateExperiment,
    refreshExperiment,
    fetchResults,
    addResult
  }
}

export function useExperimentAnalysis(experimentIds: UUID[]) {
  const [analysisResult, setAnalysisResult] = useState<any | null>(null)
  const [loading, setLoading] = useState(false)
  const [error, setError] = useState<Error | null>(null)
  
  const experimentService: ExperimentService = useMockData 
    ? new MockExperimentService()
    : new ExperimentServiceImpl()
  
  const runAnalysis = useCallback(async (parameters?: {
    compare_with?: UUID[];
    analysis_type?: string[];
    time_range?: [string, string];
  }) => {
    if (experimentIds.length === 0) return
    
    try {
      setLoading(true)
      setError(null)
      
      const result = await experimentService.analyzeExperiments(experimentIds, parameters)
      setAnalysisResult(result)
      return result
    } catch (err) {
      setError(err instanceof Error ? err : new Error('Failed to analyze experiments'))
      return null
    } finally {
      setLoading(false)
    }
  }, [experimentService, experimentIds])
  
  return {
    analysisResult,
    loading,
    error,
    runAnalysis
  }
}

export function useExperimentCreation() {
  const [creating, setCreating] = useState(false)
  const [error, setError] = useState<Error | null>(null)
  
  const experimentService: ExperimentService = useMockData 
    ? new MockExperimentService()
    : new ExperimentServiceImpl()
  
  const createExperiment = useCallback(async (experiment: Omit<Experiment, 'id' | 'provenance'>) => {
    try {
      setCreating(true)
      setError(null)
      
      const newExperiment = await experimentService.createExperiment(experiment)
      return newExperiment
    } catch (err) {
      setError(err instanceof Error ? err : new Error('Failed to create experiment'))
      return null
    } finally {
      setCreating(false)
    }
  }, [experimentService])
  
  const importExperiment = useCallback(async (file: File) => {
    try {
      setCreating(true)
      setError(null)
      
      const importedExperiment = await experimentService.importExperiment(file)
      return importedExperiment
    } catch (err) {
      setError(err instanceof Error ? err : new Error('Failed to import experiment'))
      return null
    } finally {
      setCreating(false)
    }
  }, [experimentService])
  
  return {
    creating,
    error,
    createExperiment,
    importExperiment
  }
}

export function useExperimentSearch() {
  const [results, setResults] = useState<Experiment[]>([])
  const [loading, setLoading] = useState(false)
  const [error, setError] = useState<Error | null>(null)
  const [totalCount, setTotalCount] = useState(0)
  
  const experimentService: ExperimentService = useMockData 
    ? new MockExperimentService()
    : new ExperimentServiceImpl()
  
  const searchExperiments = useCallback(async (query: string, params?: ExperimentListParams) => {
    try {
      setLoading(true)
      setError(null)
      
      const response = await experimentService.searchExperiments(query, params)
      
      setResults(response.data)
      setTotalCount(response.total)
      return response
    } catch (err) {
      setError(err instanceof Error ? err : new Error('Failed to search experiments'))
      return { data: [], total: 0 }
    } finally {
      setLoading(false)
    }
  }, [experimentService])
  
  return {
    results,
    loading,
    error,
    totalCount,
    searchExperiments
  }
}
'use client'

import { useState, useEffect, useCallback } from 'react'
import { 
  Protocol, 
  ProtocolService,
  ProtocolListParams,
  ValidationResult,
  ProtocolDiff,
  ProtocolLibraryItem
} from '../services/protocol-service'
import { ApiProtocolService } from '../services/protocol-service-api'
import { getProtocolService } from '@/features/common/api/service-factory'
import { UUID } from 'crypto'

// Now using service factory instead of direct instantiation

export function useProtocols(initialParams?: ProtocolListParams) {
  const [protocols, setProtocols] = useState<Protocol[]>([])
  const [loading, setLoading] = useState(true)
  const [error, setError] = useState<Error | null>(null)
  const [totalCount, setTotalCount] = useState(0)
  const [page, setPage] = useState(initialParams?.page || 1)
  const [perPage, setPerPage] = useState(initialParams?.per_page || 10)
  const [params, setParams] = useState<ProtocolListParams>(initialParams || {})
  
  const protocolService: ProtocolService = getProtocolService()
  
  const fetchProtocols = useCallback(async () => {
    try {
      setLoading(true)
      setError(null)
      
      const queryParams: ProtocolListParams = {
        ...params,
        page,
        per_page: perPage
      }
      
      const response = await protocolService.getProtocols(queryParams)
      
      setProtocols(response.data)
      setTotalCount(response.total)
    } catch (err) {
      setError(err instanceof Error ? err : new Error('Failed to fetch protocols'))
    } finally {
      setLoading(false)
    }
  }, [protocolService, params, page, perPage])
  
  useEffect(() => {
    fetchProtocols()
  }, [fetchProtocols])
  
  const updateParams = useCallback((newParams: Partial<ProtocolListParams>) => {
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
  
  const refreshProtocols = useCallback(() => {
    fetchProtocols()
  }, [fetchProtocols])
  
  return {
    protocols,
    loading,
    error,
    totalCount,
    page,
    perPage,
    totalPages: Math.ceil(totalCount / perPage),
    updateParams,
    changePage,
    changePerPage,
    refreshProtocols
  }
}

export function useProtocol(id: UUID | null) {
  const [protocol, setProtocol] = useState<Protocol | null>(null)
  const [loading, setLoading] = useState(false)
  const [error, setError] = useState<Error | null>(null)
  
  const protocolService: ProtocolService = getProtocolService()
  
  const fetchProtocol = useCallback(async () => {
    if (!id) return
    
    try {
      setLoading(true)
      setError(null)
      
      const response = await protocolService.getProtocol(id)
      setProtocol(response)
    } catch (err) {
      setError(err instanceof Error ? err : new Error('Failed to fetch protocol'))
    } finally {
      setLoading(false)
    }
  }, [protocolService, id])
  
  useEffect(() => {
    fetchProtocol()
  }, [fetchProtocol])
  
  const updateProtocol = useCallback(async (updates: Partial<Protocol>) => {
    if (!id || !protocol) return null
    
    try {
      setLoading(true)
      setError(null)
      
      const updatedProtocol = await protocolService.updateProtocol(id, updates)
      setProtocol(updatedProtocol)
      return updatedProtocol
    } catch (err) {
      setError(err instanceof Error ? err : new Error('Failed to update protocol'))
      return null
    } finally {
      setLoading(false)
    }
  }, [protocolService, id, protocol])
  
  const validateProtocol = useCallback(async (): Promise<ValidationResult | null> => {
    if (!protocol) return null
    
    try {
      setLoading(true)
      setError(null)
      
      const validationResult = await protocolService.validateProtocol(protocol)
      return validationResult
    } catch (err) {
      setError(err instanceof Error ? err : new Error('Failed to validate protocol'))
      return null
    } finally {
      setLoading(false)
    }
  }, [protocolService, protocol])
  
  const createVersion = useCallback(async (
    changes: Partial<Protocol>,
    versionNotes?: string
  ) => {
    if (!id) return null
    
    try {
      setLoading(true)
      setError(null)
      
      const newVersion = await protocolService.createProtocolVersion(id, changes, versionNotes)
      return newVersion
    } catch (err) {
      setError(err instanceof Error ? err : new Error('Failed to create protocol version'))
      return null
    } finally {
      setLoading(false)
    }
  }, [protocolService, id])
  
  const exportProtocol = useCallback(async (
    format: 'json' | 'yaml' | 'pdf' | 'human-readable'
  ) => {
    if (!id) return null
    
    try {
      setLoading(true)
      setError(null)
      
      const blob = await protocolService.exportProtocol(id, format)
      
      // Create a download link for the blob
      const url = URL.createObjectURL(blob)
      const a = document.createElement('a')
      a.href = url
      a.download = `protocol-${id}.${format === 'human-readable' ? 'txt' : format}`
      document.body.appendChild(a)
      a.click()
      a.remove()
      URL.revokeObjectURL(url)
      
      return true
    } catch (err) {
      setError(err instanceof Error ? err : new Error('Failed to export protocol'))
      return null
    } finally {
      setLoading(false)
    }
  }, [protocolService, id])
  
  const publishToLibrary = useCallback(async (metadata: {
    category: string;
    tags: string[];
    description?: string;
  }) => {
    if (!id) return null
    
    try {
      setLoading(true)
      setError(null)
      
      const libraryItem = await protocolService.publishToLibrary(id, metadata)
      return libraryItem
    } catch (err) {
      setError(err instanceof Error ? err : new Error('Failed to publish protocol to library'))
      return null
    } finally {
      setLoading(false)
    }
  }, [protocolService, id])
  
  const refreshProtocol = useCallback(() => {
    fetchProtocol()
  }, [fetchProtocol])
  
  return {
    protocol,
    loading,
    error,
    updateProtocol,
    validateProtocol,
    createVersion,
    exportProtocol,
    publishToLibrary,
    refreshProtocol
  }
}

export function useProtocolVersions(id: UUID | null) {
  const [versions, setVersions] = useState<{
    version: string;
    id: UUID;
    created_at: string;
    created_by: string;
    description?: string;
  }[]>([])
  const [loading, setLoading] = useState(false)
  const [error, setError] = useState<Error | null>(null)
  
  const protocolService: ProtocolService = getProtocolService()
  
  const fetchVersions = useCallback(async () => {
    if (!id) return
    
    try {
      setLoading(true)
      setError(null)
      
      const response = await protocolService.getProtocolVersions(id)
      setVersions(response.versions)
    } catch (err) {
      setError(err instanceof Error ? err : new Error('Failed to fetch protocol versions'))
    } finally {
      setLoading(false)
    }
  }, [protocolService, id])
  
  useEffect(() => {
    fetchVersions()
  }, [fetchVersions])
  
  const compareVersions = useCallback(async (
    versionId1: UUID,
    versionId2: UUID
  ): Promise<ProtocolDiff | null> => {
    try {
      setLoading(true)
      setError(null)
      
      const diff = await protocolService.compareProtocolVersions(versionId1, versionId2)
      return diff
    } catch (err) {
      setError(err instanceof Error ? err : new Error('Failed to compare protocol versions'))
      return null
    } finally {
      setLoading(false)
    }
  }, [protocolService])
  
  return {
    versions,
    loading,
    error,
    compareVersions,
    refreshVersions: fetchVersions
  }
}

export function useProtocolTemplates(initialParams?: ProtocolListParams) {
  const [templates, setTemplates] = useState<Protocol[]>([])
  const [loading, setLoading] = useState(true)
  const [error, setError] = useState<Error | null>(null)
  const [totalCount, setTotalCount] = useState(0)
  const [page, setPage] = useState(initialParams?.page || 1)
  const [perPage, setPerPage] = useState(initialParams?.per_page || 10)
  const [params, setParams] = useState<ProtocolListParams>(initialParams || {})
  
  const protocolService: ProtocolService = getProtocolService()
  
  const fetchTemplates = useCallback(async () => {
    try {
      setLoading(true)
      setError(null)
      
      const queryParams: ProtocolListParams = {
        ...params,
        page,
        per_page: perPage,
        is_template: true
      }
      
      const response = await protocolService.getTemplates(queryParams)
      
      setTemplates(response.data)
      setTotalCount(response.total)
    } catch (err) {
      setError(err instanceof Error ? err : new Error('Failed to fetch protocol templates'))
    } finally {
      setLoading(false)
    }
  }, [protocolService, params, page, perPage])
  
  useEffect(() => {
    fetchTemplates()
  }, [fetchTemplates])
  
  const updateParams = useCallback((newParams: Partial<ProtocolListParams>) => {
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
  
  const createFromTemplate = useCallback(async (
    templateId: UUID,
    customizations?: Partial<Protocol>
  ) => {
    try {
      setLoading(true)
      setError(null)
      
      const newProtocol = await protocolService.createFromTemplate(templateId, customizations)
      return newProtocol
    } catch (err) {
      setError(err instanceof Error ? err : new Error('Failed to create protocol from template'))
      return null
    } finally {
      setLoading(false)
    }
  }, [protocolService])
  
  return {
    templates,
    loading,
    error,
    totalCount,
    page,
    perPage,
    totalPages: Math.ceil(totalCount / perPage),
    updateParams,
    changePage,
    changePerPage,
    createFromTemplate,
    refreshTemplates: fetchTemplates
  }
}

export function useProtocolLibrary(initialParams?: ProtocolListParams) {
  const [libraryProtocols, setLibraryProtocols] = useState<ProtocolLibraryItem[]>([])
  const [loading, setLoading] = useState(true)
  const [error, setError] = useState<Error | null>(null)
  const [totalCount, setTotalCount] = useState(0)
  const [page, setPage] = useState(initialParams?.page || 1)
  const [perPage, setPerPage] = useState(initialParams?.per_page || 10)
  const [params, setParams] = useState<ProtocolListParams>(initialParams || {})
  
  const protocolService: ProtocolService = getProtocolService()
  
  const fetchLibraryProtocols = useCallback(async () => {
    try {
      setLoading(true)
      setError(null)
      
      const queryParams: ProtocolListParams = {
        ...params,
        page,
        per_page: perPage
      }
      
      const response = await protocolService.getLibraryProtocols(queryParams)
      
      setLibraryProtocols(response.data)
      setTotalCount(response.total)
    } catch (err) {
      setError(err instanceof Error ? err : new Error('Failed to fetch library protocols'))
    } finally {
      setLoading(false)
    }
  }, [protocolService, params, page, perPage])
  
  useEffect(() => {
    fetchLibraryProtocols()
  }, [fetchLibraryProtocols])
  
  const updateParams = useCallback((newParams: Partial<ProtocolListParams>) => {
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
  
  const getLibraryProtocol = useCallback(async (id: UUID) => {
    try {
      setLoading(true)
      setError(null)
      
      const protocol = await protocolService.getLibraryProtocol(id)
      return protocol
    } catch (err) {
      setError(err instanceof Error ? err : new Error('Failed to fetch library protocol'))
      return null
    } finally {
      setLoading(false)
    }
  }, [protocolService])
  
  return {
    libraryProtocols,
    loading,
    error,
    totalCount,
    page,
    perPage,
    totalPages: Math.ceil(totalCount / perPage),
    updateParams,
    changePage,
    changePerPage,
    getLibraryProtocol,
    refreshLibraryProtocols: fetchLibraryProtocols
  }
}

export function useProtocolCreation() {
  const [creating, setCreating] = useState(false)
  const [error, setError] = useState<Error | null>(null)
  
  const protocolService: ProtocolService = getProtocolService()
  
  const createProtocol = useCallback(async (protocol: Omit<Protocol, 'id' | 'provenance'>) => {
    try {
      setCreating(true)
      setError(null)
      
      const newProtocol = await protocolService.createProtocol(protocol)
      return newProtocol
    } catch (err) {
      setError(err instanceof Error ? err : new Error('Failed to create protocol'))
      return null
    } finally {
      setCreating(false)
    }
  }, [protocolService])
  
  const importProtocol = useCallback(async (file: File) => {
    try {
      setCreating(true)
      setError(null)
      
      const importedProtocol = await protocolService.importProtocol(file)
      return importedProtocol
    } catch (err) {
      setError(err instanceof Error ? err : new Error('Failed to import protocol'))
      return null
    } finally {
      setCreating(false)
    }
  }, [protocolService])
  
  return {
    creating,
    error,
    createProtocol,
    importProtocol
  }
}

export function useProtocolSearch() {
  const [results, setResults] = useState<Protocol[]>([])
  const [loading, setLoading] = useState(false)
  const [error, setError] = useState<Error | null>(null)
  const [totalCount, setTotalCount] = useState(0)
  
  const protocolService: ProtocolService = getProtocolService()
  
  const searchProtocols = useCallback(async (query: string, params?: ProtocolListParams) => {
    try {
      setLoading(true)
      setError(null)
      
      const response = await protocolService.searchProtocols(query, params)
      
      setResults(response.data)
      setTotalCount(response.total)
      return response
    } catch (err) {
      setError(err instanceof Error ? err : new Error('Failed to search protocols'))
      return { data: [], total: 0 }
    } finally {
      setLoading(false)
    }
  }, [protocolService])
  
  return {
    results,
    loading,
    error,
    totalCount,
    searchProtocols
  }
}

// Additional custom hooks for protocol-specific functionality

export function useProtocolDesign(mixtureId: UUID | null) {
  const [loading, setLoading] = useState(false)
  const [error, setError] = useState<Error | null>(null)
  
  const apiProtocolService = getProtocolService() as ApiProtocolService
  
  const designProtocol = useCallback(async (parameters: {
    target_concentration: number;
    sample_type: string;
    starting_temperature: number;
    target_temperature?: number;
    step_count?: number;
    custom_sensitivity?: Record<string, any>;
  }) => {
    if (!mixtureId) return null
    
    try {
      setLoading(true)
      setError(null)
      
      // Use the API-specific method
      const protocol = await apiProtocolService.designProtocol(
        mixtureId,
        parameters
      )
      return protocol
    } catch (err) {
      setError(err instanceof Error ? err : new Error('Failed to design protocol'))
      return null
    } finally {
      setLoading(false)
    }
  }, [mixtureId, apiProtocolService])
  
  const saveProtocol = useCallback(async (protocol: {
    name: string;
    description?: string;
    target_concentration: number;
    sample_type: string;
    starting_temperature: number;
    target_temperature?: number;
    step_count?: number;
    steps?: any[];
    custom_sensitivity?: Record<string, any>;
  }) => {
    if (!mixtureId) return null
    
    try {
      setLoading(true)
      setError(null)
      
      // Use the API-specific method
      const savedProtocol = await apiProtocolService.saveProtocol(
        mixtureId,
        protocol
      )
      return savedProtocol
    } catch (err) {
      setError(err instanceof Error ? err : new Error('Failed to save protocol'))
      return null
    } finally {
      setLoading(false)
    }
  }, [mixtureId, apiProtocolService])
  
  const getSensitivityProfiles = useCallback(async () => {
    try {
      setLoading(true)
      setError(null)
      
      // Use the API-specific method
      const profiles = await apiProtocolService.getSensitivityProfiles()
      return profiles
    } catch (err) {
      setError(err instanceof Error ? err : new Error('Failed to get sensitivity profiles'))
      return null
    } finally {
      setLoading(false)
    }
  }, [apiProtocolService])
  
  return {
    loading,
    error,
    designProtocol,
    saveProtocol,
    getSensitivityProfiles
  }
}

export function useProtocolsForMixture(mixtureId: UUID | null) {
  const [protocols, setProtocols] = useState<Protocol[]>([])
  const [mixtureName, setMixtureName] = useState<string>('')
  const [loading, setLoading] = useState(false)
  const [error, setError] = useState<Error | null>(null)
  
  const apiProtocolService = getProtocolService() as ApiProtocolService
  
  const fetchProtocols = useCallback(async () => {
    if (!mixtureId) return
    
    try {
      setLoading(true)
      setError(null)
      
      // Use the API-specific method
      const response = await apiProtocolService.getProtocolsForMixture(mixtureId)
      setProtocols(response.protocols)
      setMixtureName(response.mixture_name)
    } catch (err) {
      setError(err instanceof Error ? err : new Error('Failed to fetch protocols for mixture'))
    } finally {
      setLoading(false)
    }
  }, [mixtureId, apiProtocolService])
  
  useEffect(() => {
    fetchProtocols()
  }, [fetchProtocols])
  
  return {
    protocols,
    mixtureName,
    loading,
    error,
    refreshProtocols: fetchProtocols
  }
}

export function useProtocolComparison() {
  const [loading, setLoading] = useState(false)
  const [error, setError] = useState<Error | null>(null)
  
  const apiProtocolService = getProtocolService() as ApiProtocolService
  
  const compareProtocols = useCallback(async (protocolIds: UUID[]) => {
    if (!protocolIds.length || protocolIds.length < 2) {
      throw new Error('At least two protocols must be selected for comparison')
    }
    
    try {
      setLoading(true)
      setError(null)
      
      // Use the API-specific method
      const comparison = await apiProtocolService.compareProtocols(protocolIds)
      return comparison
    } catch (err) {
      setError(err instanceof Error ? err : new Error('Failed to compare protocols'))
      return null
    } finally {
      setLoading(false)
    }
  }, [apiProtocolService])
  
  return {
    loading,
    error,
    compareProtocols
  }
}
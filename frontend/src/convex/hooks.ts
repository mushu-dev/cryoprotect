'use client'

import { useQuery, useMutation } from 'convex/react'
import { useRouter } from 'next/navigation'

/**
 * Hooks for accessing and interacting with Convex data
 * This file provides convenient hooks that wrap the Convex API
 */

// NOTE: The api import will be generated when you run 'npm run convex:codegen'
// For now, we're using a type assertion as a placeholder
const api = {} as any // Replace with actual import when generated

/**
 * Hook to get molecules with pagination
 */
export function useMolecules({ limit = 10, offset = 0 } = {}) {
  const molecules = useQuery(api?.molecules?.query?.list, { 
    limit, 
    skip: offset 
  })
  
  return {
    molecules,
    loading: molecules === undefined,
    error: null // Convex throws errors directly
  }
}

/**
 * Hook to get a single molecule by ID
 */
export function useMolecule(id: string) {
  const molecule = useQuery(api?.molecules?.query?.getById, { id })
  
  return {
    molecule,
    loading: molecule === undefined,
    error: null
  }
}

/**
 * Hook to create a molecule
 */
export function useCreateMolecule() {
  const router = useRouter()
  const createMolecule = useMutation(api?.molecules?.mutations?.create)
  
  const create = async (data: any) => {
    const id = await createMolecule(data)
    return id
  }
  
  const createAndNavigate = async (data: any) => {
    const id = await createMolecule(data)
    router.push(`/molecules/${id}`)
    return id
  }
  
  return { create, createAndNavigate }
}

/**
 * Hook to update a molecule
 */
export function useUpdateMolecule() {
  const updateMolecule = useMutation(api?.molecules?.mutations?.update)
  
  return async (id: string, data: any) => {
    return await updateMolecule({ id, ...data })
  }
}

/**
 * Hook to delete a molecule
 */
export function useDeleteMolecule() {
  const router = useRouter()
  const deleteMolecule = useMutation(api?.molecules?.mutations?.remove)
  
  const remove = async (id: string) => {
    return await deleteMolecule({ id })
  }
  
  const removeAndNavigate = async (id: string, redirectUrl = '/molecules') => {
    await deleteMolecule({ id })
    router.push(redirectUrl)
  }
  
  return { remove, removeAndNavigate }
}

/**
 * Hook to search molecules
 */
export function useSearchMolecules(query: string, options = {}) {
  const results = useQuery(api?.molecules?.query?.search, { 
    query,
    ...options
  })
  
  return {
    results,
    loading: results === undefined,
    error: null
  }
}

/**
 * Hook to get properties for a molecule
 */
export function useMoleculeProperties(moleculeId: string) {
  const properties = useQuery(api?.properties?.query?.getByMoleculeId, { 
    moleculeId 
  })
  
  return {
    properties,
    loading: properties === undefined,
    error: null
  }
}

/**
 * Hook to get the current user's profile
 */
export function useUserProfile() {
  const profile = useQuery(api?.users?.query?.getCurrentUser)
  
  return {
    profile,
    loading: profile === undefined,
    error: null,
    isAuthenticated: !!profile
  }
}

/**
 * Hook to check if the current user has a specific permission
 */
export function useHasPermission(permission: string) {
  const profile = useQuery(api?.users?.query?.getCurrentUser)
  
  // This is a simplified implementation - customize based on your permission model
  const hasPermission = profile?.permissions?.includes(permission) || false
  
  return {
    hasPermission,
    loading: profile === undefined,
    error: null
  }
}

// Add more hooks as needed for your application
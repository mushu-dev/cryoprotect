import { useQuery, useMutation } from 'convex/react';
import { api } from '../../convex/_generated/api';

/**
 * Convex React hooks for data fetching and mutation.
 * 
 * These hooks provide a simple interface to interact with the Convex backend
 * for querying and mutating data. They automatically handle loading states,
 * errors, and optimistic updates.
 */

// ================ MOLECULES HOOKS ================

/**
 * Hook to fetch recent molecules
 * @param {number} limit - Number of molecules to fetch (default: 20)
 * @returns {Array} - List of molecules and loading state
 */
export function useMolecules(limit = 20) {
  return useQuery(api.molecules.query.getRecentMolecules, { limit });
}

/**
 * Hook to fetch a molecule by ID
 * @param {string} id - Molecule ID
 * @returns {Object} - Molecule data and loading state
 */
export function useMoleculeById(id) {
  return useQuery(api.molecules.query.getMolecule, { id });
}

/**
 * Hook to fetch molecules by IDs
 * @param {Array<string>} ids - Array of molecule IDs
 * @returns {Array} - List of molecules and loading state
 */
export function useMoleculesByIds(ids) {
  return useQuery(api.molecules.query.getMoleculesByIds, { ids });
}

/**
 * Hook to search molecules
 * @param {Object} filter - Search criteria
 * @param {Object} options - Search options (limit, offset, etc.)
 * @returns {Array} - List of molecules and loading state
 */
export function useSearchMolecules(filter, options = {}) {
  return useQuery(api.molecules.query.searchMolecules, { filter, options });
}

/**
 * Hook to add a new molecule
 * @returns {Function} - Mutation function to add a molecule
 */
export function useAddMolecule() {
  return useMutation(api.molecules.mutations.addMolecule);
}

/**
 * Hook to update a molecule
 * @returns {Function} - Mutation function to update a molecule
 */
export function useUpdateMolecule() {
  return useMutation(api.molecules.mutations.updateMolecule);
}

/**
 * Hook to delete a molecule
 * @returns {Function} - Mutation function to delete a molecule
 */
export function useDeleteMolecule() {
  return useMutation(api.molecules.mutations.deleteMolecule);
}

// ================ MIXTURES HOOKS ================

/**
 * Hook to fetch recent mixtures
 * @param {number} limit - Number of mixtures to fetch (default: 20)
 * @returns {Array} - List of mixtures and loading state
 */
export function useMixtures(limit = 20) {
  return useQuery(api.mixtures.mixtures.getRecentMixtures, { limit });
}

/**
 * Hook to fetch a mixture by ID
 * @param {string} id - Mixture ID
 * @returns {Object} - Mixture data and loading state
 */
export function useMixtureById(id) {
  return useQuery(api.mixtures.mixtures.getMixture, { id });
}

/**
 * Hook to fetch mixtures by IDs
 * @param {Array<string>} ids - Array of mixture IDs
 * @returns {Array} - List of mixtures and loading state
 */
export function useMixturesByIds(ids) {
  return useQuery(api.mixtures.mixtures.getMixturesByIds, { ids });
}

/**
 * Hook to search mixtures
 * @param {Object} filter - Search criteria
 * @param {Object} options - Search options (limit, offset, etc.)
 * @returns {Array} - List of mixtures and loading state
 */
export function useSearchMixtures(filter, options = {}) {
  return useQuery(api.mixtures.mixtures.searchMixtures, { filter, options });
}

/**
 * Hook to add a new mixture
 * @returns {Function} - Mutation function to add a mixture
 */
export function useAddMixture() {
  return useMutation(api.mixtures.mutations.addMixture);
}

/**
 * Hook to update a mixture
 * @returns {Function} - Mutation function to update a mixture
 */
export function useUpdateMixture() {
  return useMutation(api.mixtures.mutations.updateMixture);
}

/**
 * Hook to delete a mixture
 * @returns {Function} - Mutation function to delete a mixture
 */
export function useDeleteMixture() {
  return useMutation(api.mixtures.mutations.deleteMixture);
}

// ================ EXPERIMENTS HOOKS ================

/**
 * Hook to fetch recent experiments
 * @param {number} limit - Number of experiments to fetch (default: 20)
 * @returns {Array} - List of experiments and loading state
 */
export function useExperiments(limit = 20) {
  return useQuery(api.experiments.experiments.getRecentExperiments, { limit });
}

/**
 * Hook to fetch an experiment by ID
 * @param {string} id - Experiment ID
 * @returns {Object} - Experiment data and loading state
 */
export function useExperimentById(id) {
  return useQuery(api.experiments.experiments.getExperiment, { id });
}

// ================ PROTOCOLS HOOKS ================

/**
 * Hook to fetch recent protocols
 * @param {number} limit - Number of protocols to fetch (default: 20)
 * @returns {Array} - List of protocols and loading state
 */
export function useProtocols(limit = 20) {
  return useQuery(api.experiments.protocols.getRecentProtocols, { limit });
}

/**
 * Hook to fetch a protocol by ID
 * @param {string} id - Protocol ID
 * @returns {Object} - Protocol data and loading state
 */
export function useProtocolById(id) {
  return useQuery(api.experiments.protocols.getProtocol, { id });
}
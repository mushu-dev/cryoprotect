/**
 * Convex Experiment Data Hook
 * 
 * This hook provides real-time access to experiment data from the Convex database.
 * It supports filtering, sorting, and real-time updates for experiments.
 */
import { useState, useEffect, useCallback, useMemo } from 'react';
import { useQuery, useMutation } from 'convex/react';
import { api } from '../../../convex/_generated/api';

/**
 * Custom hook for fetching and managing experiment data from Convex
 * 
 * @param {Object} initialFilters - Initial filters to apply
 * @param {Object} options - Additional options for fetching data
 * @returns {Object} Experiment data and management functions
 */
export default function useConvexExperimentData(initialFilters = {}, options = {}) {
  // State for filters, selections, and UI
  const [filters, setFilters] = useState(initialFilters);
  const [selectedExperiments, setSelectedExperiments] = useState([]);
  const [viewMode, setViewMode] = useState('grid');
  
  // Default query options
  const queryOptions = {
    limit: options.limit || 50,
    includeResults: options.includeResults !== false,
    includeProtocol: options.includeProtocol !== false,
    includeMixture: options.includeMixture || false,
    includeTissueTypes: options.includeTissueTypes || false,
    includeEquipment: options.includeEquipment || false,
    includeTimeSeries: options.includeTimeSeries || false,
    sortBy: options.sortBy || 'updatedAt',
    sortDirection: options.sortDirection || 'desc'
  };

  // Real-time query for experiments with filters
  const experiments = useQuery(
    api.experiments.enhanced_experiments.listEnhancedExperiments,
    { filter: filters, options: queryOptions }
  ) || [];

  // Convex mutations for experiment management
  const createExperiment = useMutation(api.experiments.enhanced_experiments.createEnhancedExperiment);
  const updateExperiment = useMutation(api.experiments.enhanced_experiments.updateEnhancedExperiment);
  const deleteExperiment = useMutation(api.experiments.enhanced_experiments.deleteEnhancedExperiment);
  const updateExperimentStatus = useMutation(api.experiments.enhanced_experiments.updateEnhancedExperimentStatus);
  
  // Function to fetch a single experiment by ID
  const fetchExperiment = useCallback(async (id) => {
    // Use the query directly since this is an async function
    try {
      const experiment = await api.experiments.enhanced_experiments.getEnhancedExperiment({
        experimentId: id,
        options: {
          includeResults: true,
          includeProtocol: true,
          includeMixture: true,
          includeTissueTypes: true,
          includeEquipment: true,
          includeTimeSeries: true,
        }
      });
      
      return experiment;
    } catch (err) {
      console.error('Error fetching experiment:', err);
      throw err;
    }
  }, []);

  // Toggle experiment selection for comparison
  const toggleExperimentSelection = useCallback((experimentId) => {
    setSelectedExperiments(prev => {
      if (prev.includes(experimentId)) {
        return prev.filter(id => id !== experimentId);
      } else {
        // Limit to comparing max 4 experiments
        if (prev.length >= 4) {
          return [...prev.slice(1), experimentId];
        }
        return [...prev, experimentId];
      }
    });
  }, []);

  // Function to create a new experiment
  const handleCreateExperiment = useCallback(async (experimentData) => {
    try {
      // Transform form data to match the API schema
      const apiData = {
        name: experimentData.name || experimentData.title,
        description: experimentData.description,
        experimentTypeId: experimentData.experimentType || 'cryopreservation',
        temperature: experimentData.temperature ? parseFloat(experimentData.temperature) : undefined,
        temperatureUnit: experimentData.temperature ? '°C' : undefined,
        date: experimentData.date ? new Date(experimentData.date).getTime() : Date.now(),
        status: experimentData.status || 'planned',
        public: experimentData.public !== false,
        // Add any other fields from the form that need transformation
        parameters: {
          cellType: experimentData.cellType,
          freezingRate: experimentData.freezingRate,
          storageTemp: experimentData.storageTemperature,
          thawingMethod: experimentData.thawingMethod,
          concentration: experimentData.concentration,
          notes: experimentData.notes,
          cryoprotectants: experimentData.cryoprotectants,
        }
      };
      
      const experimentId = await createExperiment({ experiment: apiData });
      return experimentId;
    } catch (err) {
      console.error('Error creating experiment:', err);
      throw err;
    }
  }, [createExperiment]);

  // Function to update an existing experiment
  const handleUpdateExperiment = useCallback(async (experimentId, experimentData) => {
    try {
      // Transform form data to match the API schema
      const apiData = {
        name: experimentData.name || experimentData.title,
        description: experimentData.description,
        experimentTypeId: experimentData.experimentType,
        temperature: experimentData.temperature ? parseFloat(experimentData.temperature) : undefined,
        temperatureUnit: experimentData.temperature ? '°C' : undefined,
        status: experimentData.status,
        public: experimentData.public,
        // Add any other fields from the form that need transformation
        parameters: {
          cellType: experimentData.cellType,
          freezingRate: experimentData.freezingRate,
          storageTemp: experimentData.storageTemperature,
          thawingMethod: experimentData.thawingMethod,
          concentration: experimentData.concentration,
          notes: experimentData.notes,
          cryoprotectants: experimentData.cryoprotectants,
        }
      };
      
      await updateExperiment({ experimentId, update: apiData });
      return experimentId;
    } catch (err) {
      console.error('Error updating experiment:', err);
      throw err;
    }
  }, [updateExperiment]);

  // Function to delete an experiment
  const handleDeleteExperiment = useCallback(async (experimentId) => {
    try {
      await deleteExperiment({ experimentId });
      return true;
    } catch (err) {
      console.error('Error deleting experiment:', err);
      throw err;
    }
  }, [deleteExperiment]);

  // Function to update experiment status
  const handleUpdateStatus = useCallback(async (experimentId, status) => {
    try {
      await updateExperimentStatus({ experimentId, status });
      return true;
    } catch (err) {
      console.error('Error updating experiment status:', err);
      throw err;
    }
  }, [updateExperimentStatus]);

  // Get selected experiment data
  const selectedExperimentData = useMemo(() => {
    return selectedExperiments.map(id => 
      experiments.find(exp => exp._id === id)
    ).filter(Boolean);
  }, [experiments, selectedExperiments]);

  // Loading and error states
  const isLoading = experiments === undefined;
  const error = experiments === undefined ? "Failed to load experiments" : null;

  return {
    // Data
    experiments,
    filteredExperiments: experiments, // filtered on the server
    selectedExperiments,
    selectedExperimentData,
    
    // UI state
    loading: isLoading,
    error,
    filters,
    viewMode,
    
    // Functions
    setFilters,
    setViewMode,
    toggleExperimentSelection,
    fetchExperiment,
    
    // CRUD operations
    createExperiment: handleCreateExperiment,
    updateExperiment: handleUpdateExperiment,
    deleteExperiment: handleDeleteExperiment,
    updateExperimentStatus: handleUpdateStatus
  };
}
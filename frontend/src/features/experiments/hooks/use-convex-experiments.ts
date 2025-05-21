/**
 * React hooks for working with experiments using direct Convex integration
 * 
 * These hooks provide real-time data subscriptions and mutations for the
 * enhanced experimental data system, with proper error handling and
 * resilience patterns.
 */

import { useState, useCallback } from 'react';
import { useQuery, useMutation, useAction } from "convex/react";
import { api } from "../../../convex/_generated/api";
import { Id } from "../../../convex/_generated/dataModel";
import { useErrorHandler } from '../../errors/use-error-handler';
import { useToast } from '../../ui/use-toast';

// Type for enhanced experiment query options
export type EnhancedExperimentQueryOptions = {
  includeResults?: boolean;
  includeProtocol?: boolean;
  includeMixture?: boolean;
  includeTissueTypes?: boolean;
  includeEquipment?: boolean;
  includeTimeSeries?: boolean;
};

// Type for enhanced experiment filter options
export type EnhancedExperimentFilter = {
  name?: string;
  experimentTypeId?: string;
  protocolId?: Id<"protocols">;
  mixtureId?: Id<"mixtures">;
  conductedBy?: Id<"users">;
  projectId?: Id<"projects">;
  status?: string;
  dateRange?: {
    start?: number;
    end?: number;
  };
  tags?: string[];
  public?: boolean;
  tissueTypeId?: Id<"tissueTypes">;
};

// Type for list options
export type EnhancedExperimentListOptions = {
  limit?: number;
  cursor?: string;
  sortBy?: string;
  sortDirection?: "asc" | "desc";
} & EnhancedExperimentQueryOptions;

/**
 * Hook for listing enhanced experiments with real-time updates
 */
export function useEnhancedExperiments(
  filter?: EnhancedExperimentFilter,
  options?: EnhancedExperimentListOptions
) {
  const { handleError } = useErrorHandler();
  
  // Use Convex query with auto-subscription
  const experiments = useQuery(
    api.experiments.enhanced_index.listEnhancedExperiments,
    filter ? { filter, options } : { options }
  );
  
  return {
    experiments: experiments || [],
    loading: experiments === undefined,
    error: null, // Convex handles errors differently - they throw and are caught by ErrorBoundary
    pagination: {
      // If we had pagination info, we'd include it here
      hasMore: experiments?.length === (options?.limit || 50)
    }
  };
}

/**
 * Hook for searching enhanced experiments
 */
export function useSearchEnhancedExperiments() {
  const [results, setResults] = useState<any[]>([]);
  const [loading, setLoading] = useState(false);
  const { handleError } = useErrorHandler();
  const { toast } = useToast();
  
  // Use Convex query
  const searchExperiments = useAction(
    api.experiments.enhanced_index.searchEnhancedExperiments
  );
  
  const search = useCallback(async (query: string, limit?: number) => {
    setLoading(true);
    
    try {
      const searchResults = await searchExperiments({ query, limit });
      setResults(searchResults || []);
      return searchResults;
    } catch (error) {
      handleError(error, "Failed to search experiments");
      toast({
        title: "Search Failed",
        description: "There was an error searching experiments. Please try again.",
        variant: "destructive"
      });
      return [];
    } finally {
      setLoading(false);
    }
  }, [searchExperiments, handleError, toast]);
  
  return {
    search,
    results,
    loading
  };
}

/**
 * Hook for getting a single enhanced experiment by ID
 */
export function useEnhancedExperiment(
  experimentId?: Id<"enhancedExperiments">,
  options?: EnhancedExperimentQueryOptions
) {
  // Do not execute query if experimentId is undefined
  const args = experimentId ? { experimentId, options } : "skip";
  
  // Use Convex query with auto-subscription
  const experiment = useQuery(
    api.experiments.enhanced_experiments.getEnhancedExperiment,
    args
  );
  
  return {
    experiment: experiment || null,
    loading: experimentId !== undefined && experiment === undefined,
    error: null // Convex handles errors differently
  };
}

/**
 * Hook for creating a new enhanced experiment
 */
export function useCreateEnhancedExperiment() {
  const [loading, setLoading] = useState(false);
  const { handleError } = useErrorHandler();
  const { toast } = useToast();
  
  // Use Convex mutation
  const createExperimentMutation = useMutation(
    api.experiments.enhanced_experiments.createEnhancedExperiment
  );
  
  const createExperiment = useCallback(async (experimentData: any) => {
    setLoading(true);
    
    try {
      const experimentId = await createExperimentMutation({
        experiment: experimentData
      });
      
      toast({
        title: "Experiment Created",
        description: "Your experiment has been created successfully.",
      });
      
      return experimentId;
    } catch (error) {
      handleError(error, "Failed to create experiment");
      toast({
        title: "Creation Failed",
        description: "There was an error creating the experiment. Please try again.",
        variant: "destructive"
      });
      throw error;
    } finally {
      setLoading(false);
    }
  }, [createExperimentMutation, handleError, toast]);
  
  return {
    createExperiment,
    loading
  };
}

/**
 * Hook for updating an enhanced experiment
 */
export function useUpdateEnhancedExperiment() {
  const [loading, setLoading] = useState(false);
  const { handleError } = useErrorHandler();
  const { toast } = useToast();
  
  // Use Convex mutation
  const updateExperimentMutation = useMutation(
    api.experiments.enhanced_experiments.updateEnhancedExperiment
  );
  
  const updateExperiment = useCallback(async (
    experimentId: Id<"enhancedExperiments">,
    update: any
  ) => {
    setLoading(true);
    
    try {
      const result = await updateExperimentMutation({
        experimentId,
        update
      });
      
      toast({
        title: "Experiment Updated",
        description: "Your experiment has been updated successfully.",
      });
      
      return result;
    } catch (error) {
      handleError(error, "Failed to update experiment");
      toast({
        title: "Update Failed",
        description: "There was an error updating the experiment. Please try again.",
        variant: "destructive"
      });
      throw error;
    } finally {
      setLoading(false);
    }
  }, [updateExperimentMutation, handleError, toast]);
  
  return {
    updateExperiment,
    loading
  };
}

/**
 * Hook for deleting an enhanced experiment
 */
export function useDeleteEnhancedExperiment() {
  const [loading, setLoading] = useState(false);
  const { handleError } = useErrorHandler();
  const { toast } = useToast();
  
  // Use Convex mutation
  const deleteExperimentMutation = useMutation(
    api.experiments.enhanced_experiments.deleteEnhancedExperiment
  );
  
  const deleteExperiment = useCallback(async (
    experimentId: Id<"enhancedExperiments">
  ) => {
    setLoading(true);
    
    try {
      await deleteExperimentMutation({
        experimentId
      });
      
      toast({
        title: "Experiment Deleted",
        description: "The experiment has been deleted successfully.",
      });
      
      return true;
    } catch (error) {
      handleError(error, "Failed to delete experiment");
      toast({
        title: "Deletion Failed",
        description: "There was an error deleting the experiment. Please try again.",
        variant: "destructive"
      });
      throw error;
    } finally {
      setLoading(false);
    }
  }, [deleteExperimentMutation, handleError, toast]);
  
  return {
    deleteExperiment,
    loading
  };
}

/**
 * Hook for updating enhanced experiment status
 */
export function useUpdateEnhancedExperimentStatus() {
  const [loading, setLoading] = useState(false);
  const { handleError } = useErrorHandler();
  const { toast } = useToast();
  
  // Use Convex mutation
  const updateStatusMutation = useMutation(
    api.experiments.enhanced_experiments.updateEnhancedExperimentStatus
  );
  
  const updateStatus = useCallback(async (
    experimentId: Id<"enhancedExperiments">,
    status: string
  ) => {
    setLoading(true);
    
    try {
      const result = await updateStatusMutation({
        experimentId,
        status
      });
      
      toast({
        title: "Status Updated",
        description: `Experiment status updated to "${status}".`,
      });
      
      return result;
    } catch (error) {
      handleError(error, "Failed to update experiment status");
      toast({
        title: "Status Update Failed",
        description: "There was an error updating the experiment status. Please try again.",
        variant: "destructive"
      });
      throw error;
    } finally {
      setLoading(false);
    }
  }, [updateStatusMutation, handleError, toast]);
  
  return {
    updateStatus,
    loading
  };
}

/**
 * Hook for working with experiment results
 */
export function useEnhancedExperimentResults(
  experimentId?: Id<"enhancedExperiments">
) {
  // Use Convex query to get all results for an experiment
  const args = experimentId ? { experimentId } : "skip";
  
  const results = useQuery(
    api.experiments.enhanced_experiment_results.getExperimentResults,
    args
  );
  
  // Mutation for adding a new result
  const addResultMutation = useMutation(
    api.experiments.enhanced_experiment_results.addExperimentResult
  );
  
  // Mutation for updating a result
  const updateResultMutation = useMutation(
    api.experiments.enhanced_experiment_results.updateExperimentResult
  );
  
  // Mutation for deleting a result
  const deleteResultMutation = useMutation(
    api.experiments.enhanced_experiment_results.deleteExperimentResult
  );
  
  const [loading, setLoading] = useState(false);
  const { handleError } = useErrorHandler();
  const { toast } = useToast();
  
  // Function to add a new result
  const addResult = useCallback(async (resultData: any) => {
    if (!experimentId) {
      toast({
        title: "Error",
        description: "Experiment ID is required",
        variant: "destructive"
      });
      throw new Error("Experiment ID is required");
    }
    
    setLoading(true);
    
    try {
      const result = await addResultMutation({
        result: {
          ...resultData,
          experimentId
        }
      });
      
      toast({
        title: "Result Added",
        description: "The experiment result has been added successfully.",
      });
      
      return result;
    } catch (error) {
      handleError(error, "Failed to add experiment result");
      toast({
        title: "Error",
        description: "Failed to add experiment result. Please try again.",
        variant: "destructive"
      });
      throw error;
    } finally {
      setLoading(false);
    }
  }, [experimentId, addResultMutation, toast, handleError]);
  
  // Function to update a result
  const updateResult = useCallback(async (
    resultId: Id<"enhancedExperimentResults">,
    updateData: any
  ) => {
    setLoading(true);
    
    try {
      const result = await updateResultMutation({
        resultId,
        update: updateData
      });
      
      toast({
        title: "Result Updated",
        description: "The experiment result has been updated successfully.",
      });
      
      return result;
    } catch (error) {
      handleError(error, "Failed to update experiment result");
      toast({
        title: "Error",
        description: "Failed to update experiment result. Please try again.",
        variant: "destructive"
      });
      throw error;
    } finally {
      setLoading(false);
    }
  }, [updateResultMutation, toast, handleError]);
  
  // Function to delete a result
  const deleteResult = useCallback(async (
    resultId: Id<"enhancedExperimentResults">
  ) => {
    setLoading(true);
    
    try {
      await deleteResultMutation({
        resultId
      });
      
      toast({
        title: "Result Deleted",
        description: "The experiment result has been deleted successfully.",
      });
      
      return true;
    } catch (error) {
      handleError(error, "Failed to delete experiment result");
      toast({
        title: "Error",
        description: "Failed to delete experiment result. Please try again.",
        variant: "destructive"
      });
      throw error;
    } finally {
      setLoading(false);
    }
  }, [deleteResultMutation, toast, handleError]);
  
  return {
    results: results || [],
    loading: experimentId !== undefined && (results === undefined || loading),
    addResult,
    updateResult,
    deleteResult
  };
}

/**
 * Hook for working with time series data
 */
export function useExperimentTimeSeries(
  experimentId?: Id<"enhancedExperiments">
) {
  // Get all time series for an experiment
  const args = experimentId ? { experimentId } : "skip";
  
  const timeSeries = useQuery(
    api.experiments.enhanced_index.getExperimentTimeSeries,
    args
  );
  
  // Mutation for creating a new time series
  const createTimeSeriesMutation = useMutation(
    api.experiments.enhanced_index.createTimeSeries
  );
  
  // Mutation for adding data points
  const addTimeSeriesDataMutation = useMutation(
    api.experiments.enhanced_index.addTimeSeriesData
  );
  
  const [loading, setLoading] = useState(false);
  const { handleError } = useErrorHandler();
  const { toast } = useToast();
  
  // Function to create a new time series
  const createTimeSeries = useCallback(async (timeSeriesData: any) => {
    if (!experimentId) {
      toast({
        title: "Error",
        description: "Experiment ID is required",
        variant: "destructive"
      });
      throw new Error("Experiment ID is required");
    }
    
    setLoading(true);
    
    try {
      const result = await createTimeSeriesMutation({
        timeSeries: {
          ...timeSeriesData,
          experimentId
        }
      });
      
      toast({
        title: "Time Series Created",
        description: "The time series has been created successfully.",
      });
      
      return result;
    } catch (error) {
      handleError(error, "Failed to create time series");
      toast({
        title: "Error",
        description: "Failed to create time series. Please try again.",
        variant: "destructive"
      });
      throw error;
    } finally {
      setLoading(false);
    }
  }, [experimentId, createTimeSeriesMutation, toast, handleError]);
  
  // Function to add data points to a time series
  const addTimeSeriesData = useCallback(async (
    timeSeriesId: Id<"timeSeries">,
    dataPoints: any[]
  ) => {
    setLoading(true);
    
    try {
      const result = await addTimeSeriesDataMutation({
        timeSeriesId,
        dataPoints
      });
      
      toast({
        title: "Data Added",
        description: `Added ${dataPoints.length} data points to time series.`,
      });
      
      return result;
    } catch (error) {
      handleError(error, "Failed to add time series data");
      toast({
        title: "Error",
        description: "Failed to add time series data. Please try again.",
        variant: "destructive"
      });
      throw error;
    } finally {
      setLoading(false);
    }
  }, [addTimeSeriesDataMutation, toast, handleError]);
  
  // Function to get data for a specific time series
  const getTimeSeriesData = useCallback((timeSeriesId: Id<"timeSeries">) => {
    return useQuery(
      api.experiments.enhanced_index.getTimeSeriesData,
      { timeSeriesId }
    );
  }, []);
  
  return {
    timeSeries: timeSeries || [],
    loading: experimentId !== undefined && (timeSeries === undefined || loading),
    createTimeSeries,
    addTimeSeriesData,
    getTimeSeriesData
  };
}

/**
 * Hook for protocol management
 */
export function useProtocols() {
  // Get all protocols
  const protocols = useQuery(api.experiments.protocols.listProtocols);
  
  // Get template protocols
  const templateProtocols = useQuery(
    api.experiments.protocols.listProtocolTemplates
  );
  
  // Mutation for creating a protocol
  const createProtocolMutation = useMutation(
    api.experiments.protocols.createProtocol
  );
  
  // Mutation for updating a protocol
  const updateProtocolMutation = useMutation(
    api.experiments.protocols.updateProtocol
  );
  
  const [loading, setLoading] = useState(false);
  const { handleError } = useErrorHandler();
  const { toast } = useToast();
  
  // Function to get a specific protocol
  const getProtocol = useCallback((protocolId: Id<"protocols">) => {
    return useQuery(
      api.experiments.protocols.getProtocol,
      { protocolId }
    );
  }, []);
  
  // Function to create a protocol
  const createProtocol = useCallback(async (protocolData: any) => {
    setLoading(true);
    
    try {
      const result = await createProtocolMutation({
        protocol: protocolData
      });
      
      toast({
        title: "Protocol Created",
        description: "The protocol has been created successfully.",
      });
      
      return result;
    } catch (error) {
      handleError(error, "Failed to create protocol");
      toast({
        title: "Error",
        description: "Failed to create protocol. Please try again.",
        variant: "destructive"
      });
      throw error;
    } finally {
      setLoading(false);
    }
  }, [createProtocolMutation, toast, handleError]);
  
  // Function to update a protocol
  const updateProtocol = useCallback(async (
    protocolId: Id<"protocols">,
    updateData: any
  ) => {
    setLoading(true);
    
    try {
      const result = await updateProtocolMutation({
        protocolId,
        update: updateData
      });
      
      toast({
        title: "Protocol Updated",
        description: "The protocol has been updated successfully.",
      });
      
      return result;
    } catch (error) {
      handleError(error, "Failed to update protocol");
      toast({
        title: "Error",
        description: "Failed to update protocol. Please try again.",
        variant: "destructive"
      });
      throw error;
    } finally {
      setLoading(false);
    }
  }, [updateProtocolMutation, toast, handleError]);
  
  return {
    protocols: protocols || [],
    templateProtocols: templateProtocols || [],
    loading: (protocols === undefined || templateProtocols === undefined || loading),
    getProtocol,
    createProtocol,
    updateProtocol
  };
}

/**
 * Hook for lab verification functionality
 */
export function useLabVerification(
  experimentId?: Id<"enhancedExperiments">
) {
  // Get verification for a specific experiment
  const args = experimentId ? { experimentId } : "skip";
  
  const verification = useQuery(
    api.experiments.enhanced_validation.getLabVerificationByExperiment,
    args
  );
  
  // Mutation for requesting verification
  const requestVerificationMutation = useMutation(
    api.experiments.enhanced_validation.requestLabVerification
  );
  
  // Mutation for updating verification status
  const updateVerificationMutation = useMutation(
    api.experiments.enhanced_validation.updateLabVerification
  );
  
  const [loading, setLoading] = useState(false);
  const { handleError } = useErrorHandler();
  const { toast } = useToast();
  
  // Function to request verification
  const requestVerification = useCallback(async (verificationData: any) => {
    if (!experimentId) {
      toast({
        title: "Error",
        description: "Experiment ID is required",
        variant: "destructive"
      });
      throw new Error("Experiment ID is required");
    }
    
    setLoading(true);
    
    try {
      const result = await requestVerificationMutation({
        verification: {
          ...verificationData,
          experimentId,
          verificationStatus: "pending"
        }
      });
      
      toast({
        title: "Verification Requested",
        description: "The lab verification has been requested successfully.",
      });
      
      return result;
    } catch (error) {
      handleError(error, "Failed to request verification");
      toast({
        title: "Error",
        description: "Failed to request verification. Please try again.",
        variant: "destructive"
      });
      throw error;
    } finally {
      setLoading(false);
    }
  }, [experimentId, requestVerificationMutation, toast, handleError]);
  
  // Function to update verification status
  const updateVerification = useCallback(async (
    verificationId: Id<"labVerifications">,
    status: string,
    notes?: string
  ) => {
    setLoading(true);
    
    try {
      const result = await updateVerificationMutation({
        verificationId,
        status,
        notes
      });
      
      toast({
        title: "Verification Updated",
        description: `The verification status has been updated to "${status}".`,
      });
      
      return result;
    } catch (error) {
      handleError(error, "Failed to update verification");
      toast({
        title: "Error",
        description: "Failed to update verification. Please try again.",
        variant: "destructive"
      });
      throw error;
    } finally {
      setLoading(false);
    }
  }, [updateVerificationMutation, toast, handleError]);
  
  return {
    verification: verification || null,
    loading: experimentId !== undefined && (verification === undefined || loading),
    requestVerification,
    updateVerification
  };
}

/**
 * Hook for equipment tracking
 */
export function useEquipment() {
  // Get all equipment
  const equipment = useQuery(api.experiments.enhanced_index.listEquipment);
  
  // Get equipment types
  const equipmentTypes = useQuery(
    api.experiments.enhanced_index.listEquipmentTypes
  );
  
  // Mutation for adding equipment to an experiment
  const addEquipmentToExperimentMutation = useMutation(
    api.experiments.enhanced_index.addEquipmentToExperiment
  );
  
  const [loading, setLoading] = useState(false);
  const { handleError } = useErrorHandler();
  const { toast } = useToast();
  
  // Function to add equipment to an experiment
  const addEquipmentToExperiment = useCallback(async (
    experimentId: Id<"enhancedExperiments">,
    equipmentId: Id<"equipment">,
    role?: string,
    parameters?: Record<string, any>
  ) => {
    setLoading(true);
    
    try {
      const result = await addEquipmentToExperimentMutation({
        experimentId,
        equipmentId,
        role,
        parameters
      });
      
      toast({
        title: "Equipment Added",
        description: "The equipment has been added to the experiment.",
      });
      
      return result;
    } catch (error) {
      handleError(error, "Failed to add equipment to experiment");
      toast({
        title: "Error",
        description: "Failed to add equipment to experiment. Please try again.",
        variant: "destructive"
      });
      throw error;
    } finally {
      setLoading(false);
    }
  }, [addEquipmentToExperimentMutation, toast, handleError]);
  
  // Function to get equipment for a specific experiment
  const getExperimentEquipment = useCallback((experimentId: Id<"enhancedExperiments">) => {
    return useQuery(
      api.experiments.enhanced_index.getExperimentEquipment,
      { experimentId }
    );
  }, []);
  
  return {
    equipment: equipment || [],
    equipmentTypes: equipmentTypes || [],
    loading: (equipment === undefined || equipmentTypes === undefined || loading),
    addEquipmentToExperiment,
    getExperimentEquipment
  };
}
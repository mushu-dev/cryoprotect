/**
 * Enhanced Experimental Data System
 * 
 * This module exports all functionalities for the enhanced experiment system including:
 * - Protocols and protocol templates
 * - Enhanced experiments with comprehensive metadata
 * - Experimental results with uncertainty quantification
 * - Time series data
 * - Equipment tracking
 * - Validation rules
 */

// Re-export all enhanced experiment types
export * from "./enhanced_types";

// Re-export core experiment operations
export {
  createEnhancedExperiment,
  getEnhancedExperiment,
  updateEnhancedExperiment,
  deleteEnhancedExperiment,
  listEnhancedExperiments,
  searchEnhancedExperiments,
  updateEnhancedExperimentStatus
} from "./enhanced_experiments";

// Re-export experiment result operations
export {
  createEnhancedExperimentResult,
  getEnhancedExperimentResult,
  updateEnhancedExperimentResult,
  deleteEnhancedExperimentResult,
  listEnhancedExperimentResults
} from "./enhanced_experiment_results";

// Re-export protocol operations
export {
  createProtocol,
  getProtocol,
  updateProtocol,
  createProtocolVersion,
  deleteProtocol,
  listProtocols,
  getProtocolVersionHistory
} from "./protocols";

// Export validation utilities
export {
  validateCreateProtocolInput,
  validateUpdateProtocolInput,
  validateCreateTissueTypeInput,
  validateUpdateTissueTypeInput,
  validateCreateEnhancedExperimentInput,
  validateUpdateEnhancedExperimentInput,
  validateCreateEnhancedExperimentResultInput,
  validateUpdateEnhancedExperimentResultInput,
  validateUncertainty,
  validateCreateTimeSeriesInput,
  validateCreateTimeSeriesDataPointInput,
  validateCreateEquipmentInput,
  validateCreateExperimentEquipmentInput,
  validateCreateValidationRuleInput
} from "./enhanced_validation";

// Export helper functions
export {
  expandEnhancedExperimentWithResults,
  expandEnhancedExperimentWithProtocol,
  expandEnhancedExperimentWithMixture,
  expandEnhancedExperimentWithTissueTypes,
  expandEnhancedExperimentWithEquipment,
  expandEnhancedExperimentWithTimeSeries,
  expandEnhancedExperimentWithDetails,
  expandTimeSeriesWithData,
  enhancedExperimentHasResults,
  enhancedExperimentHasTimeSeries,
  enhancedExperimentHasEquipment
} from "./enhanced_helpers";
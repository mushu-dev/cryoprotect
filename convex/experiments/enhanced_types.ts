/**
 * Enhanced type definitions for experiments and related data
 */

import { Id } from "../_generated/dataModel";

/**
 * Protocol - template or instance for experiment procedures
 */
export interface Protocol {
  _id: Id<"protocols">;
  _creationTime: number;
  
  // Core fields
  name: string;
  description?: string;
  steps: ProtocolStep[];
  parameters?: Record<string, any>;
  version: number;
  parentId?: Id<"protocols">;
  isTemplate: boolean;
  category?: string;
  
  // Metadata
  createdBy?: Id<"users">;
  createdAt: number;
  updatedAt: number;
  public: boolean;
}

/**
 * Protocol step - individual procedure in a protocol
 */
export interface ProtocolStep {
  name: string;
  description?: string;
  parameters?: Record<string, any>;
  duration?: number;
  durationUnit?: string;
  temperature?: number;
  temperatureUnit?: string;
}

/**
 * Tissue type - biological sample used in experiments
 */
export interface TissueType {
  _id: Id<"tissueTypes">;
  _creationTime: number;
  
  // Core fields
  name: string;
  description?: string;
  species?: string;
  taxonomyId?: number;
  properties?: Record<string, any>;
  category?: string;
  source?: string;
  
  // Metadata
  createdBy?: Id<"users">;
  createdAt: number;
  updatedAt: number;
  public: boolean;
}

/**
 * Enhanced experiment - extended experiment record
 */
export interface EnhancedExperiment {
  _id: Id<"enhancedExperiments">;
  _creationTime: number;
  
  // Core fields
  name: string;
  description?: string;
  experimentTypeId?: string;
  protocolId?: Id<"protocols">;
  mixtureId?: Id<"mixtures">;
  temperature?: number;
  temperatureUnit?: string;
  coolingRate?: number;
  coolingRateUnit?: string;
  thawingRate?: number;
  thawingRateUnit?: string;
  pressure?: number;
  pressureUnit?: string;
  parameters?: Record<string, any>;
  version: number;
  provenance?: Record<string, any>;
  conductedBy?: Id<"users">;
  projectId?: Id<"projects">;
  date?: number;
  status: "planned" | "in-progress" | "completed" | "failed";
  tags?: string[];
  
  // Metadata
  createdBy?: Id<"users">;
  createdAt: number;
  updatedAt: number;
  public: boolean;
}

/**
 * Enhanced experiment result - extended result with uncertainty
 */
export interface EnhancedExperimentResult {
  _id: Id<"enhancedExperimentResults">;
  _creationTime: number;
  
  // Core fields
  experimentId: Id<"enhancedExperiments">;
  moleculeId?: Id<"molecules">;
  mixtureId?: Id<"mixtures">;
  tissueTypeId: Id<"tissueTypes">;
  parameterName: string;
  value: string | number | boolean | null;
  numericValue?: number;
  units?: string;
  uncertainty?: Uncertainty;
  provenance?: Provenance;
  notes?: string;
  
  // Metadata
  createdBy?: Id<"users">;
  createdAt: number;
  updatedAt: number;
}

/**
 * Uncertainty - captures measurement uncertainty
 */
export interface Uncertainty {
  type: "standard_deviation" | "range" | "confidence_interval";
  value: number | number[]; // Single value for std_dev, array [min, max] for range/CI
  confidence?: number; // For confidence intervals (e.g., 0.95 for 95% CI)
}

/**
 * Provenance - data source information
 */
export interface Provenance {
  method: string;
  reference?: string;
  timestamp: number;
  operator?: string;
}

/**
 * Equipment type - category of laboratory equipment
 */
export interface EquipmentType {
  _id: Id<"equipmentTypes">;
  _creationTime: number;
  
  // Core fields
  name: string;
  description?: string;
  manufacturer?: string;
  category?: string;
  
  // Metadata
  createdBy?: Id<"users">;
  createdAt: number;
  updatedAt: number;
}

/**
 * Equipment - specific piece of laboratory equipment
 */
export interface Equipment {
  _id: Id<"equipment">;
  _creationTime: number;
  
  // Core fields
  equipmentTypeId: Id<"equipmentTypes">;
  name: string;
  model?: string;
  serialNumber?: string;
  description?: string;
  calibrationDate?: number;
  nextCalibrationDate?: number;
  location?: string;
  parameters?: Record<string, any>;
  
  // Metadata
  createdBy?: Id<"users">;
  createdAt: number;
  updatedAt: number;
}

/**
 * Experiment equipment - links experiments to equipment used
 */
export interface ExperimentEquipment {
  _id: Id<"experimentEquipment">;
  _creationTime: number;
  
  // Core fields
  experimentId: Id<"enhancedExperiments">;
  equipmentId: Id<"equipment">;
  role?: string;
  parameters?: Record<string, any>;
  notes?: string;
  
  // Metadata
  createdBy?: Id<"users">;
  createdAt: number;
  updatedAt: number;
}

/**
 * Time series - collection of time-based measurements
 */
export interface TimeSeries {
  _id: Id<"timeSeries">;
  _creationTime: number;
  
  // Core fields
  experimentId: Id<"enhancedExperiments">;
  name: string;
  description?: string;
  parameterName: string;
  units?: string;
  metadata?: Record<string, any>;
  
  // Metadata
  createdBy?: Id<"users">;
  createdAt: number;
  updatedAt: number;
}

/**
 * Time series data point
 */
export interface TimeSeriesDataPoint {
  _id: Id<"timeSeriesData">;
  _creationTime: number;
  
  // Core fields
  timeSeriesId: Id<"timeSeries">;
  timestamp: number;
  value: number;
  uncertainty?: number;
  metadata?: Record<string, any>;
  
  // Metadata
  createdBy?: Id<"users">;
  createdAt: number;
  updatedAt: number;
}

/**
 * Validation rule for experiment data
 */
export interface ValidationRule {
  _id: Id<"validationRules">;
  _creationTime: number;
  
  // Core fields
  parameterName: string;
  ruleType: "range" | "pattern" | "comparison" | "custom";
  parameters: Record<string, any>;
  description?: string;
  severity: "error" | "warning" | "info";
  
  // Metadata
  createdBy?: Id<"users">;
  createdAt: number;
  updatedAt: number;
}

// Input types for creating/updating records

/**
 * Input for creating a protocol
 */
export interface CreateProtocolInput {
  name: string;
  description?: string;
  steps: ProtocolStep[];
  parameters?: Record<string, any>;
  version?: number;
  parentId?: Id<"protocols">;
  isTemplate?: boolean;
  category?: string;
  public?: boolean;
}

/**
 * Input for updating a protocol
 */
export interface UpdateProtocolInput {
  name?: string;
  description?: string;
  steps?: ProtocolStep[];
  parameters?: Record<string, any>;
  version?: number;
  category?: string;
  public?: boolean;
}

/**
 * Input for creating a tissue type
 */
export interface CreateTissueTypeInput {
  name: string;
  description?: string;
  species?: string;
  taxonomyId?: number;
  properties?: Record<string, any>;
  category?: string;
  source?: string;
  public?: boolean;
}

/**
 * Input for updating a tissue type
 */
export interface UpdateTissueTypeInput {
  name?: string;
  description?: string;
  species?: string;
  taxonomyId?: number;
  properties?: Record<string, any>;
  category?: string;
  source?: string;
  public?: boolean;
}

/**
 * Input for creating an enhanced experiment
 */
export interface CreateEnhancedExperimentInput {
  name: string;
  description?: string;
  experimentTypeId?: string;
  protocolId?: Id<"protocols">;
  mixtureId?: Id<"mixtures">;
  temperature?: number;
  temperatureUnit?: string;
  coolingRate?: number;
  coolingRateUnit?: string;
  thawingRate?: number;
  thawingRateUnit?: string;
  pressure?: number;
  pressureUnit?: string;
  parameters?: Record<string, any>;
  version?: number;
  provenance?: Record<string, any>;
  projectId?: Id<"projects">;
  date?: number;
  status?: "planned" | "in-progress" | "completed" | "failed";
  tags?: string[];
  public?: boolean;
}

/**
 * Input for updating an enhanced experiment
 */
export interface UpdateEnhancedExperimentInput {
  name?: string;
  description?: string;
  experimentTypeId?: string;
  protocolId?: Id<"protocols">;
  mixtureId?: Id<"mixtures">;
  temperature?: number;
  temperatureUnit?: string;
  coolingRate?: number;
  coolingRateUnit?: string;
  thawingRate?: number;
  thawingRateUnit?: string;
  pressure?: number;
  pressureUnit?: string;
  parameters?: Record<string, any>;
  version?: number;
  provenance?: Record<string, any>;
  projectId?: Id<"projects">;
  date?: number;
  status?: "planned" | "in-progress" | "completed" | "failed";
  tags?: string[];
  public?: boolean;
}

/**
 * Input for creating an enhanced experiment result
 */
export interface CreateEnhancedExperimentResultInput {
  experimentId: Id<"enhancedExperiments">;
  moleculeId?: Id<"molecules">;
  mixtureId?: Id<"mixtures">;
  tissueTypeId: Id<"tissueTypes">;
  parameterName: string;
  value: string | number | boolean | null;
  numericValue?: number;
  units?: string;
  uncertainty?: Uncertainty;
  provenance?: Provenance;
  notes?: string;
}

/**
 * Input for updating an enhanced experiment result
 */
export interface UpdateEnhancedExperimentResultInput {
  parameterName?: string;
  value?: string | number | boolean | null;
  numericValue?: number;
  units?: string;
  uncertainty?: Uncertainty;
  provenance?: Provenance;
  notes?: string;
}

/**
 * Enhanced experiment with expanded data
 */
export interface EnhancedExperimentWithDetails extends EnhancedExperiment {
  results?: EnhancedExperimentResult[];
  protocol?: Protocol;
  mixture?: {
    _id: Id<"mixtures">;
    name: string;
    components?: Array<{
      moleculeId: Id<"molecules">;
      moleculeName: string;
      concentration: number;
      units: string;
    }>;
  };
  tissueTypes?: TissueType[];
  equipment?: Array<{
    _id: Id<"equipment">;
    name: string;
    type: string;
    role?: string;
  }>;
  timeSeries?: TimeSeries[];
}

/**
 * Enhanced filter options for experiment queries
 */
export interface EnhancedExperimentFilter {
  name?: string;
  experimentTypeId?: string;
  protocolId?: Id<"protocols">;
  mixtureId?: Id<"mixtures">;
  conductedBy?: Id<"users">;
  projectId?: Id<"projects">;
  status?: "planned" | "in-progress" | "completed" | "failed";
  dateRange?: {
    start?: number;
    end?: number;
  };
  tags?: string[];
  public?: boolean;
  parameterName?: string;
  parameterValue?: string | number | boolean;
  tissueTypeId?: Id<"tissueTypes">;
}

/**
 * Query options for enhanced experiments
 */
export interface EnhancedExperimentQueryOptions {
  limit?: number;
  cursor?: string;
  includeResults?: boolean;
  includeProtocol?: boolean;
  includeMixture?: boolean;
  includeTissueTypes?: boolean;
  includeEquipment?: boolean;
  includeTimeSeries?: boolean;
  sortBy?: "name" | "date" | "status" | "updatedAt";
  sortDirection?: "asc" | "desc";
}

/**
 * Time series with data points
 */
export interface TimeSeriesWithData extends TimeSeries {
  dataPoints: TimeSeriesDataPoint[];
}

/**
 * Input for creating a time series
 */
export interface CreateTimeSeriesInput {
  experimentId: Id<"enhancedExperiments">;
  name: string;
  description?: string;
  parameterName: string;
  units?: string;
  metadata?: Record<string, any>;
}

/**
 * Input for creating a time series data point
 */
export interface CreateTimeSeriesDataPointInput {
  timeSeriesId: Id<"timeSeries">;
  timestamp: number;
  value: number;
  uncertainty?: number;
  metadata?: Record<string, any>;
}

/**
 * Input for creating equipment
 */
export interface CreateEquipmentInput {
  equipmentTypeId: Id<"equipmentTypes">;
  name: string;
  model?: string;
  serialNumber?: string;
  description?: string;
  calibrationDate?: number;
  nextCalibrationDate?: number;
  location?: string;
  parameters?: Record<string, any>;
}

/**
 * Input for linking equipment to an experiment
 */
export interface CreateExperimentEquipmentInput {
  experimentId: Id<"enhancedExperiments">;
  equipmentId: Id<"equipment">;
  role?: string;
  parameters?: Record<string, any>;
  notes?: string;
}

/**
 * Input for creating a validation rule
 */
export interface CreateValidationRuleInput {
  parameterName: string;
  ruleType: "range" | "pattern" | "comparison" | "custom";
  parameters: Record<string, any>;
  description?: string;
  severity: "error" | "warning" | "info";
}
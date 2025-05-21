/**
 * Type definitions for predictions and scientific models
 */

import { Id } from "../_generated/dataModel";

/**
 * Scientific model document in the database
 */
export interface ScientificModel {
  _id: Id<"scientificModels">;
  _creationTime: number;
  
  // Core fields
  name: string;
  description?: string;
  version: string;
  type: string;
  parameters?: string;
  createdBy?: Id<"users">;
  
  // Metadata
  createdAt: number;
  updatedAt: number;
  public: boolean;
}

/**
 * Input for creating a new scientific model
 */
export interface CreateScientificModelInput {
  name: string;
  description?: string;
  version: string;
  type: string;
  parameters?: string;
  public?: boolean;
}

/**
 * Input for updating a scientific model
 */
export interface UpdateScientificModelInput {
  name?: string;
  description?: string;
  version?: string;
  type?: string;
  parameters?: string;
  public?: boolean;
}

/**
 * Prediction document in the database
 */
export interface Prediction {
  _id: Id<"predictions">;
  _creationTime: number;
  
  // Core fields
  modelId: Id<"scientificModels">;
  moleculeId?: Id<"molecules">;
  mixtureId?: Id<"mixtures">;
  parameterName: string;
  value: string | number | boolean | null;
  numericValue?: number;
  units?: string;
  confidence?: number;
  
  // Metadata
  createdAt: number;
  updatedAt: number;
}

/**
 * Input for creating a new prediction
 */
export interface CreatePredictionInput {
  modelId: Id<"scientificModels">;
  moleculeId?: Id<"molecules">;
  mixtureId?: Id<"mixtures">;
  parameterName: string;
  value: string | number | boolean | null;
  numericValue?: number;
  units?: string;
  confidence?: number;
}

/**
 * Input for updating a prediction
 */
export interface UpdatePredictionInput {
  parameterName?: string;
  value?: string | number | boolean | null;
  numericValue?: number;
  units?: string;
  confidence?: number;
}

/**
 * Prediction with expanded model data
 */
export interface PredictionWithModel extends Prediction {
  model: {
    _id: Id<"scientificModels">;
    name: string;
    version: string;
    type: string;
  };
}

/**
 * Prediction with expanded target (molecule or mixture) data
 */
export interface PredictionWithTarget extends Prediction {
  model: {
    _id: Id<"scientificModels">;
    name: string;
    version: string;
    type: string;
  };
  molecule?: {
    _id: Id<"molecules">;
    name: string;
    formula?: string;
  };
  mixture?: {
    _id: Id<"mixtures">;
    name: string;
  };
}

/**
 * Filter conditions for querying scientific models
 */
export interface ScientificModelFilter {
  name?: string;
  type?: string;
  version?: string;
  createdBy?: Id<"users">;
  public?: boolean;
}

/**
 * Filter conditions for querying predictions
 */
export interface PredictionFilter {
  modelId?: Id<"scientificModels">;
  moleculeId?: Id<"molecules">;
  mixtureId?: Id<"mixtures">;
  parameterName?: string;
  confidenceThreshold?: number;
  valueEquals?: string | number | boolean;
  valueRange?: {
    min?: number;
    max?: number;
  };
}

/**
 * Options for scientific model queries
 */
export interface ScientificModelQueryOptions {
  limit?: number;
  cursor?: string;
  sortBy?: "name" | "type" | "updatedAt";
  sortDirection?: "asc" | "desc";
}

/**
 * Options for prediction queries
 */
export interface PredictionQueryOptions {
  limit?: number;
  cursor?: string;
  includeModel?: boolean;
  includeTarget?: boolean;
}
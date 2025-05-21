/**
 * Type definitions for experiments and experiment results
 */

import { Id } from "../_generated/dataModel";

/**
 * Experiment document in the database
 */
export interface Experiment {
  _id: Id<"experiments">;
  _creationTime: number;
  
  // Core fields
  name: string;
  description?: string;
  mixtureId?: Id<"mixtures">;
  protocol?: string;
  conductedBy?: Id<"users">;
  projectId?: Id<"projects">;
  date?: number;
  status: "planned" | "in-progress" | "completed" | "failed";
  
  // Metadata
  createdAt: number;
  updatedAt: number;
  public: boolean;
}

/**
 * Input for creating a new experiment
 */
export interface CreateExperimentInput {
  name: string;
  description?: string;
  mixtureId?: Id<"mixtures">;
  protocol?: string;
  projectId?: Id<"projects">;
  date?: number;
  status?: "planned" | "in-progress" | "completed" | "failed";
  public?: boolean;
}

/**
 * Input for updating an experiment
 */
export interface UpdateExperimentInput {
  name?: string;
  description?: string;
  mixtureId?: Id<"mixtures">;
  protocol?: string;
  projectId?: Id<"projects">;
  date?: number;
  status?: "planned" | "in-progress" | "completed" | "failed";
  public?: boolean;
}

/**
 * Experiment result document in the database
 */
export interface ExperimentResult {
  _id: Id<"experimentResults">;
  _creationTime: number;
  
  // Core fields
  experimentId: Id<"experiments">;
  parameterName: string;
  value: string | number | boolean | null;
  numericValue?: number;
  units?: string;
  notes?: string;
  
  // Metadata
  createdAt: number;
  updatedAt: number;
}

/**
 * Input for creating a new experiment result
 */
export interface CreateExperimentResultInput {
  experimentId: Id<"experiments">;
  parameterName: string;
  value: string | number | boolean | null;
  numericValue?: number;
  units?: string;
  notes?: string;
}

/**
 * Input for updating an experiment result
 */
export interface UpdateExperimentResultInput {
  parameterName?: string;
  value?: string | number | boolean | null;
  numericValue?: number;
  units?: string;
  notes?: string;
}

/**
 * Experiment with expanded results data
 */
export interface ExperimentWithResults extends Experiment {
  results: ExperimentResult[];
}

/**
 * Experiment with mixture details
 */
export interface ExperimentWithMixture extends Experiment {
  mixture?: {
    _id: Id<"mixtures">;
    name: string;
    componentCount?: number;
  };
}

/**
 * Experiment with complete details (results and mixture)
 */
export interface ExperimentComplete extends Experiment {
  results: ExperimentResult[];
  mixture?: {
    _id: Id<"mixtures">;
    name: string;
    components: Array<{
      moleculeId: Id<"molecules">;
      moleculeName: string;
      concentration: number;
      units: string;
    }>;
  };
}

/**
 * Filter conditions for querying experiments
 */
export interface ExperimentFilter {
  name?: string;
  mixtureId?: Id<"mixtures">;
  conductedBy?: Id<"users">;
  projectId?: Id<"projects">;
  status?: "planned" | "in-progress" | "completed" | "failed";
  dateRange?: {
    start?: number;
    end?: number;
  };
  public?: boolean;
}

/**
 * Filter conditions for querying experiment results
 */
export interface ExperimentResultFilter {
  experimentId?: Id<"experiments">;
  parameterName?: string;
  valueEquals?: string | number | boolean;
  valueRange?: {
    min?: number;
    max?: number;
  };
}

/**
 * Options for experiment queries
 */
export interface ExperimentQueryOptions {
  limit?: number;
  cursor?: string;
  includeResults?: boolean;
  includeMixture?: boolean;
  sortBy?: "name" | "date" | "status" | "updatedAt";
  sortDirection?: "asc" | "desc";
}

/**
 * Options for experiment result queries
 */
export interface ExperimentResultQueryOptions {
  limit?: number;
  cursor?: string;
}
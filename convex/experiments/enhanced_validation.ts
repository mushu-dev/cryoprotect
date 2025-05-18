/**
 * Enhanced validation functions for experimental data
 */

import { ConvexError } from "convex/values";
import { Id } from "../_generated/dataModel";
import { DatabaseReader } from "../_generated/server";

import {
  CreateProtocolInput,
  UpdateProtocolInput,
  CreateTissueTypeInput,
  UpdateTissueTypeInput,
  CreateEnhancedExperimentInput,
  UpdateEnhancedExperimentInput,
  CreateEnhancedExperimentResultInput,
  UpdateEnhancedExperimentResultInput,
  Uncertainty,
  ProtocolStep,
  CreateTimeSeriesInput,
  CreateTimeSeriesDataPointInput,
  CreateEquipmentInput,
  CreateExperimentEquipmentInput,
  CreateValidationRuleInput
} from "./enhanced_types";

/**
 * Validates protocol steps for required fields and valid values
 */
export function validateProtocolSteps(steps: ProtocolStep[]): void {
  if (!steps || !Array.isArray(steps) || steps.length === 0) {
    throw new ConvexError("Protocol must have at least one step");
  }
  
  for (let i = 0; i < steps.length; i++) {
    const step = steps[i];
    if (!step.name || step.name.trim() === "") {
      throw new ConvexError(`Step ${i + 1} must have a name`);
    }
    
    if (step.duration !== undefined && step.duration < 0) {
      throw new ConvexError(`Step ${i + 1} cannot have negative duration`);
    }
    
    if (step.duration !== undefined && !step.durationUnit) {
      throw new ConvexError(`Step ${i + 1} with duration must specify durationUnit`);
    }
    
    if (step.temperature !== undefined && !step.temperatureUnit) {
      throw new ConvexError(`Step ${i + 1} with temperature must specify temperatureUnit`);
    }
  }
}

/**
 * Validates input for creating a new protocol
 */
export function validateCreateProtocolInput(input: CreateProtocolInput): void {
  if (!input.name || input.name.trim() === "") {
    throw new ConvexError("Protocol name is required");
  }
  
  validateProtocolSteps(input.steps);
  
  if (input.version !== undefined && input.version < 1) {
    throw new ConvexError("Protocol version must be a positive number");
  }
}

/**
 * Validates input for updating a protocol
 */
export function validateUpdateProtocolInput(input: UpdateProtocolInput): void {
  if (input.name !== undefined && input.name.trim() === "") {
    throw new ConvexError("Protocol name cannot be empty");
  }
  
  if (input.steps !== undefined) {
    validateProtocolSteps(input.steps);
  }
  
  if (input.version !== undefined && input.version < 1) {
    throw new ConvexError("Protocol version must be a positive number");
  }
}

/**
 * Validates input for creating a new tissue type
 */
export function validateCreateTissueTypeInput(input: CreateTissueTypeInput): void {
  if (!input.name || input.name.trim() === "") {
    throw new ConvexError("Tissue type name is required");
  }
  
  if (input.taxonomyId !== undefined && input.taxonomyId < 1) {
    throw new ConvexError("Taxonomy ID must be a positive integer");
  }
}

/**
 * Validates input for updating a tissue type
 */
export function validateUpdateTissueTypeInput(input: UpdateTissueTypeInput): void {
  if (input.name !== undefined && input.name.trim() === "") {
    throw new ConvexError("Tissue type name cannot be empty");
  }
  
  if (input.taxonomyId !== undefined && input.taxonomyId < 1) {
    throw new ConvexError("Taxonomy ID must be a positive integer");
  }
}

/**
 * Validates uncertainty object for proper structure and values
 */
export function validateUncertainty(uncertainty: Uncertainty): void {
  const validTypes = ["standard_deviation", "range", "confidence_interval"];
  if (!validTypes.includes(uncertainty.type)) {
    throw new ConvexError(`Uncertainty type must be one of: ${validTypes.join(", ")}`);
  }
  
  if (uncertainty.type === "standard_deviation") {
    if (typeof uncertainty.value !== "number" || uncertainty.value < 0) {
      throw new ConvexError("Standard deviation must be a non-negative number");
    }
  } else if (uncertainty.type === "range" || uncertainty.type === "confidence_interval") {
    if (!Array.isArray(uncertainty.value) || uncertainty.value.length !== 2) {
      throw new ConvexError(`${uncertainty.type} must have array value with [min, max]`);
    }
    
    if (uncertainty.value[0] >= uncertainty.value[1]) {
      throw new ConvexError(`${uncertainty.type} minimum value must be less than maximum`);
    }
  }
  
  if (uncertainty.type === "confidence_interval") {
    if (!uncertainty.confidence) {
      throw new ConvexError("Confidence interval must specify confidence level");
    }
    
    if (uncertainty.confidence <= 0 || uncertainty.confidence >= 1) {
      throw new ConvexError("Confidence level must be between 0 and 1");
    }
  }
}

/**
 * Validates input for creating a new enhanced experiment
 */
export function validateCreateEnhancedExperimentInput(input: CreateEnhancedExperimentInput): void {
  if (!input.name || input.name.trim() === "") {
    throw new ConvexError("Experiment name is required");
  }
  
  if (input.temperature !== undefined && !input.temperatureUnit) {
    throw new ConvexError("Temperature value requires temperature unit");
  }
  
  if (input.coolingRate !== undefined && !input.coolingRateUnit) {
    throw new ConvexError("Cooling rate value requires cooling rate unit");
  }
  
  if (input.thawingRate !== undefined && !input.thawingRateUnit) {
    throw new ConvexError("Thawing rate value requires thawing rate unit");
  }
  
  if (input.pressure !== undefined && !input.pressureUnit) {
    throw new ConvexError("Pressure value requires pressure unit");
  }
  
  const validStatuses = ["planned", "in-progress", "completed", "failed"];
  if (input.status && !validStatuses.includes(input.status)) {
    throw new ConvexError(`Status must be one of: ${validStatuses.join(", ")}`);
  }
  
  if (input.version !== undefined && input.version < 1) {
    throw new ConvexError("Experiment version must be a positive number");
  }
}

/**
 * Validates input for updating an enhanced experiment
 */
export function validateUpdateEnhancedExperimentInput(input: UpdateEnhancedExperimentInput): void {
  if (input.name !== undefined && input.name.trim() === "") {
    throw new ConvexError("Experiment name cannot be empty");
  }
  
  if (input.temperature !== undefined && !input.temperatureUnit) {
    throw new ConvexError("Temperature value requires temperature unit");
  }
  
  if (input.coolingRate !== undefined && !input.coolingRateUnit) {
    throw new ConvexError("Cooling rate value requires cooling rate unit");
  }
  
  if (input.thawingRate !== undefined && !input.thawingRateUnit) {
    throw new ConvexError("Thawing rate value requires thawing rate unit");
  }
  
  if (input.pressure !== undefined && !input.pressureUnit) {
    throw new ConvexError("Pressure value requires pressure unit");
  }
  
  const validStatuses = ["planned", "in-progress", "completed", "failed"];
  if (input.status && !validStatuses.includes(input.status)) {
    throw new ConvexError(`Status must be one of: ${validStatuses.join(", ")}`);
  }
  
  if (input.version !== undefined && input.version < 1) {
    throw new ConvexError("Experiment version must be a positive number");
  }
}

/**
 * Validates input for creating a new enhanced experiment result
 */
export function validateCreateEnhancedExperimentResultInput(input: CreateEnhancedExperimentResultInput): void {
  if (!input.parameterName || input.parameterName.trim() === "") {
    throw new ConvexError("Parameter name is required");
  }
  
  if (input.value === undefined || input.value === null) {
    throw new ConvexError("Result value is required");
  }
  
  if (input.uncertainty) {
    validateUncertainty(input.uncertainty);
  }
  
  // Validate check for molecule OR mixture, but not both
  if (input.moleculeId && input.mixtureId) {
    throw new ConvexError("Result can only be associated with either a molecule OR a mixture, not both");
  }
}

/**
 * Validates input for updating an enhanced experiment result
 */
export function validateUpdateEnhancedExperimentResultInput(input: UpdateEnhancedExperimentResultInput): void {
  if (input.parameterName !== undefined && input.parameterName.trim() === "") {
    throw new ConvexError("Parameter name cannot be empty");
  }
  
  if (input.uncertainty) {
    validateUncertainty(input.uncertainty);
  }
}

/**
 * Validates input for creating a time series
 */
export function validateCreateTimeSeriesInput(input: CreateTimeSeriesInput): void {
  if (!input.name || input.name.trim() === "") {
    throw new ConvexError("Time series name is required");
  }
  
  if (!input.parameterName || input.parameterName.trim() === "") {
    throw new ConvexError("Parameter name is required");
  }
}

/**
 * Validates input for creating a time series data point
 */
export function validateCreateTimeSeriesDataPointInput(input: CreateTimeSeriesDataPointInput): void {
  if (input.timestamp === undefined) {
    throw new ConvexError("Timestamp is required");
  }
  
  if (input.value === undefined) {
    throw new ConvexError("Value is required");
  }
  
  if (input.uncertainty !== undefined && input.uncertainty < 0) {
    throw new ConvexError("Uncertainty cannot be negative");
  }
}

/**
 * Validates input for creating equipment
 */
export function validateCreateEquipmentInput(input: CreateEquipmentInput): void {
  if (!input.name || input.name.trim() === "") {
    throw new ConvexError("Equipment name is required");
  }
}

/**
 * Validates input for linking equipment to an experiment
 */
export function validateCreateExperimentEquipmentInput(input: CreateExperimentEquipmentInput): void {
  // All required fields are handled by type checking
}

/**
 * Validates input for creating a validation rule
 */
export function validateCreateValidationRuleInput(input: CreateValidationRuleInput): void {
  if (!input.parameterName || input.parameterName.trim() === "") {
    throw new ConvexError("Parameter name is required");
  }
  
  const validRuleTypes = ["range", "pattern", "comparison", "custom"];
  if (!validRuleTypes.includes(input.ruleType)) {
    throw new ConvexError(`Rule type must be one of: ${validRuleTypes.join(", ")}`);
  }
  
  if (!input.parameters || Object.keys(input.parameters).length === 0) {
    throw new ConvexError("Rule parameters are required");
  }
  
  const validSeverities = ["error", "warning", "info"];
  if (!validSeverities.includes(input.severity)) {
    throw new ConvexError(`Severity must be one of: ${validSeverities.join(", ")}`);
  }
  
  // Validate parameters based on rule type
  if (input.ruleType === "range") {
    if (input.parameters.min === undefined || input.parameters.max === undefined) {
      throw new ConvexError("Range rule must have min and max parameters");
    }
    
    if (input.parameters.min >= input.parameters.max) {
      throw new ConvexError("Range min must be less than max");
    }
  } else if (input.ruleType === "pattern") {
    if (!input.parameters.pattern) {
      throw new ConvexError("Pattern rule must have pattern parameter");
    }
  } else if (input.ruleType === "comparison") {
    if (!input.parameters.operator || !input.parameters.value) {
      throw new ConvexError("Comparison rule must have operator and value parameters");
    }
    
    const validOperators = ["eq", "ne", "gt", "lt", "gte", "lte"];
    if (!validOperators.includes(input.parameters.operator)) {
      throw new ConvexError(`Comparison operator must be one of: ${validOperators.join(", ")}`);
    }
  }
}

/**
 * Checks if a protocol exists
 */
export async function validateProtocolExists(db: DatabaseReader, protocolId: Id<"protocols">): Promise<void> {
  const protocol = await db.get(protocolId);
  if (!protocol) {
    throw new ConvexError(`Protocol with ID ${protocolId} not found`);
  }
}

/**
 * Checks if a tissue type exists
 */
export async function validateTissueTypeExists(db: DatabaseReader, tissueTypeId: Id<"tissueTypes">): Promise<void> {
  const tissueType = await db.get(tissueTypeId);
  if (!tissueType) {
    throw new ConvexError(`Tissue type with ID ${tissueTypeId} not found`);
  }
}

/**
 * Checks if an enhanced experiment exists
 */
export async function validateEnhancedExperimentExists(db: DatabaseReader, experimentId: Id<"enhancedExperiments">): Promise<void> {
  const experiment = await db.get(experimentId);
  if (!experiment) {
    throw new ConvexError(`Experiment with ID ${experimentId} not found`);
  }
}

/**
 * Checks if user has access to an enhanced experiment
 */
export async function validateEnhancedExperimentAccess(
  db: DatabaseReader,
  experimentId: Id<"enhancedExperiments">,
  userId: Id<"users"> | undefined,
  isPublic: boolean,
  conductedBy: Id<"users"> | undefined
): Promise<void> {
  // Public experiments are accessible to all
  if (isPublic) {
    return;
  }
  
  // If not public, only the creator or project members can access
  if (!userId) {
    throw new ConvexError("Authentication required to access this experiment");
  }
  
  // Check if user is the creator
  if (conductedBy && userId.equals(conductedBy)) {
    return;
  }
  
  // Check if experiment has a project and user is a member
  const experiment = await db.get(experimentId);
  if (!experiment) {
    throw new ConvexError(`Experiment with ID ${experimentId} not found`);
  }
  
  if (experiment.projectId) {
    const teamMember = await db
      .query("teamMembers")
      .withIndex("by_user", q => q.eq("userId", userId))
      .filter(q => q.eq(q.field("projectId"), experiment.projectId!))
      .first();
    
    if (teamMember) {
      return;
    }
  }
  
  throw new ConvexError("You do not have permission to access this experiment");
}
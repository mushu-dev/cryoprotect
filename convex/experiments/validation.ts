/**
 * Validation functions for experiments and experiment results
 */

import { ConvexError } from "convex/values";
import { Id } from "../_generated/dataModel";
import { 
  CreateExperimentInput, 
  UpdateExperimentInput,
  CreateExperimentResultInput,
  UpdateExperimentResultInput 
} from "./types";

/**
 * Validates the input for creating an experiment
 */
export function validateCreateExperimentInput(input: CreateExperimentInput): void {
  // Name is required and should be a non-empty string
  if (!input.name || input.name.trim() === "") {
    throw new ConvexError("Experiment name is required");
  }
  
  // Description should be a string if provided
  if (input.description !== undefined && typeof input.description !== "string") {
    throw new ConvexError("Experiment description must be a string");
  }
  
  // MixtureId should be a valid ID if provided
  if (input.mixtureId !== undefined && typeof input.mixtureId !== "object") {
    throw new ConvexError("Mixture ID must be a valid mixture reference");
  }
  
  // Protocol should be a string if provided
  if (input.protocol !== undefined && typeof input.protocol !== "string") {
    throw new ConvexError("Protocol must be a string");
  }
  
  // ProjectId should be a valid ID if provided
  if (input.projectId !== undefined && typeof input.projectId !== "object") {
    throw new ConvexError("Project ID must be a valid project reference");
  }
  
  // Date should be a valid timestamp if provided
  if (input.date !== undefined && typeof input.date !== "number") {
    throw new ConvexError("Date must be a valid timestamp");
  }
  
  // Status should be a valid value if provided
  if (input.status !== undefined) {
    const validStatuses = ["planned", "in-progress", "completed", "failed"];
    if (!validStatuses.includes(input.status)) {
      throw new ConvexError(
        `Status must be one of: ${validStatuses.join(", ")}`
      );
    }
  }
  
  // Public should be a boolean if provided
  if (input.public !== undefined && typeof input.public !== "boolean") {
    throw new ConvexError("Public flag must be a boolean");
  }
}

/**
 * Validates the input for updating an experiment
 */
export function validateUpdateExperimentInput(input: UpdateExperimentInput): void {
  // Name should be a non-empty string if provided
  if (input.name !== undefined && (typeof input.name !== "string" || input.name.trim() === "")) {
    throw new ConvexError("Experiment name must be a non-empty string");
  }
  
  // Description should be a string if provided
  if (input.description !== undefined && typeof input.description !== "string") {
    throw new ConvexError("Experiment description must be a string");
  }
  
  // MixtureId should be a valid ID if provided
  if (input.mixtureId !== undefined && typeof input.mixtureId !== "object") {
    throw new ConvexError("Mixture ID must be a valid mixture reference");
  }
  
  // Protocol should be a string if provided
  if (input.protocol !== undefined && typeof input.protocol !== "string") {
    throw new ConvexError("Protocol must be a string");
  }
  
  // ProjectId should be a valid ID if provided
  if (input.projectId !== undefined && typeof input.projectId !== "object") {
    throw new ConvexError("Project ID must be a valid project reference");
  }
  
  // Date should be a valid timestamp if provided
  if (input.date !== undefined && typeof input.date !== "number") {
    throw new ConvexError("Date must be a valid timestamp");
  }
  
  // Status should be a valid value if provided
  if (input.status !== undefined) {
    const validStatuses = ["planned", "in-progress", "completed", "failed"];
    if (!validStatuses.includes(input.status)) {
      throw new ConvexError(
        `Status must be one of: ${validStatuses.join(", ")}`
      );
    }
  }
  
  // Public should be a boolean if provided
  if (input.public !== undefined && typeof input.public !== "boolean") {
    throw new ConvexError("Public flag must be a boolean");
  }
}

/**
 * Validates the input for creating an experiment result
 */
export function validateCreateExperimentResultInput(input: CreateExperimentResultInput): void {
  // ExperimentId is required
  if (!input.experimentId) {
    throw new ConvexError("Experiment ID is required");
  }
  
  // ParameterName is required and should be a non-empty string
  if (!input.parameterName || typeof input.parameterName !== "string" || input.parameterName.trim() === "") {
    throw new ConvexError("Parameter name is required");
  }
  
  // Value type should be valid
  if (input.value !== null && 
      typeof input.value !== "string" && 
      typeof input.value !== "number" && 
      typeof input.value !== "boolean") {
    throw new ConvexError("Value must be a string, number, boolean, or null");
  }
  
  // NumericValue should be a number if provided
  if (input.numericValue !== undefined && typeof input.numericValue !== "number") {
    throw new ConvexError("Numeric value must be a number");
  }
  
  // Units should be a string if provided
  if (input.units !== undefined && typeof input.units !== "string") {
    throw new ConvexError("Units must be a string");
  }
  
  // Notes should be a string if provided
  if (input.notes !== undefined && typeof input.notes !== "string") {
    throw new ConvexError("Notes must be a string");
  }
  
  // If value is numeric, numericValue should be provided for efficient querying
  if (typeof input.value === "number" && input.numericValue === undefined) {
    throw new ConvexError("numericValue must be provided when value is a number");
  }
}

/**
 * Validates the input for updating an experiment result
 */
export function validateUpdateExperimentResultInput(input: UpdateExperimentResultInput): void {
  // ParameterName should be a non-empty string if provided
  if (input.parameterName !== undefined && 
      (typeof input.parameterName !== "string" || input.parameterName.trim() === "")) {
    throw new ConvexError("Parameter name must be a non-empty string");
  }
  
  // Value type should be valid if provided
  if (input.value !== undefined && 
      input.value !== null && 
      typeof input.value !== "string" && 
      typeof input.value !== "number" && 
      typeof input.value !== "boolean") {
    throw new ConvexError("Value must be a string, number, boolean, or null");
  }
  
  // NumericValue should be a number if provided
  if (input.numericValue !== undefined && typeof input.numericValue !== "number") {
    throw new ConvexError("Numeric value must be a number");
  }
  
  // Units should be a string if provided
  if (input.units !== undefined && typeof input.units !== "string") {
    throw new ConvexError("Units must be a string");
  }
  
  // Notes should be a string if provided
  if (input.notes !== undefined && typeof input.notes !== "string") {
    throw new ConvexError("Notes must be a string");
  }
  
  // If value is numeric, numericValue should be provided for efficient querying
  if (typeof input.value === "number" && input.numericValue === undefined) {
    throw new ConvexError("numericValue must be provided when value is a number");
  }
}

/**
 * Validates that an experiment ID exists
 */
export function validateExperimentExists(experimentId: Id<"experiments"> | undefined): void {
  if (!experimentId) {
    throw new ConvexError("Experiment ID is required");
  }
}

/**
 * Validates that an experiment result ID exists
 */
export function validateExperimentResultExists(resultId: Id<"experimentResults"> | undefined): void {
  if (!resultId) {
    throw new ConvexError("Experiment result ID is required");
  }
}

/**
 * Validates that a user has access to an experiment
 */
export function validateExperimentAccess(
  experimentId: Id<"experiments">, 
  userId: Id<"users"> | null, 
  isPublic: boolean, 
  conductedById?: Id<"users">
): void {
  // If experiment is public, anyone can access
  if (isPublic) {
    return;
  }
  
  // If user is not logged in and experiment is not public, deny access
  if (!userId) {
    throw new ConvexError("Authentication required to access this experiment");
  }
  
  // If experiment belongs to user, allow access
  if (conductedById && userId.equals(conductedById)) {
    return;
  }
  
  // Otherwise deny access
  // Note: In a real implementation, we would check project membership here
  throw new ConvexError("You do not have permission to access this experiment");
}
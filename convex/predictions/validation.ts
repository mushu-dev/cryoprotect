/**
 * Validation functions for predictions and scientific models
 */

import { ConvexError } from "convex/values";
import { Id } from "../_generated/dataModel";
import { 
  CreateScientificModelInput, 
  UpdateScientificModelInput,
  CreatePredictionInput,
  UpdatePredictionInput 
} from "./types";

/**
 * Validates the input for creating a scientific model
 */
export function validateCreateScientificModelInput(input: CreateScientificModelInput): void {
  // Name is required and should be a non-empty string
  if (!input.name || input.name.trim() === "") {
    throw new ConvexError("Model name is required");
  }
  
  // Description should be a string if provided
  if (input.description !== undefined && typeof input.description !== "string") {
    throw new ConvexError("Model description must be a string");
  }
  
  // Version is required and should be a non-empty string
  if (!input.version || typeof input.version !== "string" || input.version.trim() === "") {
    throw new ConvexError("Model version is required");
  }
  
  // Type is required and should be a non-empty string
  if (!input.type || typeof input.type !== "string" || input.type.trim() === "") {
    throw new ConvexError("Model type is required");
  }
  
  // Parameters should be a string if provided
  if (input.parameters !== undefined && typeof input.parameters !== "string") {
    throw new ConvexError("Model parameters must be a string");
  }
  
  // Public should be a boolean if provided
  if (input.public !== undefined && typeof input.public !== "boolean") {
    throw new ConvexError("Public flag must be a boolean");
  }
}

/**
 * Validates the input for updating a scientific model
 */
export function validateUpdateScientificModelInput(input: UpdateScientificModelInput): void {
  // Name should be a non-empty string if provided
  if (input.name !== undefined && (typeof input.name !== "string" || input.name.trim() === "")) {
    throw new ConvexError("Model name must be a non-empty string");
  }
  
  // Description should be a string if provided
  if (input.description !== undefined && typeof input.description !== "string") {
    throw new ConvexError("Model description must be a string");
  }
  
  // Version should be a non-empty string if provided
  if (input.version !== undefined && (typeof input.version !== "string" || input.version.trim() === "")) {
    throw new ConvexError("Model version must be a non-empty string");
  }
  
  // Type should be a non-empty string if provided
  if (input.type !== undefined && (typeof input.type !== "string" || input.type.trim() === "")) {
    throw new ConvexError("Model type must be a non-empty string");
  }
  
  // Parameters should be a string if provided
  if (input.parameters !== undefined && typeof input.parameters !== "string") {
    throw new ConvexError("Model parameters must be a string");
  }
  
  // Public should be a boolean if provided
  if (input.public !== undefined && typeof input.public !== "boolean") {
    throw new ConvexError("Public flag must be a boolean");
  }
}

/**
 * Validates the input for creating a prediction
 */
export function validateCreatePredictionInput(input: CreatePredictionInput): void {
  // ModelId is required
  if (!input.modelId) {
    throw new ConvexError("Model ID is required");
  }
  
  // Either moleculeId or mixtureId must be provided, but not both
  if ((!input.moleculeId && !input.mixtureId) || (input.moleculeId && input.mixtureId)) {
    throw new ConvexError("Either molecule ID or mixture ID must be provided, but not both");
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
  
  // Confidence should be a number between 0 and 1 if provided
  if (input.confidence !== undefined) {
    if (typeof input.confidence !== "number") {
      throw new ConvexError("Confidence must be a number");
    }
    
    if (input.confidence < 0 || input.confidence > 1) {
      throw new ConvexError("Confidence must be between 0 and 1");
    }
  }
  
  // If value is numeric, numericValue should be provided for efficient querying
  if (typeof input.value === "number" && input.numericValue === undefined) {
    throw new ConvexError("numericValue must be provided when value is a number");
  }
}

/**
 * Validates the input for updating a prediction
 */
export function validateUpdatePredictionInput(input: UpdatePredictionInput): void {
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
  
  // Confidence should be a number between 0 and 1 if provided
  if (input.confidence !== undefined) {
    if (typeof input.confidence !== "number") {
      throw new ConvexError("Confidence must be a number");
    }
    
    if (input.confidence < 0 || input.confidence > 1) {
      throw new ConvexError("Confidence must be between 0 and 1");
    }
  }
  
  // If value is numeric, numericValue should be provided for efficient querying
  if (typeof input.value === "number" && input.numericValue === undefined) {
    throw new ConvexError("numericValue must be provided when value is a number");
  }
}

/**
 * Validates that a scientific model ID exists
 */
export function validateScientificModelExists(modelId: Id<"scientificModels"> | undefined): void {
  if (!modelId) {
    throw new ConvexError("Scientific model ID is required");
  }
}

/**
 * Validates that a prediction ID exists
 */
export function validatePredictionExists(predictionId: Id<"predictions"> | undefined): void {
  if (!predictionId) {
    throw new ConvexError("Prediction ID is required");
  }
}

/**
 * Validates that a user has access to a scientific model
 */
export function validateScientificModelAccess(
  modelId: Id<"scientificModels">, 
  userId: Id<"users"> | null, 
  isPublic: boolean, 
  createdById?: Id<"users">
): void {
  // If model is public, anyone can access
  if (isPublic) {
    return;
  }
  
  // If user is not logged in and model is not public, deny access
  if (!userId) {
    throw new ConvexError("Authentication required to access this scientific model");
  }
  
  // If model belongs to user, allow access
  if (createdById && userId.equals(createdById)) {
    return;
  }
  
  // Otherwise deny access
  throw new ConvexError("You do not have permission to access this scientific model");
}

/**
 * Validates that target (molecule or mixture) exists
 */
export function validatePredictionTarget(
  moleculeId?: Id<"molecules">,
  mixtureId?: Id<"mixtures">
): void {
  // Either moleculeId or mixtureId must be provided, but not both
  if ((!moleculeId && !mixtureId) || (moleculeId && mixtureId)) {
    throw new ConvexError("Either molecule ID or mixture ID must be provided, but not both");
  }
}
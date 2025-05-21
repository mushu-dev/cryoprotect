/**
 * Validation functions for mixtures and mixture components
 */

import { ConvexError } from "convex/values";
import { Id } from "../_generated/dataModel";
import { 
  CreateMixtureInput, 
  UpdateMixtureInput,
  CreateMixtureComponentInput,
  UpdateMixtureComponentInput 
} from "./types";

/**
 * Validates the input for creating a mixture
 */
export function validateCreateMixtureInput(input: CreateMixtureInput): void {
  // Name is required and should be a non-empty string
  if (!input.name || input.name.trim() === "") {
    throw new ConvexError("Mixture name is required");
  }
  
  // Type should be a valid string if provided
  if (input.type !== undefined && typeof input.type !== "string") {
    throw new ConvexError("Mixture type must be a string");
  }
  
  // Description should be a string if provided
  if (input.description !== undefined && typeof input.description !== "string") {
    throw new ConvexError("Mixture description must be a string");
  }
  
  // ProjectId should be a valid ID if provided
  if (input.projectId !== undefined && typeof input.projectId !== "object") {
    throw new ConvexError("Project ID must be a valid project reference");
  }
  
  // Public should be a boolean if provided
  if (input.public !== undefined && typeof input.public !== "boolean") {
    throw new ConvexError("Public flag must be a boolean");
  }
}

/**
 * Validates the input for updating a mixture
 */
export function validateUpdateMixtureInput(input: UpdateMixtureInput): void {
  // Name should be a non-empty string if provided
  if (input.name !== undefined && (typeof input.name !== "string" || input.name.trim() === "")) {
    throw new ConvexError("Mixture name must be a non-empty string");
  }
  
  // Type should be a valid string if provided
  if (input.type !== undefined && typeof input.type !== "string") {
    throw new ConvexError("Mixture type must be a string");
  }
  
  // Description should be a string if provided
  if (input.description !== undefined && typeof input.description !== "string") {
    throw new ConvexError("Mixture description must be a string");
  }
  
  // ProjectId should be a valid ID if provided
  if (input.projectId !== undefined && typeof input.projectId !== "object") {
    throw new ConvexError("Project ID must be a valid project reference");
  }
  
  // Public should be a boolean if provided
  if (input.public !== undefined && typeof input.public !== "boolean") {
    throw new ConvexError("Public flag must be a boolean");
  }
}

/**
 * Validates the input for creating a mixture component
 */
export function validateCreateMixtureComponentInput(input: CreateMixtureComponentInput): void {
  // MixtureId is required
  if (!input.mixtureId) {
    throw new ConvexError("Mixture ID is required");
  }
  
  // MoleculeId is required
  if (!input.moleculeId) {
    throw new ConvexError("Molecule ID is required");
  }
  
  // Concentration is required and should be a positive number
  if (typeof input.concentration !== "number" || input.concentration < 0) {
    throw new ConvexError("Concentration must be a positive number");
  }
  
  // Units is required and should be a non-empty string
  if (!input.units || typeof input.units !== "string" || input.units.trim() === "") {
    throw new ConvexError("Units are required");
  }
  
  // Role should be a string if provided
  if (input.role !== undefined && typeof input.role !== "string") {
    throw new ConvexError("Role must be a string");
  }
  
  // Notes should be a string if provided
  if (input.notes !== undefined && typeof input.notes !== "string") {
    throw new ConvexError("Notes must be a string");
  }
}

/**
 * Validates the input for updating a mixture component
 */
export function validateUpdateMixtureComponentInput(input: UpdateMixtureComponentInput): void {
  // Concentration should be a positive number if provided
  if (input.concentration !== undefined && (typeof input.concentration !== "number" || input.concentration < 0)) {
    throw new ConvexError("Concentration must be a positive number");
  }
  
  // Units should be a non-empty string if provided
  if (input.units !== undefined && (typeof input.units !== "string" || input.units.trim() === "")) {
    throw new ConvexError("Units must be a non-empty string");
  }
  
  // Role should be a string if provided
  if (input.role !== undefined && typeof input.role !== "string") {
    throw new ConvexError("Role must be a string");
  }
  
  // Notes should be a string if provided
  if (input.notes !== undefined && typeof input.notes !== "string") {
    throw new ConvexError("Notes must be a string");
  }
}

/**
 * Validates that a mixture ID exists
 */
export function validateMixtureExists(mixtureId: Id<"mixtures"> | undefined): void {
  if (!mixtureId) {
    throw new ConvexError("Mixture ID is required");
  }
}

/**
 * Validates that a mixture component ID exists
 */
export function validateMixtureComponentExists(componentId: Id<"mixtureComponents"> | undefined): void {
  if (!componentId) {
    throw new ConvexError("Mixture component ID is required");
  }
}

/**
 * Validates that a user has access to a mixture
 */
export function validateMixtureAccess(
  mixtureId: Id<"mixtures">, 
  userId: Id<"users"> | null, 
  isPublic: boolean, 
  createdById?: Id<"users">
): void {
  // If mixture is public, anyone can access
  if (isPublic) {
    return;
  }
  
  // If user is not logged in and mixture is not public, deny access
  if (!userId) {
    throw new ConvexError("Authentication required to access this mixture");
  }
  
  // If mixture belongs to user, allow access
  if (createdById && userId.equals(createdById)) {
    return;
  }
  
  // Otherwise deny access
  // Note: In a real implementation, we would check project membership here
  throw new ConvexError("You do not have permission to access this mixture");
}
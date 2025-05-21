/**
 * Validation functions for molecular cross references and synonyms
 */

import { ConvexError } from "convex/values";
import { Id } from "../_generated/dataModel";
import { 
  CreateCrossReferenceInput, 
  UpdateCrossReferenceInput,
  CreateSynonymInput,
  UpdateSynonymInput
} from "./types";

/**
 * Validates the input for creating a cross reference
 */
export function validateCreateCrossReferenceInput(input: CreateCrossReferenceInput): void {
  // MoleculeId is required
  if (!input.moleculeId) {
    throw new ConvexError("Molecule ID is required");
  }
  
  // DatabaseName is required and should be a non-empty string
  if (!input.databaseName || input.databaseName.trim() === "") {
    throw new ConvexError("Database name is required");
  }
  
  // Identifier is required and should be a non-empty string
  if (!input.identifier || input.identifier.trim() === "") {
    throw new ConvexError("Identifier is required");
  }
  
  // URL should be a string if provided
  if (input.url !== undefined && typeof input.url !== "string") {
    throw new ConvexError("URL must be a string");
  }
  
  // Validate URL format if provided
  if (input.url && !isValidUrl(input.url)) {
    throw new ConvexError("URL must be valid");
  }
}

/**
 * Validates the input for updating a cross reference
 */
export function validateUpdateCrossReferenceInput(input: UpdateCrossReferenceInput): void {
  // DatabaseName should be a non-empty string if provided
  if (input.databaseName !== undefined && 
      (typeof input.databaseName !== "string" || input.databaseName.trim() === "")) {
    throw new ConvexError("Database name must be a non-empty string");
  }
  
  // Identifier should be a non-empty string if provided
  if (input.identifier !== undefined && 
      (typeof input.identifier !== "string" || input.identifier.trim() === "")) {
    throw new ConvexError("Identifier must be a non-empty string");
  }
  
  // URL should be a string if provided
  if (input.url !== undefined && typeof input.url !== "string") {
    throw new ConvexError("URL must be a string");
  }
  
  // Validate URL format if provided
  if (input.url && !isValidUrl(input.url)) {
    throw new ConvexError("URL must be valid");
  }
}

/**
 * Validates the input for creating a synonym
 */
export function validateCreateSynonymInput(input: CreateSynonymInput): void {
  // MoleculeId is required
  if (!input.moleculeId) {
    throw new ConvexError("Molecule ID is required");
  }
  
  // Name is required and should be a non-empty string
  if (!input.name || input.name.trim() === "") {
    throw new ConvexError("Synonym name is required");
  }
  
  // Type should be a string if provided
  if (input.type !== undefined && typeof input.type !== "string") {
    throw new ConvexError("Type must be a string");
  }
  
  // Source should be a valid ID if provided
  if (input.source !== undefined && typeof input.source !== "object") {
    throw new ConvexError("Source must be a valid data source reference");
  }
}

/**
 * Validates the input for updating a synonym
 */
export function validateUpdateSynonymInput(input: UpdateSynonymInput): void {
  // Name should be a non-empty string if provided
  if (input.name !== undefined && 
      (typeof input.name !== "string" || input.name.trim() === "")) {
    throw new ConvexError("Name must be a non-empty string");
  }
  
  // Type should be a string if provided
  if (input.type !== undefined && typeof input.type !== "string") {
    throw new ConvexError("Type must be a string");
  }
  
  // Source should be a valid ID if provided
  if (input.source !== undefined && typeof input.source !== "object") {
    throw new ConvexError("Source must be a valid data source reference");
  }
}

/**
 * Validates that a cross reference ID exists
 */
export function validateCrossReferenceExists(crossRefId: Id<"moleculeCrossReferences"> | undefined): void {
  if (!crossRefId) {
    throw new ConvexError("Cross reference ID is required");
  }
}

/**
 * Validates that a synonym ID exists
 */
export function validateSynonymExists(synonymId: Id<"moleculeSynonyms"> | undefined): void {
  if (!synonymId) {
    throw new ConvexError("Synonym ID is required");
  }
}

/**
 * Validates that a molecule ID exists
 */
export function validateMoleculeExists(moleculeId: Id<"molecules"> | undefined): void {
  if (!moleculeId) {
    throw new ConvexError("Molecule ID is required");
  }
}

/**
 * Helper to validate URL format
 */
function isValidUrl(url: string): boolean {
  try {
    // Simple URL validation
    new URL(url);
    return true;
  } catch (error) {
    return false;
  }
}
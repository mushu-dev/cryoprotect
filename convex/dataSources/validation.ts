/**
 * Validation functions for data sources
 */

import { ConvexError } from "convex/values";
import { Id } from "../_generated/dataModel";
import { 
  CreateDataSourceInput, 
  UpdateDataSourceInput 
} from "./types";

/**
 * Validates the input for creating a data source
 */
export function validateCreateDataSourceInput(input: CreateDataSourceInput): void {
  // Name is required and should be a non-empty string
  if (!input.name || input.name.trim() === "") {
    throw new ConvexError("Data source name is required");
  }
  
  // Description should be a string if provided
  if (input.description !== undefined && typeof input.description !== "string") {
    throw new ConvexError("Description must be a string");
  }
  
  // URL should be a string if provided
  if (input.url !== undefined && typeof input.url !== "string") {
    throw new ConvexError("URL must be a string");
  }
  
  // Validate URL format if provided
  if (input.url && !isValidUrl(input.url)) {
    throw new ConvexError("URL must be valid");
  }
  
  // Type is required and should be a valid value
  const validTypes = ["database", "publication", "experiment", "calculation"];
  if (!input.type || !validTypes.includes(input.type)) {
    throw new ConvexError(`Type must be one of: ${validTypes.join(", ")}`);
  }
  
  // Version should be a string if provided
  if (input.version !== undefined && typeof input.version !== "string") {
    throw new ConvexError("Version must be a string");
  }
}

/**
 * Validates the input for updating a data source
 */
export function validateUpdateDataSourceInput(input: UpdateDataSourceInput): void {
  // Name should be a non-empty string if provided
  if (input.name !== undefined && (typeof input.name !== "string" || input.name.trim() === "")) {
    throw new ConvexError("Name must be a non-empty string");
  }
  
  // Description should be a string if provided
  if (input.description !== undefined && typeof input.description !== "string") {
    throw new ConvexError("Description must be a string");
  }
  
  // URL should be a string if provided
  if (input.url !== undefined && typeof input.url !== "string") {
    throw new ConvexError("URL must be a string");
  }
  
  // Validate URL format if provided
  if (input.url && !isValidUrl(input.url)) {
    throw new ConvexError("URL must be valid");
  }
  
  // Type should be a valid value if provided
  if (input.type !== undefined) {
    const validTypes = ["database", "publication", "experiment", "calculation"];
    if (!validTypes.includes(input.type)) {
      throw new ConvexError(`Type must be one of: ${validTypes.join(", ")}`);
    }
  }
  
  // Version should be a string if provided
  if (input.version !== undefined && typeof input.version !== "string") {
    throw new ConvexError("Version must be a string");
  }
}

/**
 * Validates that a data source ID exists
 */
export function validateDataSourceExists(dataSourceId: Id<"dataSources"> | undefined): void {
  if (!dataSourceId) {
    throw new ConvexError("Data source ID is required");
  }
}

/**
 * Validates that a data source can be deleted (no references)
 */
export function validateDataSourceNotReferenced(
  dataSourceId: Id<"dataSources">,
  moleculeCount: number,
  propertyCount: number,
  experimentCount: number
): void {
  const totalReferences = moleculeCount + propertyCount + experimentCount;
  
  if (totalReferences > 0) {
    throw new ConvexError(
      `Cannot delete data source that is referenced by ${totalReferences} entities. ` +
      `(${moleculeCount} molecules, ${propertyCount} properties, ${experimentCount} experiments)`
    );
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
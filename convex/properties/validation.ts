/**
 * Validation functions for property data
 */

import { CreatePropertyTypeInput, UpdatePropertyTypeInput, CreateMolecularPropertyInput, UpdateMolecularPropertyInput } from "./types";

/**
 * Validate a property type name
 * @param name - The property type name to validate
 * @returns True if the name is valid
 */
export const isValidPropertyTypeName = (name: string): boolean => {
  return /^[a-zA-Z0-9_]+$/.test(name) && name.trim().length > 0;
};

/**
 * Validate a property type display name
 * @param displayName - The display name to validate
 * @returns True if the display name is valid
 */
export const isValidDisplayName = (displayName: string): boolean => {
  return displayName.trim().length > 0;
};

/**
 * Validate a property data type
 * @param dataType - The data type to validate
 * @returns True if the data type is valid
 */
export const isValidDataType = (dataType: string): boolean => {
  return ["string", "number", "boolean"].includes(dataType);
};

/**
 * Validate a property value based on its data type
 * @param value - The value to validate
 * @param dataType - The expected data type
 * @returns True if the value matches the data type
 */
export const isValidPropertyValue = (
  value: any,
  dataType: "string" | "number" | "boolean"
): boolean => {
  if (value === null) {
    return true; // Null is allowed for all types
  }

  switch (dataType) {
    case "string":
      return typeof value === "string";
    case "number":
      return typeof value === "number" && !isNaN(value);
    case "boolean":
      return typeof value === "boolean";
    default:
      return false;
  }
};

/**
 * Validate a property type creation input
 * @param input - The property type creation input to validate
 * @returns An object with validation result and optional error message
 */
export const validateCreatePropertyTypeInput = (
  input: CreatePropertyTypeInput
): { valid: boolean; error?: string } => {
  if (!isValidPropertyTypeName(input.name)) {
    return { valid: false, error: "Invalid property type name" };
  }

  if (!isValidDisplayName(input.displayName)) {
    return { valid: false, error: "Invalid display name" };
  }

  if (!isValidDataType(input.dataType)) {
    return { valid: false, error: "Invalid data type" };
  }

  return { valid: true };
};

/**
 * Validate a property type update input
 * @param input - The property type update input to validate
 * @returns An object with validation result and optional error message
 */
export const validateUpdatePropertyTypeInput = (
  input: UpdatePropertyTypeInput
): { valid: boolean; error?: string } => {
  if (input.displayName && !isValidDisplayName(input.displayName)) {
    return { valid: false, error: "Invalid display name" };
  }

  return { valid: true };
};

/**
 * Validate a molecular property creation input
 * @param input - The molecular property creation input to validate
 * @param propertyDataType - The data type of the property
 * @returns An object with validation result and optional error message
 */
export const validateCreateMolecularPropertyInput = (
  input: CreateMolecularPropertyInput,
  propertyDataType: "string" | "number" | "boolean"
): { valid: boolean; error?: string } => {
  if (!isValidPropertyValue(input.value, propertyDataType)) {
    return { 
      valid: false, 
      error: `Invalid value for property with data type ${propertyDataType}` 
    };
  }

  if (
    propertyDataType === "number" &&
    input.numericValue !== undefined &&
    typeof input.numericValue !== "number"
  ) {
    return { valid: false, error: "numericValue must be a number" };
  }

  if (
    input.confidence !== undefined &&
    (typeof input.confidence !== "number" ||
      input.confidence < 0 ||
      input.confidence > 1)
  ) {
    return {
      valid: false,
      error: "Confidence must be a number between 0 and 1",
    };
  }

  return { valid: true };
};

/**
 * Validate a molecular property update input
 * @param input - The molecular property update input to validate
 * @param propertyDataType - The data type of the property
 * @returns An object with validation result and optional error message
 */
export const validateUpdateMolecularPropertyInput = (
  input: UpdateMolecularPropertyInput,
  propertyDataType: "string" | "number" | "boolean"
): { valid: boolean; error?: string } => {
  if (
    input.value !== undefined &&
    !isValidPropertyValue(input.value, propertyDataType)
  ) {
    return {
      valid: false,
      error: `Invalid value for property with data type ${propertyDataType}`,
    };
  }

  if (
    propertyDataType === "number" &&
    input.numericValue !== undefined &&
    typeof input.numericValue !== "number"
  ) {
    return { valid: false, error: "numericValue must be a number" };
  }

  if (
    input.confidence !== undefined &&
    (typeof input.confidence !== "number" ||
      input.confidence < 0 ||
      input.confidence > 1)
  ) {
    return {
      valid: false,
      error: "Confidence must be a number between 0 and 1",
    };
  }

  return { valid: true };
};
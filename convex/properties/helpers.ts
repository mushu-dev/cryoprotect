/**
 * Helper functions for working with property data
 */

import { DatabaseReader } from "../_generated/server";
import { Id } from "../_generated/dataModel";
import { PropertyType, MolecularProperty } from "./types";
import { getCurrentTimestamp } from "../utils/common";

/**
 * Get a property type by ID
 * @param db - Convex database reader
 * @param propertyTypeId - The ID of the property type to retrieve
 * @returns The property type or null if not found
 */
export const getPropertyTypeById = async (
  db: DatabaseReader,
  propertyTypeId: Id<"propertyTypes">
): Promise<PropertyType | null> => {
  return await db.get(propertyTypeId) as PropertyType | null;
};

/**
 * Get a property type by name
 * @param db - Convex database reader
 * @param name - The name of the property type to retrieve
 * @returns The property type or null if not found
 */
export const getPropertyTypeByName = async (
  db: DatabaseReader,
  name: string
): Promise<PropertyType | null> => {
  return await db
    .query("propertyTypes")
    .withIndex("by_name", (q) => q.eq("name", name))
    .first() as PropertyType | null;
};

/**
 * Get a molecular property by ID
 * @param db - Convex database reader
 * @param propertyId - The ID of the molecular property to retrieve
 * @returns The molecular property or null if not found
 */
export const getMolecularPropertyById = async (
  db: DatabaseReader,
  propertyId: Id<"molecularProperties">
): Promise<MolecularProperty | null> => {
  return await db.get(propertyId) as MolecularProperty | null;
};

/**
 * Check if a molecular property with the given molecule ID and property type ID exists
 * @param db - Convex database reader
 * @param moleculeId - The ID of the molecule
 * @param propertyTypeId - The ID of the property type
 * @returns True if a property exists for the molecule and property type
 */
export const moleculePropertyExists = async (
  db: DatabaseReader,
  moleculeId: Id<"molecules">,
  propertyTypeId: Id<"propertyTypes">
): Promise<boolean> => {
  const property = await db
    .query("molecularProperties")
    .withIndex("by_molecule_property", (q) => 
      q.eq("moleculeId", moleculeId).eq("propertyTypeId", propertyTypeId)
    )
    .first();
  
  return !!property;
};

/**
 * Convert a value to the specified data type
 * @param value - The value to convert
 * @param dataType - The target data type
 * @returns The converted value
 */
export const convertValueToType = (
  value: any,
  dataType: "string" | "number" | "boolean"
): string | number | boolean | null => {
  if (value === null || value === undefined) {
    return null;
  }

  switch (dataType) {
    case "string":
      return String(value);
    case "number":
      const num = Number(value);
      return isNaN(num) ? null : num;
    case "boolean":
      if (typeof value === "string") {
        return value.toLowerCase() === "true" || value === "1";
      }
      return Boolean(value);
    default:
      return null;
  }
};

/**
 * Extract numeric value from a property value
 * @param value - The property value
 * @returns The numeric value or undefined if not applicable
 */
export const extractNumericValue = (value: any): number | undefined => {
  if (typeof value === "number" && !isNaN(value)) {
    return value;
  }
  
  if (typeof value === "string") {
    // Try to extract number from string (remove units, etc.)
    const numberMatch = value.match(/^[-+]?[0-9]*\.?[0-9]+/);
    if (numberMatch) {
      const num = parseFloat(numberMatch[0]);
      if (!isNaN(num)) {
        return num;
      }
    }
  }
  
  return undefined;
};

/**
 * Get all properties for a molecule
 * @param db - Convex database reader
 * @param moleculeId - The ID of the molecule
 * @returns Array of properties for the molecule
 */
export const getMoleculeProperties = async (
  db: DatabaseReader,
  moleculeId: Id<"molecules">
): Promise<Array<MolecularProperty & { propertyType: PropertyType | null }>> => {
  const properties = await db
    .query("molecularProperties")
    .withIndex("by_molecule", (q) => q.eq("moleculeId", moleculeId))
    .collect() as MolecularProperty[];
  
  // Get all property types in a single batch
  const propertyTypeIds = [...new Set(properties.map(p => p.propertyTypeId))];
  const propertyTypes = await Promise.all(
    propertyTypeIds.map(id => db.get(id))
  );
  
  // Create a lookup map for property types
  const propertyTypeMap = new Map();
  propertyTypes.forEach(type => {
    if (type) {
      propertyTypeMap.set(type._id.toString(), type);
    }
  });
  
  // Join property types with properties
  return properties.map(property => ({
    ...property,
    propertyType: propertyTypeMap.get(property.propertyTypeId.toString()) || null
  }));
};
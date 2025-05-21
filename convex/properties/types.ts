/**
 * Type definitions for molecular properties
 */

import { Id } from "../_generated/dataModel";

/**
 * Property type document in the database
 */
export interface PropertyType {
  _id: Id<"propertyTypes">;
  _creationTime: number;
  
  // Core fields
  name: string;
  displayName: string;
  description?: string;
  dataType: "string" | "number" | "boolean";
  units?: string;
  defaultUnits?: string;
  category?: string;
  isCalculated?: boolean;
  
  // Metadata
  createdAt: number;
  updatedAt: number;
}

/**
 * Input for creating a new property type
 */
export interface CreatePropertyTypeInput {
  name: string;
  displayName: string;
  description?: string;
  dataType: "string" | "number" | "boolean";
  units?: string;
  defaultUnits?: string;
  category?: string;
  isCalculated?: boolean;
}

/**
 * Input for updating a property type
 */
export interface UpdatePropertyTypeInput {
  displayName?: string;
  description?: string;
  units?: string;
  defaultUnits?: string;
  category?: string;
  isCalculated?: boolean;
}

/**
 * Molecular property document in the database
 */
export interface MolecularProperty {
  _id: Id<"molecularProperties">;
  _creationTime: number;
  
  // Core fields
  moleculeId: Id<"molecules">;
  propertyTypeId: Id<"propertyTypes">;
  value: string | number | boolean | null;
  numericValue?: number;
  units?: string;
  source?: Id<"dataSources">;
  calculationMethod?: string;
  confidence?: number;
  
  // Metadata
  createdAt: number;
  updatedAt: number;
}

/**
 * Input for creating a new molecular property
 */
export interface CreateMolecularPropertyInput {
  moleculeId: Id<"molecules">;
  propertyTypeId: Id<"propertyTypes">;
  value: string | number | boolean | null;
  numericValue?: number;
  units?: string;
  source?: Id<"dataSources">;
  calculationMethod?: string;
  confidence?: number;
}

/**
 * Input for updating a molecular property
 */
export interface UpdateMolecularPropertyInput {
  value?: string | number | boolean | null;
  numericValue?: number;
  units?: string;
  source?: Id<"dataSources">;
  calculationMethod?: string;
  confidence?: number;
}

/**
 * Filter conditions for querying property types
 */
export interface PropertyTypeFilter {
  name?: string;
  dataType?: "string" | "number" | "boolean";
  category?: string;
  isCalculated?: boolean;
}

/**
 * Filter conditions for querying molecular properties
 */
export interface MolecularPropertyFilter {
  moleculeId?: Id<"molecules">;
  propertyTypeId?: Id<"propertyTypes">;
  propertyTypeName?: string;
  valueEquals?: string | number | boolean;
  valueRange?: {
    min?: number;
    max?: number;
  };
}

/**
 * Options for property type queries
 */
export interface PropertyTypeQueryOptions {
  limit?: number;
  cursor?: string;
}

/**
 * Options for molecular property queries
 */
export interface MolecularPropertyQueryOptions {
  limit?: number;
  cursor?: string;
}
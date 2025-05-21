/**
 * Type definitions for the molecule data model
 */

import { Id } from "../_generated/dataModel";

/**
 * Molecule document in the database
 */
export interface Molecule {
  _id: Id<"molecules">;
  _creationTime: number;
  
  // Core fields
  name: string;
  pubchemCid?: string;
  canonicalSmiles?: string;
  inchiKey?: string;
  formula?: string;
  
  // Status fields
  status: "active" | "deprecated" | "consolidated";
  consolidated?: boolean;
  consolidatedWith?: Id<"molecules">;
  
  // Metadata
  createdAt: number;
  updatedAt: number;
  dataSource?: Id<"dataSources">;
  sourceId?: string;
}

/**
 * Input for creating a new molecule
 */
export interface CreateMoleculeInput {
  name: string;
  pubchemCid?: string;
  canonicalSmiles?: string;
  inchiKey?: string;
  formula?: string;
  dataSource?: Id<"dataSources">;
  sourceId?: string;
}

/**
 * Input for updating a molecule
 */
export interface UpdateMoleculeInput {
  name?: string;
  pubchemCid?: string;
  canonicalSmiles?: string;
  inchiKey?: string;
  formula?: string;
  status?: "active" | "deprecated" | "consolidated";
  consolidated?: boolean;
  consolidatedWith?: Id<"molecules">;
  dataSource?: Id<"dataSources">;
  sourceId?: string;
}

/**
 * Molecule with additional derived data
 */
export interface MoleculeWithProperties extends Molecule {
  properties?: Array<{
    propertyTypeId: Id<"propertyTypes">;
    propertyTypeName: string;
    value: any;
    units?: string;
  }>;
}

/**
 * Filter conditions for querying molecules
 */
export interface MoleculeFilter {
  name?: string;
  pubchemCid?: string;
  inchiKey?: string;
  formula?: string;
  status?: "active" | "deprecated" | "consolidated";
  consolidated?: boolean;
  includeDeprecated?: boolean;
  
  // Property filters
  hasProperty?: Id<"propertyTypes">;
  propertyEquals?: {
    propertyTypeId: Id<"propertyTypes">;
    value: any;
  };
  propertyRange?: {
    propertyTypeId: Id<"propertyTypes">;
    min?: number;
    max?: number;
  };
}

/**
 * Result of a molecule consolidation operation
 */
export interface MoleculeConsolidationResult {
  primaryId: Id<"molecules">;
  consolidatedIds: Id<"molecules">[];
  totalProperties: number;
  mergedProperties: number;
}

/**
 * Options for molecule queries
 */
export interface MoleculeQueryOptions {
  limit?: number;
  cursor?: string;
  includeProperties?: boolean;
  propertiesFilter?: Id<"propertyTypes">[];
}
/**
 * Type definitions for molecular cross references and synonyms
 */

import { Id } from "../_generated/dataModel";

/**
 * Molecule cross reference document in the database
 */
export interface MoleculeCrossReference {
  _id: Id<"moleculeCrossReferences">;
  _creationTime: number;
  
  // Core fields
  moleculeId: Id<"molecules">;
  databaseName: string;
  identifier: string;
  url?: string;
  
  // Metadata
  createdAt: number;
  updatedAt: number;
}

/**
 * Input for creating a new molecule cross reference
 */
export interface CreateCrossReferenceInput {
  moleculeId: Id<"molecules">;
  databaseName: string;
  identifier: string;
  url?: string;
}

/**
 * Input for updating a molecule cross reference
 */
export interface UpdateCrossReferenceInput {
  databaseName?: string;
  identifier?: string;
  url?: string;
}

/**
 * Molecule synonym document in the database
 */
export interface MoleculeSynonym {
  _id: Id<"moleculeSynonyms">;
  _creationTime: number;
  
  // Core fields
  moleculeId: Id<"molecules">;
  name: string;
  type?: string;  // common, iupac, trade, etc.
  source?: Id<"dataSources">;
  
  // Metadata
  createdAt: number;
  updatedAt: number;
}

/**
 * Input for creating a new molecule synonym
 */
export interface CreateSynonymInput {
  moleculeId: Id<"molecules">;
  name: string;
  type?: string;
  source?: Id<"dataSources">;
}

/**
 * Input for updating a molecule synonym
 */
export interface UpdateSynonymInput {
  name?: string;
  type?: string;
  source?: Id<"dataSources">;
}

/**
 * Cross reference with expanded molecule data
 */
export interface CrossReferenceWithMolecule extends MoleculeCrossReference {
  molecule: {
    _id: Id<"molecules">;
    name: string;
    formula?: string;
  };
}

/**
 * Synonym with expanded molecule data
 */
export interface SynonymWithMolecule extends MoleculeSynonym {
  molecule: {
    _id: Id<"molecules">;
    name: string;
    formula?: string;
  };
}

/**
 * Filter conditions for querying cross references
 */
export interface CrossReferenceFilter {
  moleculeId?: Id<"molecules">;
  databaseName?: string;
  identifier?: string;
}

/**
 * Filter conditions for querying synonyms
 */
export interface SynonymFilter {
  moleculeId?: Id<"molecules">;
  name?: string;
  type?: string;
  source?: Id<"dataSources">;
}

/**
 * Options for cross reference queries
 */
export interface CrossReferenceQueryOptions {
  limit?: number;
  cursor?: string;
  includeMolecule?: boolean;
}

/**
 * Options for synonym queries
 */
export interface SynonymQueryOptions {
  limit?: number;
  cursor?: string;
  includeMolecule?: boolean;
}
/**
 * Type definitions for mixtures and mixture components
 */

import { Id } from "../_generated/dataModel";

/**
 * Mixture document in the database
 */
export interface Mixture {
  _id: Id<"mixtures">;
  _creationTime: number;
  
  // Core fields
  name: string;
  description?: string;
  type?: string;
  
  // Ownership
  createdBy?: Id<"users">;
  projectId?: Id<"projects">;
  
  // Metadata
  createdAt: number;
  updatedAt: number;
  public: boolean;
}

/**
 * Input for creating a new mixture
 */
export interface CreateMixtureInput {
  name: string;
  description?: string;
  type?: string;
  projectId?: Id<"projects">;
  public?: boolean;
}

/**
 * Input for updating a mixture
 */
export interface UpdateMixtureInput {
  name?: string;
  description?: string;
  type?: string;
  projectId?: Id<"projects">;
  public?: boolean;
}

/**
 * Mixture component document in the database
 */
export interface MixtureComponent {
  _id: Id<"mixtureComponents">;
  _creationTime: number;
  
  // Core fields
  mixtureId: Id<"mixtures">;
  moleculeId: Id<"molecules">;
  concentration: number;
  units: string;
  role?: string;
  notes?: string;
  
  // Metadata
  createdAt: number;
  updatedAt: number;
}

/**
 * Input for creating a new mixture component
 */
export interface CreateMixtureComponentInput {
  mixtureId: Id<"mixtures">;
  moleculeId: Id<"molecules">;
  concentration: number;
  units: string;
  role?: string;
  notes?: string;
}

/**
 * Input for updating a mixture component
 */
export interface UpdateMixtureComponentInput {
  concentration?: number;
  units?: string;
  role?: string;
  notes?: string;
}

/**
 * Mixture with expanded components data
 */
export interface MixtureWithComponents extends Mixture {
  components: Array<MixtureComponentWithDetails>;
}

/**
 * Mixture component with expanded molecule data
 */
export interface MixtureComponentWithDetails extends MixtureComponent {
  molecule?: {
    _id: Id<"molecules">;
    name: string;
    formula?: string;
    pubchemCid?: string;
  };
}

/**
 * Filter conditions for querying mixtures
 */
export interface MixtureFilter {
  name?: string;
  type?: string;
  createdBy?: Id<"users">;
  projectId?: Id<"projects">;
  public?: boolean;
  includeComponents?: boolean;
  
  // Component filters
  containsMolecule?: Id<"molecules">;
  moleculeConcentrationRange?: {
    moleculeId: Id<"molecules">;
    min?: number;
    max?: number;
  };
}

/**
 * Filter conditions for querying mixture components
 */
export interface MixtureComponentFilter {
  mixtureId?: Id<"mixtures">;
  moleculeId?: Id<"molecules">;
  role?: string;
  concentrationRange?: {
    min?: number;
    max?: number;
  };
}

/**
 * Options for mixture queries
 */
export interface MixtureQueryOptions {
  limit?: number;
  cursor?: string;
  includeComponents?: boolean;
}

/**
 * Options for mixture component queries
 */
export interface MixtureComponentQueryOptions {
  limit?: number;
  cursor?: string;
  includeMoleculeDetails?: boolean;
}
/**
 * Convex schema definitions for the frontend
 * This file defines the TypeScript types and interfaces for Convex data models
 */

import { Infer, v } from 'convex/values';
import { Doc, Id } from './_generated/dataModel';

// Helper types for Convex ID fields
export type IdOf<T extends string> = Id<T>;
export type WithId<T> = T & { _id: Id<any> };

/**
 * User schema - corresponds to the user table in Convex
 */
export const userSchema = {
  email: v.string(),
  name: v.optional(v.string()),
  role: v.optional(v.string()),
  roles: v.optional(v.array(v.string())),
  permissions: v.optional(v.array(v.string())),
  lastLoginAt: v.optional(v.number()),
  createdAt: v.number(),
  updatedAt: v.number(),
};

export type UserSchema = Infer<typeof userSchema>;
export type User = Doc<"users">;

/**
 * Molecule schema - corresponds to the molecules table in Convex
 */
export const moleculeSchema = {
  name: v.string(),
  inchikey: v.string(),
  smiles: v.string(),
  molecularFormula: v.string(),
  molecularWeight: v.number(),
  isConsolidated: v.boolean(),
  moleculeStatus: v.string(), // 'original' | 'primary' | 'duplicate'
  primaryMoleculeId: v.optional(v.id("molecules")),
  primaryMoleculeName: v.optional(v.string()),
  duplicateCount: v.optional(v.number()),
  createdAt: v.number(),
  updatedAt: v.number(),
};

export type MoleculeSchema = Infer<typeof moleculeSchema>;
export type Molecule = Doc<"molecules">;

/**
 * Molecule Property schema - corresponds to the moleculeProperties table in Convex
 */
export const moleculePropertySchema = {
  moleculeId: v.id("molecules"),
  propertyTypeId: v.string(),
  propertyName: v.string(),
  propertyType: v.string(),
  numericValue: v.optional(v.number()),
  textValue: v.optional(v.string()),
  unit: v.optional(v.string()),
  source: v.string(),
  createdAt: v.number(),
  updatedAt: v.number(),
};

export type MoleculePropertySchema = Infer<typeof moleculePropertySchema>;
export type MoleculeProperty = Doc<"moleculeProperties">;

/**
 * Mixture schema - corresponds to the mixtures table in Convex
 */
export const mixtureSchema = {
  name: v.string(),
  description: v.string(),
  componentCount: v.optional(v.number()),
  createdBy: v.optional(v.id("users")),
  createdAt: v.number(),
  updatedAt: v.number(),
};

export type MixtureSchema = Infer<typeof mixtureSchema>;
export type Mixture = Doc<"mixtures">;

/**
 * Mixture Component schema - corresponds to the mixtureComponents table in Convex
 */
export const mixtureComponentSchema = {
  mixtureId: v.id("mixtures"),
  moleculeId: v.id("molecules"),
  moleculeName: v.string(),
  concentration: v.number(),
  concentrationUnit: v.string(),
  role: v.optional(v.string()),
  createdAt: v.number(),
  updatedAt: v.number(),
};

export type MixtureComponentSchema = Infer<typeof mixtureComponentSchema>;
export type MixtureComponent = Doc<"mixtureComponents">;

/**
 * Experiment schema - corresponds to the experiments table in Convex
 */
export const experimentSchema = {
  name: v.string(),
  description: v.optional(v.string()),
  protocolId: v.id("protocols"),
  protocolVersion: v.optional(v.string()),
  tissueTypeId: v.id("tissueTypes"),
  experimentType: v.string(),
  startDate: v.string(),
  endDate: v.optional(v.string()),
  status: v.string(), // 'planned' | 'in_progress' | 'completed' | 'aborted' | 'failed'
  researcher: v.string(),
  labId: v.optional(v.string()),
  equipment: v.optional(v.array(v.string())),
  environmentalConditions: v.optional(v.map(v.string(), v.any())),
  notes: v.optional(v.string()),
  tags: v.optional(v.array(v.string())),
  createdBy: v.optional(v.id("users")),
  createdAt: v.number(),
  updatedAt: v.number(),
};

export type ExperimentSchema = Infer<typeof experimentSchema>;
export type Experiment = Doc<"experiments">;

/**
 * Experiment Result schema - corresponds to the experimentResults table in Convex
 */
export const experimentResultSchema = {
  experimentId: v.id("experiments"),
  tissueTypeId: v.id("tissueTypes"),
  moleculeId: v.optional(v.id("molecules")),
  mixtureId: v.optional(v.id("mixtures")),
  concentration: v.optional(v.number()),
  concentrationUnit: v.optional(v.string()),
  viabilityPercentage: v.optional(v.number()),
  recoveryRate: v.optional(v.number()),
  functionalityScore: v.optional(v.number()),
  resultDetails: v.optional(v.map(v.string(), v.any())),
  notes: v.optional(v.string()),
  protocolStepId: v.optional(v.id("protocolSteps")),
  timestamp: v.string(),
  createdBy: v.optional(v.id("users")),
  createdAt: v.number(),
  updatedAt: v.number(),
};

export type ExperimentResultSchema = Infer<typeof experimentResultSchema>;
export type ExperimentResult = Doc<"experimentResults">;

/**
 * Protocol schema - corresponds to the protocols table in Convex
 */
export const protocolSchema = {
  name: v.string(),
  description: v.optional(v.string()),
  version: v.string(),
  parentVersionId: v.optional(v.id("protocols")),
  parameters: v.optional(v.map(v.string(), v.any())),
  tags: v.optional(v.array(v.string())),
  isTemplate: v.boolean(),
  createdBy: v.optional(v.id("users")),
  createdAt: v.number(),
  updatedAt: v.number(),
};

export type ProtocolSchema = Infer<typeof protocolSchema>;
export type Protocol = Doc<"protocols">;

/**
 * Protocol Step schema - corresponds to the protocolSteps table in Convex
 */
export const protocolStepSchema = {
  protocolId: v.id("protocols"),
  name: v.string(),
  description: v.optional(v.string()),
  order: v.number(),
  duration: v.optional(v.number()),
  durationUnit: v.optional(v.string()),
  temperature: v.optional(v.number()),
  temperatureUnit: v.optional(v.string()),
  parameters: v.optional(v.map(v.string(), v.any())),
  equipment: v.optional(v.array(v.string())),
  createdBy: v.optional(v.id("users")),
  createdAt: v.number(),
  updatedAt: v.number(),
};

export type ProtocolStepSchema = Infer<typeof protocolStepSchema>;
export type ProtocolStep = Doc<"protocolSteps">;

/**
 * Tissue Type schema - corresponds to the tissueTypes table in Convex
 */
export const tissueTypeSchema = {
  name: v.string(),
  description: v.optional(v.string()),
  species: v.optional(v.string()),
  category: v.optional(v.string()),
  properties: v.optional(v.map(v.string(), v.any())),
  createdBy: v.optional(v.id("users")),
  createdAt: v.number(),
  updatedAt: v.number(),
};

export type TissueTypeSchema = Infer<typeof tissueTypeSchema>;
export type TissueType = Doc<"tissueTypes">;

/**
 * Time Series schema - corresponds to the timeSeries table in Convex
 */
export const timeSeriesSchema = {
  experimentId: v.id("experiments"),
  resultId: v.optional(v.id("experimentResults")),
  parameter: v.string(),
  unit: v.string(),
  startTime: v.string(),
  endTime: v.string(),
  notes: v.optional(v.string()),
  createdBy: v.optional(v.id("users")),
  createdAt: v.number(),
  updatedAt: v.number(),
};

export type TimeSeriesSchema = Infer<typeof timeSeriesSchema>;
export type TimeSeries = Doc<"timeSeries">;

/**
 * Time Series Data Point schema - corresponds to the timeSeriesDataPoints table in Convex
 */
export const timeSeriesDataPointSchema = {
  timeSeriesId: v.id("timeSeries"),
  time: v.number(),
  value: v.number(),
  uncertainty: v.optional(v.number()),
  metadata: v.optional(v.map(v.string(), v.any())),
  createdBy: v.optional(v.id("users")),
  createdAt: v.number(),
  updatedAt: v.number(),
};

export type TimeSeriesDataPointSchema = Infer<typeof timeSeriesDataPointSchema>;
export type TimeSeriesDataPoint = Doc<"timeSeriesDataPoints">;

// Export all schemas for use in other parts of the application
export const schemas = {
  user: userSchema,
  molecule: moleculeSchema,
  moleculeProperty: moleculePropertySchema,
  mixture: mixtureSchema,
  mixtureComponent: mixtureComponentSchema,
  experiment: experimentSchema,
  experimentResult: experimentResultSchema,
  protocol: protocolSchema,
  protocolStep: protocolStepSchema,
  tissueType: tissueTypeSchema,
  timeSeries: timeSeriesSchema,
  timeSeriesDataPoint: timeSeriesDataPointSchema,
};
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
 * Enhanced Experiment schema - corresponds to the enhancedExperiments table in Convex
 */
export const enhancedExperimentSchema = {
  name: v.string(),
  description: v.optional(v.string()),
  protocolId: v.optional(v.id("protocols")),
  protocolVersion: v.optional(v.string()),
  status: v.string(), // 'planned' | 'in_progress' | 'completed' | 'aborted' | 'failed'
  date: v.optional(v.number()),
  conductedBy: v.optional(v.id("users")),
  temperature: v.optional(v.number()),
  temperatureUnit: v.optional(v.string()),
  pressure: v.optional(v.number()),
  pressureUnit: v.optional(v.string()),
  coolingRate: v.optional(v.number()),
  coolingRateUnit: v.optional(v.string()),
  concentration: v.optional(v.number()),
  concentrationUnit: v.optional(v.string()),
  ph: v.optional(v.number()),
  mixingRate: v.optional(v.number()),
  mixingRateUnit: v.optional(v.string()),
  parameters: v.optional(v.map(v.string(), v.any())),
  equipment: v.optional(v.array(v.string())),
  notes: v.optional(v.string()),
  tags: v.optional(v.array(v.string())),
  isTemplate: v.optional(v.boolean()),
  createdAt: v.number(),
  updatedAt: v.number(),
};

export type EnhancedExperimentSchema = Infer<typeof enhancedExperimentSchema>;
export type EnhancedExperiment = Doc<"enhancedExperiments">;

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
 * Protocol Template schema - corresponds to the protocolTemplates table in Convex
 */
export const protocolTemplateSchema = {
  name: v.string(),
  description: v.optional(v.string()),
  version: v.string(),
  parentVersionId: v.optional(v.id("protocolTemplates")),
  parameters: v.optional(v.map(v.string(), v.any())),
  steps: v.array(v.object({
    name: v.string(),
    description: v.optional(v.string()),
    order: v.number(),
    duration: v.optional(v.number()),
    durationUnit: v.optional(v.string()),
    temperature: v.optional(v.number()),
    temperatureUnit: v.optional(v.string()),
    parameters: v.optional(v.map(v.string(), v.any())),
    equipment: v.optional(v.array(v.string())),
  })),
  category: v.optional(v.string()),
  tags: v.optional(v.array(v.string())),
  createdBy: v.optional(v.id("users")),
  createdAt: v.number(),
  updatedAt: v.number(),
};

export type ProtocolTemplateSchema = Infer<typeof protocolTemplateSchema>;
export type ProtocolTemplate = Doc<"protocolTemplates">;

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
  name: v.string(),
  description: v.optional(v.string()),
  experimentId: v.optional(v.id("enhancedExperiments")),
  units: v.optional(v.string()),
  startTimestamp: v.optional(v.number()),
  endTimestamp: v.optional(v.number()),
  parameter: v.string(),
  type: v.optional(v.string()),
  tags: v.optional(v.array(v.string())),
  metadata: v.optional(v.map(v.string(), v.any())),
  createdBy: v.optional(v.id("users")),
  createdAt: v.number(),
  updatedAt: v.number(),
};

export type TimeSeriesSchema = Infer<typeof timeSeriesSchema>;
export type TimeSeries = Doc<"timeSeries">;

/**
 * Time Series Data schema - corresponds to the timeSeriesData table in Convex
 */
export const timeSeriesDataSchema = {
  timeSeriesId: v.id("timeSeries"),
  timestamp: v.number(),
  value: v.number(),
  uncertainty: v.optional(v.number()),
  metadata: v.optional(v.map(v.string(), v.any())),
  createdAt: v.number(),
  updatedAt: v.number(),
};

export type TimeSeriesDataSchema = Infer<typeof timeSeriesDataSchema>;
export type TimeSeriesData = Doc<"timeSeriesData">;

/**
 * Uncertainty Analysis schema - corresponds to the uncertaintyAnalysis table in Convex
 */
export const uncertaintyAnalysisSchema = {
  timeSeriesId: v.id("timeSeries"),
  analysisType: v.string(), // 'confidence_interval', 'monte_carlo', 'sensitivity', 'error_propagation'
  results: v.any(),
  parameters: v.optional(v.object({
    confidenceLevel: v.optional(v.number()),
    monteCarloSamples: v.optional(v.number()),
    sensitivityParameters: v.optional(v.array(v.string())),
  })),
  createdAt: v.number(),
  updatedAt: v.number(),
};

export type UncertaintyAnalysisSchema = Infer<typeof uncertaintyAnalysisSchema>;
export type UncertaintyAnalysis = Doc<"uncertaintyAnalysis">;

/**
 * Lab Verification schema - corresponds to the labVerifications table in Convex
 */
export const labVerificationSchema = {
  experimentId: v.id("enhancedExperiments"),
  verificationStatus: v.string(), // 'pending', 'verified', 'rejected', 'needs_revision'
  requestDate: v.number(),
  verificationDate: v.optional(v.number()),
  verifiedBy: v.optional(v.id("users")),
  requestedBy: v.id("users"),
  equipmentUsed: v.union(v.string(), v.array(v.string())),
  methodologyDescription: v.optional(v.string()),
  controlProcedures: v.optional(v.string()),
  requestNotes: v.optional(v.string()),
  verifierNotes: v.optional(v.string()),
  reproducibilityRating: v.optional(v.number()),
  qualityRating: v.optional(v.number()),
  documentationRating: v.optional(v.number()),
  overallRating: v.optional(v.number()),
  evidenceUrls: v.optional(v.array(v.string())),
  createdAt: v.number(),
  updatedAt: v.number(),
};

export type LabVerificationSchema = Infer<typeof labVerificationSchema>;
export type LabVerification = Doc<"labVerifications">;

// Export all schemas for use in other parts of the application
export const schemas = {
  user: userSchema,
  molecule: moleculeSchema,
  moleculeProperty: moleculePropertySchema,
  mixture: mixtureSchema,
  mixtureComponent: mixtureComponentSchema,
  experiment: experimentSchema,
  enhancedExperiment: enhancedExperimentSchema,
  experimentResult: experimentResultSchema,
  protocol: protocolSchema,
  protocolStep: protocolStepSchema,
  protocolTemplate: protocolTemplateSchema,
  tissueType: tissueTypeSchema,
  timeSeries: timeSeriesSchema,
  timeSeriesData: timeSeriesDataSchema,
  uncertaintyAnalysis: uncertaintyAnalysisSchema,
  labVerification: labVerificationSchema,
};
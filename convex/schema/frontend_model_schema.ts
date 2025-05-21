/**
 * Frontend Model Schema for Convex
 * 
 * This schema implements the data models required by the frontend application,
 * ensuring alignment between the frontend types and the Convex database.
 */

import { defineSchema, defineTable } from "convex/server";
import { v } from "convex/values";

export default defineSchema({
  // User-related tables
  users: defineTable({
    email: v.string(),
    name: v.optional(v.string()),
    role: v.optional(v.string()),
    roles: v.optional(v.array(v.string())),
    permissions: v.optional(v.array(v.string())),
    lastLoginAt: v.optional(v.number()),
    createdAt: v.number(),
    updatedAt: v.number()
  })
    .index("by_email", ["email"])
    .index("by_role", ["role"]),

  // Molecule-related tables
  molecules: defineTable({
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
    updatedAt: v.number()
  })
    .index("by_name", ["name"])
    .index("by_inchikey", ["inchikey"])
    .index("by_status", ["moleculeStatus"])
    .index("by_primary", ["primaryMoleculeId"])
    .index("by_consolidated", ["isConsolidated"]),

  moleculeProperties: defineTable({
    moleculeId: v.id("molecules"),
    propertyTypeId: v.string(),
    propertyName: v.string(),
    propertyType: v.string(),
    numericValue: v.optional(v.number()),
    textValue: v.optional(v.string()),
    unit: v.optional(v.string()),
    source: v.string(),
    createdAt: v.number(),
    updatedAt: v.number()
  })
    .index("by_molecule", ["moleculeId"])
    .index("by_property_name", ["propertyName"])
    .index("by_property_type", ["propertyType"]),

  // Mixture-related tables
  mixtures: defineTable({
    name: v.string(),
    description: v.string(),
    componentCount: v.optional(v.number()),
    createdBy: v.optional(v.id("users")),
    createdAt: v.number(),
    updatedAt: v.number()
  })
    .index("by_name", ["name"])
    .index("by_creator", ["createdBy"]),

  mixtureComponents: defineTable({
    mixtureId: v.id("mixtures"),
    moleculeId: v.id("molecules"),
    moleculeName: v.string(),
    concentration: v.number(),
    concentrationUnit: v.string(),
    role: v.optional(v.string()),
    createdAt: v.number(),
    updatedAt: v.number()
  })
    .index("by_mixture", ["mixtureId"])
    .index("by_molecule", ["moleculeId"]),

  // Tissue types table
  tissueTypes: defineTable({
    name: v.string(),
    description: v.optional(v.string()),
    species: v.optional(v.string()),
    category: v.optional(v.string()),
    properties: v.optional(v.map(v.string(), v.any())),
    createdBy: v.optional(v.id("users")),
    createdAt: v.number(),
    updatedAt: v.number()
  })
    .index("by_name", ["name"])
    .index("by_species", ["species"])
    .index("by_category", ["category"]),

  // Protocol-related tables
  protocols: defineTable({
    name: v.string(),
    description: v.optional(v.string()),
    version: v.string(),
    parentVersionId: v.optional(v.id("protocols")),
    parameters: v.optional(v.map(v.string(), v.any())),
    tags: v.optional(v.array(v.string())),
    isTemplate: v.boolean(),
    createdBy: v.optional(v.id("users")),
    createdAt: v.number(),
    updatedAt: v.number()
  })
    .index("by_name", ["name"])
    .index("by_version", ["name", "version"])
    .index("by_template", ["isTemplate"])
    .index("by_tags", ["tags"]),

  protocolSteps: defineTable({
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
    updatedAt: v.number()
  })
    .index("by_protocol", ["protocolId"])
    .index("by_order", ["protocolId", "order"]),

  // Experiment-related tables
  experiments: defineTable({
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
    updatedAt: v.number()
  })
    .index("by_name", ["name"])
    .index("by_protocol", ["protocolId"])
    .index("by_tissue", ["tissueTypeId"])
    .index("by_type", ["experimentType"])
    .index("by_status", ["status"])
    .index("by_date", ["startDate"])
    .index("by_researcher", ["researcher"])
    .index("by_creator", ["createdBy"])
    .index("by_tags", ["tags"]),

  experimentResults: defineTable({
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
    updatedAt: v.number()
  })
    .index("by_experiment", ["experimentId"])
    .index("by_tissue", ["tissueTypeId"])
    .index("by_molecule", ["moleculeId"])
    .index("by_mixture", ["mixtureId"])
    .index("by_timestamp", ["timestamp"]),

  // Time series data tables
  timeSeries: defineTable({
    experimentId: v.id("experiments"),
    resultId: v.optional(v.id("experimentResults")),
    parameter: v.string(),
    unit: v.string(),
    startTime: v.string(),
    endTime: v.string(),
    notes: v.optional(v.string()),
    createdBy: v.optional(v.id("users")),
    createdAt: v.number(),
    updatedAt: v.number()
  })
    .index("by_experiment", ["experimentId"])
    .index("by_result", ["resultId"])
    .index("by_parameter", ["parameter"]),

  timeSeriesDataPoints: defineTable({
    timeSeriesId: v.id("timeSeries"),
    time: v.number(),
    value: v.number(),
    uncertainty: v.optional(v.number()),
    metadata: v.optional(v.map(v.string(), v.any())),
    createdBy: v.optional(v.id("users")),
    createdAt: v.number(),
    updatedAt: v.number()
  })
    .index("by_series", ["timeSeriesId"])
    .index("by_time", ["timeSeriesId", "time"]),
});
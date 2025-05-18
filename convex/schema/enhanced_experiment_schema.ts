/**
 * Enhanced Experimental Data Schema for Convex
 * 
 * This schema implements the enhanced experimental data model with support for:
 * - Protocol templates and versioning
 * - Time-series data
 * - Equipment tracking
 * - Validation rules
 * - Uncertainty quantification
 */

import { defineSchema, defineTable } from "convex/server";
import { v } from "convex/values";

/**
 * Enhanced experimental schema for scientific data management
 */
export default defineSchema({
  // EXPERIMENT TYPE TABLES

  /**
   * Protocol templates - reusable experiment protocols with versioning
   */
  protocols: defineTable({
    name: v.string(),
    description: v.optional(v.string()),
    steps: v.array(v.object({
      name: v.string(),
      description: v.optional(v.string()),
      parameters: v.optional(v.map(v.string(), v.any())),
      duration: v.optional(v.number()),
      durationUnit: v.optional(v.string()),
      temperature: v.optional(v.number()),
      temperatureUnit: v.optional(v.string())
    })),
    parameters: v.optional(v.map(v.string(), v.any())),
    version: v.number(),
    parentId: v.optional(v.id("protocols")),
    isTemplate: v.boolean(),
    category: v.optional(v.string()),
    createdBy: v.optional(v.id("users")),
    createdAt: v.number(),
    updatedAt: v.number(),
    public: v.boolean()
  })
    .index("by_name", ["name"])
    .index("by_category", ["category"])
    .index("by_version", ["name", "version"])
    .index("by_creator", ["createdBy"])
    .index("by_template", ["isTemplate"])
    .index("by_public", ["public"]),

  /**
   * Tissue types - biological samples used in experiments
   */
  tissueTypes: defineTable({
    name: v.string(),
    description: v.optional(v.string()),
    species: v.optional(v.string()),
    taxonomyId: v.optional(v.number()),
    properties: v.optional(v.map(v.string(), v.any())),
    category: v.optional(v.string()),
    source: v.optional(v.string()),
    createdBy: v.optional(v.id("users")),
    createdAt: v.number(),
    updatedAt: v.number(),
    public: v.boolean()
  })
    .index("by_name", ["name"])
    .index("by_species", ["species"])
    .index("by_category", ["category"])
    .index("by_creator", ["createdBy"])
    .index("by_public", ["public"]),

  /**
   * Experiments - enhanced with protocol and equipment tracking
   */
  enhancedExperiments: defineTable({
    name: v.string(),
    description: v.optional(v.string()),
    experimentTypeId: v.optional(v.string()),
    protocolId: v.optional(v.id("protocols")),
    mixtureId: v.optional(v.id("mixtures")),
    temperature: v.optional(v.number()),
    temperatureUnit: v.optional(v.string()),
    coolingRate: v.optional(v.number()),
    coolingRateUnit: v.optional(v.string()),
    thawingRate: v.optional(v.number()),
    thawingRateUnit: v.optional(v.string()),
    pressure: v.optional(v.number()),
    pressureUnit: v.optional(v.string()),
    parameters: v.optional(v.map(v.string(), v.any())),
    version: v.number(),
    provenance: v.optional(v.map(v.string(), v.any())),
    conductedBy: v.optional(v.id("users")),
    projectId: v.optional(v.id("projects")),
    date: v.optional(v.number()),
    status: v.string(), // planned, in-progress, completed, failed
    tags: v.optional(v.array(v.string())),
    createdBy: v.optional(v.id("users")),
    createdAt: v.number(),
    updatedAt: v.number(),
    public: v.boolean()
  })
    .index("by_name", ["name"])
    .index("by_protocol", ["protocolId"])
    .index("by_mixture", ["mixtureId"])
    .index("by_user", ["conductedBy"])
    .index("by_project", ["projectId"])
    .index("by_status", ["status"])
    .index("by_date", ["date"])
    .index("by_public", ["public"])
    .index("by_tags", ["tags"]),

  /**
   * Enhanced experiment results with uncertainty tracking
   */
  enhancedExperimentResults: defineTable({
    experimentId: v.id("enhancedExperiments"),
    moleculeId: v.optional(v.id("molecules")),
    mixtureId: v.optional(v.id("mixtures")),
    tissueTypeId: v.id("tissueTypes"),
    parameterName: v.string(),
    value: v.union(v.string(), v.number(), v.boolean(), v.null()),
    numericValue: v.optional(v.number()),
    units: v.optional(v.string()),
    uncertainty: v.optional(v.object({
      type: v.string(), // standard_deviation, range, confidence_interval
      value: v.union(v.number(), v.array(v.number())),
      confidence: v.optional(v.number())
    })),
    provenance: v.optional(v.object({
      method: v.string(),
      reference: v.optional(v.string()),
      timestamp: v.number(),
      operator: v.optional(v.string())
    })),
    notes: v.optional(v.string()),
    createdBy: v.optional(v.id("users")),
    createdAt: v.number(),
    updatedAt: v.number()
  })
    .index("by_experiment", ["experimentId"])
    .index("by_molecule", ["moleculeId"])
    .index("by_mixture", ["mixtureId"])
    .index("by_tissue", ["tissueTypeId"])
    .index("by_parameter", ["parameterName"]),

  /**
   * Equipment types - categories of laboratory equipment
   */
  equipmentTypes: defineTable({
    name: v.string(),
    description: v.optional(v.string()),
    manufacturer: v.optional(v.string()),
    category: v.optional(v.string()),
    createdBy: v.optional(v.id("users")),
    createdAt: v.number(),
    updatedAt: v.number()
  })
    .index("by_name", ["name"])
    .index("by_category", ["category"])
    .index("by_manufacturer", ["manufacturer"]),

  /**
   * Equipment - specific equipment instances
   */
  equipment: defineTable({
    equipmentTypeId: v.id("equipmentTypes"),
    name: v.string(),
    model: v.optional(v.string()),
    serialNumber: v.optional(v.string()),
    description: v.optional(v.string()),
    calibrationDate: v.optional(v.number()),
    nextCalibrationDate: v.optional(v.number()),
    location: v.optional(v.string()),
    parameters: v.optional(v.map(v.string(), v.any())),
    createdBy: v.optional(v.id("users")),
    createdAt: v.number(),
    updatedAt: v.number()
  })
    .index("by_type", ["equipmentTypeId"])
    .index("by_name", ["name"])
    .index("by_location", ["location"]),

  /**
   * Experiment equipment - links experiments to equipment used
   */
  experimentEquipment: defineTable({
    experimentId: v.id("enhancedExperiments"),
    equipmentId: v.id("equipment"),
    role: v.optional(v.string()),
    parameters: v.optional(v.map(v.string(), v.any())),
    notes: v.optional(v.string()),
    createdBy: v.optional(v.id("users")),
    createdAt: v.number(),
    updatedAt: v.number()
  })
    .index("by_experiment", ["experimentId"])
    .index("by_equipment", ["equipmentId"]),

  /**
   * Time series - time-series data for experiments
   */
  timeSeries: defineTable({
    experimentId: v.id("enhancedExperiments"),
    name: v.string(),
    description: v.optional(v.string()),
    parameterName: v.string(),
    units: v.optional(v.string()),
    metadata: v.optional(v.map(v.string(), v.any())),
    createdBy: v.optional(v.id("users")),
    createdAt: v.number(),
    updatedAt: v.number()
  })
    .index("by_experiment", ["experimentId"])
    .index("by_parameter", ["parameterName"]),

  /**
   * Time series data - individual data points for time series
   */
  timeSeriesData: defineTable({
    timeSeriesId: v.id("timeSeries"),
    timestamp: v.number(),
    value: v.number(),
    uncertainty: v.optional(v.number()),
    metadata: v.optional(v.map(v.string(), v.any())),
    createdBy: v.optional(v.id("users")),
    createdAt: v.number(),
    updatedAt: v.number()
  })
    .index("by_series", ["timeSeriesId"])
    .index("by_timestamp", ["timeSeriesId", "timestamp"]),

  /**
   * Validation rules - rules for validating experimental data
   */
  validationRules: defineTable({
    parameterName: v.string(),
    ruleType: v.string(), // range, pattern, comparison, custom
    parameters: v.map(v.string(), v.any()),
    description: v.optional(v.string()),
    severity: v.string(), // error, warning, info
    createdBy: v.optional(v.id("users")),
    createdAt: v.number(),
    updatedAt: v.number()
  })
    .index("by_parameter", ["parameterName"])
    .index("by_type", ["ruleType"])
    .index("by_severity", ["severity"]),
});
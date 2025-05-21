/**
 * CryoProtect Convex Schema Definition
 * 
 * This file defines the document schema for the CryoProtect database in Convex.
 * It includes all tables, fields, types, and indexing patterns.
 */

import { defineSchema, defineTable } from "convex/server";
import { v } from "convex/values";

/**
 * Schema definition for the CryoProtect Convex database.
 * 
 * The schema is designed to optimize for document-oriented storage while
 * maintaining the integrity and relationships of the scientific data.
 */
export default defineSchema({
  // CORE TABLES

  /**
   * Molecules table - stores basic molecule information.
   * Each molecule has unique identifiers, names, and structure information.
   */
  molecules: defineTable({
    // Identifiers
    name: v.string(),
    pubchemCid: v.optional(v.string()),
    canonicalSmiles: v.optional(v.string()),
    inchiKey: v.optional(v.string()),
    formula: v.optional(v.string()),
    
    // Status fields
    status: v.string(), // active, deprecated, consolidated
    consolidated: v.optional(v.boolean()),
    consolidatedWith: v.optional(v.id("molecules")),
    
    // Metadata
    createdAt: v.number(),
    updatedAt: v.number(),
    dataSource: v.optional(v.id("dataSources")),
    sourceId: v.optional(v.string())
  })
    .index("by_pubchemCid", ["pubchemCid"])
    .index("by_inchiKey", ["inchiKey"])
    .index("by_name", ["name"])
    .index("by_status", ["status"])
    .index("by_consolidated", ["consolidated"]),

  /**
   * Molecular properties table - stores properties of molecules
   * Properties can have different types, units, and sources
   */
  molecularProperties: defineTable({
    moleculeId: v.id("molecules"),
    propertyTypeId: v.id("propertyTypes"),
    value: v.union(v.string(), v.number(), v.boolean(), v.null()),
    numericValue: v.optional(v.number()),
    units: v.optional(v.string()),
    source: v.optional(v.id("dataSources")),
    calculationMethod: v.optional(v.string()),
    confidence: v.optional(v.number()),
    createdAt: v.number(),
    updatedAt: v.number()
  })
    .index("by_molecule", ["moleculeId"])
    .index("by_property_type", ["propertyTypeId"])
    .index("by_molecule_property", ["moleculeId", "propertyTypeId"]),

  /**
   * Property types - defines the available property types and their characteristics
   */
  propertyTypes: defineTable({
    name: v.string(),
    displayName: v.string(),
    description: v.optional(v.string()),
    dataType: v.string(), // string, number, boolean
    units: v.optional(v.string()),
    defaultUnits: v.optional(v.string()),
    category: v.optional(v.string()),
    isCalculated: v.optional(v.boolean()),
    createdAt: v.number(),
    updatedAt: v.number()
  })
    .index("by_name", ["name"])
    .index("by_category", ["category"]),

  /**
   * Mixtures - defines combinations of molecules in specific concentrations
   */
  mixtures: defineTable({
    name: v.string(),
    description: v.optional(v.string()),
    type: v.optional(v.string()),
    createdBy: v.optional(v.id("users")),
    projectId: v.optional(v.id("projects")),
    createdAt: v.number(),
    updatedAt: v.number(),
    public: v.boolean()
  })
    .index("by_creator", ["createdBy"])
    .index("by_project", ["projectId"])
    .index("by_public", ["public"]),

  /**
   * Mixture Components - defines the molecules in a mixture and their proportions
   */
  mixtureComponents: defineTable({
    mixtureId: v.id("mixtures"),
    moleculeId: v.id("molecules"),
    concentration: v.number(),
    units: v.string(),
    role: v.optional(v.string()),
    notes: v.optional(v.string()),
    createdAt: v.number(),
    updatedAt: v.number()
  })
    .index("by_mixture", ["mixtureId"])
    .index("by_molecule", ["moleculeId"]),

  /**
   * Experiments - tracks laboratory experiments
   */
  experiments: defineTable({
    name: v.string(),
    description: v.optional(v.string()),
    mixtureId: v.optional(v.id("mixtures")),
    protocol: v.optional(v.string()),
    conductedBy: v.optional(v.id("users")),
    projectId: v.optional(v.id("projects")),
    date: v.optional(v.number()),
    status: v.string(), // planned, in-progress, completed, failed
    createdAt: v.number(),
    updatedAt: v.number(),
    public: v.boolean()
  })
    .index("by_mixture", ["mixtureId"])
    .index("by_user", ["conductedBy"])
    .index("by_project", ["projectId"])
    .index("by_status", ["status"])
    .index("by_public", ["public"]),

  /**
   * Enhanced Experiments - tracks experiments with improved data structure
   */
  enhancedExperiments: defineTable({
    title: v.string(),
    description: v.optional(v.string()),
    experimentTypeId: v.string(),
    protocolId: v.optional(v.id("protocols")),
    datePerformed: v.optional(v.number()),
    temperature: v.optional(v.number()),
    coolingRate: v.optional(v.number()),
    thawingRate: v.optional(v.number()),
    parameters: v.optional(v.map(v.string(), v.any())),
    version: v.number(),
    status: v.string(), // planned, in-progress, completed, failed
    projectId: v.optional(v.id("projects")),
    createdBy: v.optional(v.id("users")),
    createdAt: v.number(),
    updatedAt: v.number(),
    public: v.boolean()
  })
    .index("by_title", ["title"])
    .index("by_experiment_type", ["experimentTypeId"])
    .index("by_protocol", ["protocolId"])
    .index("by_project", ["projectId"])
    .index("by_creator", ["createdBy"])
    .index("by_status", ["status"])
    .index("by_public", ["public"]),

  /**
   * Protocols - defines experimental protocols
   */
  protocols: defineTable({
    name: v.string(),
    description: v.optional(v.string()),
    steps: v.array(v.object({
      id: v.string(),
      name: v.string(),
      description: v.optional(v.string()),
      parameters: v.optional(v.map(v.string(), v.any())),
      duration: v.optional(v.number()),
      durationUnit: v.optional(v.string()),
      temperature: v.optional(v.number()),
      temperatureUnit: v.optional(v.string())
    })),
    version: v.optional(v.number()),
    parentId: v.optional(v.id("protocols")),
    category: v.optional(v.string()),
    isTemplate: v.optional(v.boolean()),
    parameters: v.optional(v.map(v.string(), v.any())),
    createdBy: v.optional(v.id("users")),
    createdAt: v.number(),
    updatedAt: v.number(),
    public: v.optional(v.boolean())
  })
    .index("by_name", ["name"])
    .index("by_creator", ["createdBy"])
    .index("by_category", ["category"])
    .index("by_template", ["isTemplate"])
    .index("by_public", ["public"]),

  /**
   * Protocol-Experiment Links - connects protocols to experiments
   */
  protocol_experiment_links: defineTable({
    protocolId: v.id("protocols"),
    experimentId: v.string(),
    createdBy: v.optional(v.id("users")),
    createdAt: v.number(),
    metadata: v.optional(v.map(v.string(), v.any()))
  })
    .index("by_protocol", ["protocolId"])
    .index("by_experiment", ["experimentId"])
    .index("by_protocol_experiment", ["protocolId", "experimentId"]),

  /**
   * Protocol Step Results - stores results from executing protocol steps
   */
  protocol_step_results: defineTable({
    linkId: v.id("protocol_experiment_links"),
    stepId: v.string(),
    protocolId: v.id("protocols"),
    experimentId: v.string(),
    results: v.any(),
    createdBy: v.optional(v.id("users")),
    createdAt: v.number(),
    updatedAt: v.number(),
    updatedBy: v.optional(v.id("users"))
  })
    .index("by_link", ["linkId"])
    .index("by_link_step", ["linkId", "stepId"])
    .index("by_experiment", ["experimentId"])
    .index("by_protocol_step", ["protocolId", "stepId"]),

  /**
   * Experiment Type - defines different types of experiments
   */
  experimentTypes: defineTable({
    name: v.string(),
    description: v.optional(v.string()),
    category: v.optional(v.string()),
    defaultProtocolId: v.optional(v.id("protocols")),
    parameters: v.optional(v.map(v.string(), v.any())),
    createdBy: v.optional(v.id("users")),
    createdAt: v.number(),
    updatedAt: v.number()
  })
    .index("by_name", ["name"])
    .index("by_category", ["category"]),

  /**
   * Tissue Types - defines biological tissue types used in experiments
   */
  tissueTypes: defineTable({
    name: v.string(),
    description: v.optional(v.string()),
    species: v.optional(v.string()),
    taxonomyId: v.optional(v.number()),
    properties: v.optional(v.map(v.string(), v.any())),
    createdBy: v.optional(v.id("users")),
    createdAt: v.number(),
    updatedAt: v.number()
  })
    .index("by_name", ["name"])
    .index("by_species", ["species"]),

  /**
   * Experiment Results - stores results from experiments
   */
  experimentResults: defineTable({
    experimentId: v.id("experiments"),
    parameterName: v.string(),
    value: v.union(v.string(), v.number(), v.boolean(), v.null()),
    numericValue: v.optional(v.number()),
    units: v.optional(v.string()),
    notes: v.optional(v.string()),
    createdAt: v.number(),
    updatedAt: v.number()
  })
    .index("by_experiment", ["experimentId"]),

  /**
   * Enhanced Experiment Results - stores results with improved structure
   */
  enhancedExperimentResults: defineTable({
    experimentId: v.union(v.string(), v.id("enhancedExperiments")),
    tissueTypeId: v.union(v.string(), v.id("tissueTypes")),
    moleculeId: v.optional(v.id("molecules")),
    mixtureId: v.optional(v.id("mixtures")),
    concentration: v.optional(v.number()),
    concentrationUnit: v.optional(v.string()),
    viabilityPercentage: v.optional(v.number()),
    recoveryRate: v.optional(v.number()),
    functionalityScore: v.optional(v.number()),
    uncertainty: v.optional(v.map(v.string(), v.map(v.string(), v.any()))),
    resultDetails: v.optional(v.map(v.string(), v.any())),
    notes: v.optional(v.string()),
    createdBy: v.optional(v.id("users")),
    createdAt: v.number(),
    updatedAt: v.number()
  })
    .index("by_experiment", ["experimentId"])
    .index("by_tissue_type", ["tissueTypeId"])
    .index("by_molecule", ["moleculeId"])
    .index("by_mixture", ["mixtureId"]),

  /**
   * Predictions - stores model predictions for molecules or mixtures
   */
  predictions: defineTable({
    modelId: v.id("scientificModels"),
    moleculeId: v.optional(v.id("molecules")),
    mixtureId: v.optional(v.id("mixtures")),
    parameterName: v.string(),
    value: v.union(v.string(), v.number(), v.boolean(), v.null()),
    numericValue: v.optional(v.number()),
    units: v.optional(v.string()),
    confidence: v.optional(v.number()),
    createdAt: v.number(),
    updatedAt: v.number()
  })
    .index("by_molecule", ["moleculeId"])
    .index("by_mixture", ["mixtureId"])
    .index("by_model", ["modelId"]),

  /**
   * Scientific Models - defines the prediction models used
   */
  scientificModels: defineTable({
    name: v.string(),
    description: v.optional(v.string()),
    version: v.string(),
    type: v.string(),
    parameters: v.optional(v.string()),
    createdBy: v.optional(v.id("users")),
    createdAt: v.number(),
    updatedAt: v.number(),
    public: v.boolean()
  })
    .index("by_type", ["type"])
    .index("by_creator", ["createdBy"])
    .index("by_version", ["name", "version"]),

  // SUPPORTING TABLES

  /**
   * Data Sources - defines sources of molecule data
   */
  dataSources: defineTable({
    name: v.string(),
    description: v.optional(v.string()),
    url: v.optional(v.string()),
    type: v.string(), // database, publication, experiment, calculation
    version: v.optional(v.string()),
    createdAt: v.number(),
    updatedAt: v.number()
  })
    .index("by_name", ["name"])
    .index("by_type", ["type"]),

  /**
   * Molecule Synonyms - stores alternative names for molecules
   */
  moleculeSynonyms: defineTable({
    moleculeId: v.id("molecules"),
    name: v.string(),
    type: v.optional(v.string()), // common, iupac, trade, etc.
    source: v.optional(v.id("dataSources")),
    createdAt: v.number(),
    updatedAt: v.number()
  })
    .index("by_molecule", ["moleculeId"])
    .index("by_name", ["name"]),

  /**
   * Molecule Cross References - links to external databases
   */
  moleculeCrossReferences: defineTable({
    moleculeId: v.id("molecules"),
    databaseName: v.string(), // ChEMBL, DrugBank, etc.
    identifier: v.string(),
    url: v.optional(v.string()),
    createdAt: v.number(),
    updatedAt: v.number()
  })
    .index("by_molecule", ["moleculeId"])
    .index("by_database", ["databaseName"]),

  /**
   * Toxicity Data - stores toxicity information for molecules
   */
  toxicityData: defineTable({
    moleculeId: v.id("molecules"),
    assayId: v.id("toxicityAssays"),
    result: v.union(v.string(), v.number(), v.boolean()),
    numericValue: v.optional(v.number()),
    units: v.optional(v.string()),
    source: v.optional(v.id("dataSources")),
    createdAt: v.number(),
    updatedAt: v.number()
  })
    .index("by_molecule", ["moleculeId"])
    .index("by_assay", ["assayId"]),

  /**
   * Toxicity Assays - defines the toxicity tests available
   */
  toxicityAssays: defineTable({
    name: v.string(),
    description: v.optional(v.string()),
    type: v.string(),
    endpoint: v.string(),
    species: v.optional(v.string()),
    createdAt: v.number(),
    updatedAt: v.number()
  })
    .index("by_name", ["name"])
    .index("by_type", ["type"]),

  /**
   * Users - stores user information
   */
  users: defineTable({
    email: v.string(),
    name: v.optional(v.string()),
    role: v.string(), // admin, scientist, viewer
    lastLogin: v.optional(v.number()),
    createdAt: v.number(),
    updatedAt: v.number()
  })
    .index("by_email", ["email"])
    .index("by_role", ["role"]),

  /**
   * Projects - organizes work into projects
   */
  projects: defineTable({
    name: v.string(),
    description: v.optional(v.string()),
    ownerId: v.id("users"),
    status: v.string(), // active, archived, completed
    createdAt: v.number(),
    updatedAt: v.number(),
    public: v.boolean()
  })
    .index("by_owner", ["ownerId"])
    .index("by_status", ["status"])
    .index("by_public", ["public"]),

  /**
   * Team Members - associates users with projects
   */
  teamMembers: defineTable({
    projectId: v.id("projects"),
    userId: v.id("users"),
    role: v.string(), // owner, editor, viewer
    createdAt: v.number(),
    updatedAt: v.number()
  })
    .index("by_project", ["projectId"])
    .index("by_user", ["userId"])
    .index("by_project_user", ["projectId", "userId"]),

  /**
   * Scientific Data Audit - tracks changes to scientific data
   */
  scientificDataAudit: defineTable({
    table: v.string(),
    documentId: v.id(), // Generic ID reference
    operation: v.string(), // create, update, delete
    userId: v.optional(v.id("users")),
    previousValue: v.optional(v.any()),
    newValue: v.optional(v.any()),
    timestamp: v.number(),
    reason: v.optional(v.string())
  })
    .index("by_document", ["table", "documentId"])
    .index("by_user", ["userId"])
    .index("by_operation", ["operation"])
    .index("by_timestamp", ["timestamp"]),
});
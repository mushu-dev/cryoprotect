import { defineSchema, defineTable } from "convex/server";
import { v } from "convex/values";

export default defineSchema({
  molecules: defineTable({
    name: v.string(),
    pubchemCid: v.optional(v.string()),
    canonicalSmiles: v.optional(v.string()),
    inchiKey: v.optional(v.string()),
    formula: v.optional(v.string()),
    status: v.string(),
    createdAt: v.optional(v.string()),
    updatedAt: v.optional(v.string()),
  })
    .index("by_pubchem_cid", ["pubchemCid"])
    .index("by_status", ["status"])
    .index("by_name", ["name"]),

  // Molecular properties calculated by RDKit
  molecularProperties: defineTable({
    moleculeId: v.id("molecules"),
    
    // Basic descriptors
    molecularWeight: v.optional(v.number()),
    exactMass: v.optional(v.number()),
    logP: v.optional(v.number()),
    tpsa: v.optional(v.number()), // Topological Polar Surface Area
    
    // H-bond properties
    hbondDonors: v.optional(v.number()),
    hbondAcceptors: v.optional(v.number()),
    
    // Rotatable bonds and ring info
    rotatableBonds: v.optional(v.number()),
    aromaticRings: v.optional(v.number()),
    aliphaticRings: v.optional(v.number()),
    
    // Complexity metrics
    complexity: v.optional(v.number()),
    heavyAtomCount: v.optional(v.number()),
    
    // Calculated fingerprints (as strings for storage)
    morganFingerprint: v.optional(v.string()),
    rdkitFingerprint: v.optional(v.string()),
    
    // Calculation metadata
    calculatedAt: v.optional(v.string()),
    calculationVersion: v.optional(v.string()),
    status: v.optional(v.string()), // 'calculated', 'error', 'pending'
    
    // Legacy fields (for backward compatibility)
    calculationMethod: v.optional(v.string()),
    numericValue: v.optional(v.number()),
    propertyTypeId: v.optional(v.string()),
    source: v.optional(v.union(v.string(), v.null())),
    units: v.optional(v.string()),
    value: v.optional(v.number()),
  })
    .index("by_molecule", ["moleculeId"])
    .index("by_status", ["status"]),

  // Cryoprotectant scoring system
  cryoprotectantScores: defineTable({
    moleculeId: v.id("molecules"),
    
    // Core scoring metrics
    glassTempScore: v.optional(v.number()), // Glass transition temperature effectiveness
    viscosityScore: v.optional(v.number()), // Viscosity at freezing temps
    permeabilityScore: v.optional(v.number()), // Cell membrane permeability
    toxicityScore: v.optional(v.number()), // Cytotoxicity (lower is better)
    
    // Combined scores
    overallScore: v.number(), // Weighted combination of all metrics
    category: v.string(), // 'excellent', 'good', 'fair', 'poor'
    
    // Scoring metadata
    scoringAlgorithmVersion: v.optional(v.string()),
    calculatedAt: v.optional(v.string()),
    confidence: v.optional(v.number()), // 0-1 confidence in score
  })
    .index("by_molecule", ["moleculeId"])
    .index("by_overall_score", ["overallScore"])
    .index("by_category", ["category"]),

  // Experimental data linkage
  experimentalData: defineTable({
    moleculeId: v.id("molecules"),
    
    // Experimental conditions
    temperature: v.optional(v.number()),
    concentration: v.optional(v.number()),
    concentrationUnit: v.optional(v.string()),
    
    // Measured properties
    measuredViscosity: v.optional(v.number()),
    measuredGlassTemp: v.optional(v.number()),
    cellViability: v.optional(v.number()),
    
    // Study metadata
    studyReference: v.optional(v.string()),
    methodology: v.optional(v.string()),
    qualityScore: v.optional(v.number()),
    
    createdAt: v.string(),
    updatedAt: v.string(),
  })
    .index("by_molecule", ["moleculeId"])
    .index("by_temperature", ["temperature"]),

  // Calculation jobs tracking
  calculationJobs: defineTable({
    moleculeId: v.id("molecules"),
    jobType: v.string(), // 'properties', 'scoring', 'similarity'
    status: v.string(), // 'pending', 'running', 'completed', 'failed'
    
    // Job parameters
    parameters: v.optional(v.string()), // JSON string of parameters
    
    // Results
    result: v.optional(v.string()), // JSON string of results
    errorMessage: v.optional(v.string()),
    
    // Timing
    startedAt: v.optional(v.string()),
    completedAt: v.optional(v.string()),
    createdAt: v.string(),
  })
    .index("by_molecule", ["moleculeId"])
    .index("by_status", ["status"])
    .index("by_job_type", ["jobType"]),

  // Molecular similarity cache
  molecularSimilarity: defineTable({
    molecule1Id: v.id("molecules"),
    molecule2Id: v.id("molecules"),
    
    // Similarity metrics
    tanimotoSimilarity: v.optional(v.number()),
    diceCoefficient: v.optional(v.number()),
    cosineSimilarity: v.optional(v.number()),
    
    // Calculation metadata
    calculatedAt: v.string(),
    fingerprintType: v.string(), // 'morgan', 'rdkit', etc.
  })
    .index("by_molecule1", ["molecule1Id"])
    .index("by_molecule2", ["molecule2Id"])
    .index("by_similarity", ["tanimotoSimilarity"]),
});
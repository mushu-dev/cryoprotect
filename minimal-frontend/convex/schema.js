import { defineSchema, defineTable } from "convex/server";
import { v } from "convex/values";

/**
 * Schema definition for the CryoProtect Convex database - minimal version.
 * 
 * This is a simplified version of the schema for demonstration purposes.
 */
export default defineSchema({
  // CORE TABLES

  /**
   * Molecules table - stores basic molecule information.
   */
  molecules: defineTable({
    // Identifiers
    name: v.string(),
    pubchemCid: v.optional(v.string()),
    formula: v.optional(v.string()),
    smiles: v.optional(v.string()),
    
    // Properties
    molecularWeight: v.optional(v.number()),
    isCryoprotectant: v.optional(v.boolean()),
    description: v.optional(v.string()),
    
    // Metadata
    createdAt: v.number(),
    updatedAt: v.number(),
  })
    .index("by_name", ["name"])
    .index("by_pubchemCid", ["pubchemCid"]),

  /**
   * Molecular properties table - stores properties of molecules
   */
  molecularProperties: defineTable({
    moleculeId: v.id("molecules"),
    name: v.string(),
    value: v.union(v.string(), v.number(), v.boolean(), v.null()),
    units: v.optional(v.string()),
    createdAt: v.number(),
  })
    .index("by_molecule", ["moleculeId"]),

  /**
   * Mixtures - defines combinations of molecules in specific concentrations
   */
  mixtures: defineTable({
    name: v.string(),
    description: v.optional(v.string()),
    freezingPoint: v.optional(v.number()),
    createdAt: v.number(),
    updatedAt: v.number(),
  })
    .index("by_name", ["name"]),

  /**
   * Mixture Components - defines the molecules in a mixture and their proportions
   */
  mixtureComponents: defineTable({
    mixtureId: v.id("mixtures"),
    moleculeId: v.id("molecules"),
    concentration: v.number(),
    units: v.string(),
    role: v.optional(v.string()),
    createdAt: v.number(),
  })
    .index("by_mixture", ["mixtureId"])
    .index("by_molecule", ["moleculeId"]),
});
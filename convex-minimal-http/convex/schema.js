import { defineSchema, defineTable } from "convex/server";
import { v } from "convex/values";

/**
 * Minimal schema definition for the CryoProtect Convex database.
 */
export default defineSchema({
  /**
   * Molecules table - stores basic molecule information.
   */
  molecules: defineTable({
    name: v.string(),
    pubchemCid: v.optional(v.string()),
    formula: v.optional(v.string()),
    smiles: v.optional(v.string()),
    molecularWeight: v.optional(v.number()),
    isCryoprotectant: v.optional(v.boolean()),
    description: v.optional(v.string()),
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
    createdAt: v.number(),
  })
    .index("by_mixture", ["mixtureId"])
    .index("by_molecule", ["moleculeId"]),
});
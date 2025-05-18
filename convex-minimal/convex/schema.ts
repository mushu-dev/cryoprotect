import { defineSchema, defineTable } from "convex/server";
import { v } from "convex/values";

// Define a minimal schema for our API functions
export default defineSchema({
  // Define example tables
  molecules: defineTable({
    name: v.string(),
    formula: v.optional(v.string()),
    smiles: v.optional(v.string()),
    inchiKey: v.optional(v.string()),
    pubchemCid: v.optional(v.string()),
    status: v.string(),
    createdAt: v.number(),
    updatedAt: v.number()
  }),
  users: defineTable({
    email: v.string(),
    name: v.optional(v.string()),
    passwordHash: v.optional(v.string()),
    role: v.optional(v.string()),
    createdAt: v.number(),
    updatedAt: v.number()
  })
});
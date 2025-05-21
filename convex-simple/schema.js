// Basic schema for testing deployment
import { defineSchema, defineTable } from "convex/server";
import { v } from "convex/values";

export default defineSchema({
  // Simple molecules table for testing
  molecules: defineTable({
    name: v.string(),
    formula: v.optional(v.string()),
    createdAt: v.number(),
  })
    .index("by_name", ["name"])
});
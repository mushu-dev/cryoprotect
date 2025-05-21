/**
 * Define the schema for migrations in Convex.
 * 
 * This file provides schema definitions that are used by the Convex migration system.
 * Import this into your main schema.ts file.
 */

import { defineSchema, defineTable } from "convex/server";
import { v } from "convex/values";

/**
 * Schema for the migrations table
 */
export const migrations = defineTable({
  // Version number (e.g., "001", "002", etc.)
  version: v.string(),
  // Name of the migration
  name: v.string(),
  // Whether the migration has been successfully applied
  applied: v.boolean(),
  // When the migration was applied
  appliedAt: v.number(),
  // Any error information if the migration failed
  error: v.optional(v.string()),
  // Optional metadata for the migration
  metadata: v.optional(v.object({}))
})
  .index("by_version", ["version"])
  .index("by_applied", ["applied"]);

/**
 * Use this to include the migration schema in your main schema
 */
export function includeMigrationSchema(schema: any) {
  return {
    ...schema,
    migrations
  };
}
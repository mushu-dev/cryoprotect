/**
 * Convex migration tracking module.
 * 
 * This module provides functions for tracking applied migrations in Convex.
 */

import { mutation, query } from "../_generated/server";
import { v } from "convex/values";
import { Id } from "../_generated/dataModel";

/**
 * Schema for our migrations table
 */
export const migrationsTableSchema = {
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
};

/**
 * Get all applied migrations sorted by version
 */
export const getAppliedMigrations = query({
  handler: async (ctx) => {
    // Get all applied migrations sorted by version
    return await ctx.db
      .query("migrations")
      .filter((q) => q.eq(q.field("applied"), true))
      .order("asc", "version")
      .collect();
  },
});

/**
 * Get all migrations (applied and pending)
 */
export const getAllMigrations = query({
  handler: async (ctx) => {
    // Get all migrations sorted by version
    return await ctx.db
      .query("migrations")
      .order("asc", "version")
      .collect();
  },
});

/**
 * Check if a migration has been applied
 */
export const isMigrationApplied = query({
  args: { version: v.string() },
  handler: async (ctx, args) => {
    const migration = await ctx.db
      .query("migrations")
      .filter((q) => q.eq(q.field("version"), args.version))
      .first();
    
    return migration !== null && migration.applied;
  },
});

/**
 * Record a migration as applied
 */
export const recordMigration = mutation({
  args: { 
    version: v.string(),
    name: v.string(),
    metadata: v.optional(v.object({}))
  },
  handler: async (ctx, args) => {
    // Check if this migration already exists
    const existing = await ctx.db
      .query("migrations")
      .filter((q) => q.eq(q.field("version"), args.version))
      .first();
    
    const timestamp = Date.now();
    
    if (existing) {
      // Update the existing record
      return await ctx.db.patch(existing._id, {
        applied: true,
        appliedAt: timestamp,
        error: undefined,
        metadata: args.metadata
      });
    } else {
      // Create a new record
      return await ctx.db.insert("migrations", {
        version: args.version,
        name: args.name,
        applied: true,
        appliedAt: timestamp,
        metadata: args.metadata
      });
    }
  },
});

/**
 * Record a failed migration
 */
export const recordFailedMigration = mutation({
  args: { 
    version: v.string(),
    name: v.string(),
    error: v.string(),
    metadata: v.optional(v.object({}))
  },
  handler: async (ctx, args) => {
    // Check if this migration already exists
    const existing = await ctx.db
      .query("migrations")
      .filter((q) => q.eq(q.field("version"), args.version))
      .first();
    
    const timestamp = Date.now();
    
    if (existing) {
      // Update the existing record
      return await ctx.db.patch(existing._id, {
        applied: false,
        appliedAt: timestamp,
        error: args.error,
        metadata: args.metadata
      });
    } else {
      // Create a new record
      return await ctx.db.insert("migrations", {
        version: args.version,
        name: args.name,
        applied: false,
        appliedAt: timestamp,
        error: args.error,
        metadata: args.metadata
      });
    }
  },
});

/**
 * Remove a migration record (for rollback)
 */
export const removeMigrationRecord = mutation({
  args: { version: v.string() },
  handler: async (ctx, args) => {
    const existing = await ctx.db
      .query("migrations")
      .filter((q) => q.eq(q.field("version"), args.version))
      .first();
    
    if (existing) {
      await ctx.db.delete(existing._id);
      return true;
    }
    
    return false;
  },
});

/**
 * Initialize migration tracking
 * Makes sure the migrations table exists
 */
export const initializeMigrationTracking = mutation({
  handler: async (ctx) => {
    // In Convex, tables are created automatically when you insert data
    // We can just insert a dummy record and then delete it
    const id = await ctx.db.insert("migrations", {
      version: "000",
      name: "initialize_migration_tracking",
      applied: true,
      appliedAt: Date.now(),
    });
    
    await ctx.db.delete(id);
    
    return true;
  },
});
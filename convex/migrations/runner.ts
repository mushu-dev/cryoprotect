/**
 * Convex migration runner.
 * 
 * This module provides functions for applying and rolling back migrations in Convex.
 */

import { mutation, query, action } from "../_generated/server";
import { v } from "convex/values";
import { Id } from "../_generated/dataModel";
import { internal } from "../_generated/api";
import { ConvexError } from "convex/values";
import fs from "fs";
import path from "path";

/**
 * Get status of all migrations with comparison to filesystem
 */
export const getMigrationStatus = query({
  args: { migrationsDir: v.optional(v.string()) },
  handler: async (ctx, args) => {
    // Get all migrations from the database
    const dbMigrations = await ctx.db
      .query("migrations")
      .collect();
    
    // Create a mapping of known migrations
    const migrations: Record<string, any> = {};
    
    // Add database migrations to the mapping
    for (const migration of dbMigrations) {
      migrations[migration.version] = {
        version: migration.version,
        name: migration.name,
        applied: migration.applied,
        appliedAt: migration.appliedAt,
        error: migration.error,
        source: 'database',
        id: migration._id
      };
    }
    
    return { migrations };
  },
});

/**
 * Apply a single migration
 */
export const applyMigration = mutation({
  args: { 
    version: v.string(),
    name: v.string(),
    migrationJson: v.string()
  },
  handler: async (ctx, args) => {
    try {
      // Parse the migration actions
      const migration = JSON.parse(args.migrationJson);
      
      // Apply schema changes if present
      if (migration.schema) {
        for (const [tableName, fields] of Object.entries(migration.schema)) {
          // Unfortunately, we can't dynamically create types in the Convex runtime
          // For schema migrations, we'll need to do these in the schema.ts file
          // This is just a placeholder for tracking the migration
        }
      }
      
      // Apply data migrations if present
      if (migration.data) {
        for (const action of migration.data) {
          if (action.type === 'insert') {
            await ctx.db.insert(action.table, action.data);
          } else if (action.type === 'update') {
            // For updates, we need to find the records first
            const records = await ctx.db
              .query(action.table)
              .filter(q => {
                let query = q;
                for (const [field, value] of Object.entries(action.filter)) {
                  query = query.eq(q.field(field), value);
                }
                return query;
              })
              .collect();
            
            for (const record of records) {
              await ctx.db.patch(record._id, action.data);
            }
          } else if (action.type === 'delete') {
            // For deletes, we need to find the records first
            const records = await ctx.db
              .query(action.table)
              .filter(q => {
                let query = q;
                for (const [field, value] of Object.entries(action.filter)) {
                  query = query.eq(q.field(field), value);
                }
                return query;
              })
              .collect();
            
            for (const record of records) {
              await ctx.db.delete(record._id);
            }
          } else if (action.type === 'function') {
            // For custom functions, we would need to dispatch to the right function
            // This is a complex topic and would require more infrastructure
            throw new ConvexError("Custom function migrations not supported yet");
          }
        }
      }
      
      // Record the successful migration
      await ctx.db.insert("migrations", {
        version: args.version,
        name: args.name,
        applied: true,
        appliedAt: Date.now(),
        metadata: { migrationJson: args.migrationJson }
      });
      
      return {
        success: true,
        message: `Migration ${args.version} (${args.name}) applied successfully`
      };
    } catch (error) {
      // Record the failed migration
      await ctx.db.insert("migrations", {
        version: args.version,
        name: args.name,
        applied: false,
        appliedAt: Date.now(),
        error: error instanceof Error ? error.message : String(error),
        metadata: { migrationJson: args.migrationJson }
      });
      
      throw new ConvexError(`Failed to apply migration ${args.version}: ${error}`);
    }
  },
});

/**
 * Roll back a single migration
 */
export const rollbackMigration = mutation({
  args: { version: v.string() },
  handler: async (ctx, args) => {
    try {
      // Find the migration in the database
      const migration = await ctx.db
        .query("migrations")
        .filter(q => q.eq(q.field("version"), args.version))
        .first();
      
      if (!migration) {
        throw new ConvexError(`Migration ${args.version} not found`);
      }
      
      if (!migration.applied) {
        throw new ConvexError(`Migration ${args.version} is not applied`);
      }
      
      // Check if we have the migrationJson in metadata
      if (!migration.metadata || !migration.metadata.migrationJson) {
        throw new ConvexError(`Migration ${args.version} does not have rollback information`);
      }
      
      // Parse the migration actions
      const migrationData = JSON.parse(migration.metadata.migrationJson);
      
      // Apply rollback actions (This is reverse of the apply logic)
      if (migrationData.data) {
        // Process in reverse order
        for (const action of [...migrationData.data].reverse()) {
          if (action.type === 'insert') {
            // For inserts, we delete on rollback
            const records = await ctx.db
              .query(action.table)
              .filter(q => {
                let query = q;
                for (const [field, value] of Object.entries(action.data)) {
                  query = query.eq(q.field(field), value);
                }
                return query;
              })
              .collect();
            
            for (const record of records) {
              await ctx.db.delete(record._id);
            }
          } else if (action.type === 'update') {
            // For updates, we update with the original values on rollback
            if (!action.originalData) {
              throw new ConvexError(`Cannot rollback update for ${action.table}: missing originalData`);
            }
            
            const records = await ctx.db
              .query(action.table)
              .filter(q => {
                let query = q;
                for (const [field, value] of Object.entries(action.filter)) {
                  query = query.eq(q.field(field), value);
                }
                return query;
              })
              .collect();
            
            for (const record of records) {
              await ctx.db.patch(record._id, action.originalData);
            }
          } else if (action.type === 'delete') {
            // For deletes, we insert the original data on rollback
            if (!action.originalData) {
              throw new ConvexError(`Cannot rollback delete for ${action.table}: missing originalData`);
            }
            
            await ctx.db.insert(action.table, action.originalData);
          } else if (action.type === 'function') {
            throw new ConvexError("Custom function rollbacks not supported yet");
          }
        }
      }
      
      // Schema changes can't be automatically rolled back in Convex
      // You'd need to implement a new migration to roll back schema changes
      
      // Remove the migration record
      await ctx.db.delete(migration._id);
      
      return {
        success: true,
        message: `Migration ${args.version} (${migration.name}) rolled back successfully`
      };
    } catch (error) {
      throw new ConvexError(`Failed to roll back migration ${args.version}: ${error}`);
    }
  },
});

/**
 * Apply migrations up to target version
 */
export const applyMigrations = action({
  args: { 
    targetVersion: v.optional(v.string()),
    dryRun: v.optional(v.boolean()),
    migrationsDir: v.optional(v.string())
  },
  handler: async (ctx, args) => {
    const dryRun = args.dryRun ?? false;
    const results = [];
    
    // Validate args
    if (!dryRun) {
      // Initialize migration tracking
      await ctx.runMutation(internal.migrations.tracker.initializeMigrationTracking, {});
    }
    
    // Get applied migrations
    const appliedMigrations = await ctx.runQuery(internal.migrations.tracker.getAppliedMigrations, {});
    const appliedVersions = new Set(appliedMigrations.map(m => m.version));
    
    // Placeholder for migration files
    // In a real implementation, you'd need to load these from filesystem
    // Currently we simulate this with hardcoded migrations based on a pattern
    // In a real project, these would be loaded from actual files or from config
    
    // Sort migrations by version
    const migrationFiles = [
      { version: "001", name: "initial_schema", path: "simulated_path_1" },
      { version: "002", name: "add_team_data", path: "simulated_path_2" },
      { version: "003", name: "add_indexes", path: "simulated_path_3" }
    ];
    
    // Filter to pending migrations up to target version
    const pendingMigrations = migrationFiles.filter(m => {
      const isApplied = appliedVersions.has(m.version);
      const isInRange = !args.targetVersion || m.version <= args.targetVersion;
      return !isApplied && isInRange;
    });
    
    if (pendingMigrations.length === 0) {
      return { applied: results, message: "No pending migrations to apply" };
    }
    
    // Apply migrations
    for (const migration of pendingMigrations) {
      results.push({
        version: migration.version,
        name: migration.name,
        status: dryRun ? "would_apply" : "applying"
      });
      
      if (!dryRun) {
        try {
          // Each migration would be a JSON file with schema and data changes
          // Here, we mock the migration JSON for the example
          const migrationJson = JSON.stringify({
            schema: {
              // Schema changes would be defined here
            },
            data: [
              // Data migrations would be defined here
              // Example: { type: 'insert', table: 'users', data: { name: 'Admin' } }
            ]
          });
          
          const result = await ctx.runMutation(internal.migrations.runner.applyMigration, {
            version: migration.version,
            name: migration.name,
            migrationJson
          });
          
          results[results.length - 1].status = "applied";
        } catch (error) {
          results[results.length - 1].status = "failed";
          results[results.length - 1].error = error instanceof Error ? error.message : String(error);
          break;
        }
      }
    }
    
    return {
      applied: results,
      message: dryRun 
        ? `Would apply ${results.length} migrations`
        : `Applied ${results.filter(r => r.status === "applied").length} migrations`
    };
  },
});

/**
 * Roll back migrations down to target version
 */
export const rollbackMigrations = action({
  args: { 
    targetVersion: v.optional(v.string()),
    dryRun: v.optional(v.boolean())
  },
  handler: async (ctx, args) => {
    const dryRun = args.dryRun ?? false;
    const results = [];
    
    // Get applied migrations
    const appliedMigrations = await ctx.runQuery(internal.migrations.tracker.getAppliedMigrations, {});
    
    // Filter migrations to roll back
    const migrationsToRollback = appliedMigrations
      .filter(m => !args.targetVersion || m.version > args.targetVersion)
      .sort((a, b) => b.version.localeCompare(a.version));  // Sort in reverse order
    
    if (migrationsToRollback.length === 0) {
      return { rolledBack: results, message: "No migrations to roll back" };
    }
    
    // Roll back migrations
    for (const migration of migrationsToRollback) {
      results.push({
        version: migration.version,
        name: migration.name,
        status: dryRun ? "would_roll_back" : "rolling_back"
      });
      
      if (!dryRun) {
        try {
          await ctx.runMutation(internal.migrations.runner.rollbackMigration, {
            version: migration.version
          });
          
          results[results.length - 1].status = "rolled_back";
        } catch (error) {
          results[results.length - 1].status = "failed";
          results[results.length - 1].error = error instanceof Error ? error.message : String(error);
          break;
        }
      }
    }
    
    return {
      rolledBack: results,
      message: dryRun 
        ? `Would roll back ${results.length} migrations`
        : `Rolled back ${results.filter(r => r.status === "rolled_back").length} migrations`
    };
  },
});
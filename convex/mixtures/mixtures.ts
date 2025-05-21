/**
 * CRUD operations for mixtures
 */

import { mutation, query } from "../_generated/server";
import { v } from "convex/values";
import { Id } from "../_generated/dataModel";
import { 
  Mixture, 
  CreateMixtureInput, 
  UpdateMixtureInput,
  MixtureFilter,
  MixtureQueryOptions,
  MixtureWithComponents
} from "./types";
import { 
  validateCreateMixtureInput, 
  validateUpdateMixtureInput,
  validateMixtureExists,
  validateMixtureAccess
} from "./validation";
import { expandMixtureWithComponents } from "./helpers";
import { getCurrentUser } from "../auth/users";

/**
 * Create a new mixture
 */
export const createMixture = mutation({
  args: {
    mixture: v.object({
      name: v.string(),
      description: v.optional(v.string()),
      type: v.optional(v.string()),
      projectId: v.optional(v.id("projects")),
      public: v.optional(v.boolean())
    })
  },
  handler: async (ctx, args) => {
    // Validate input
    validateCreateMixtureInput(args.mixture);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Prepare data for insertion
    const now = Date.now();
    const mixtureData = {
      ...args.mixture,
      createdBy: userId,
      createdAt: now,
      updatedAt: now,
      public: args.mixture.public ?? false
    };
    
    // Insert mixture
    const mixtureId = await ctx.db.insert("mixtures", mixtureData);
    
    // Create audit log entry
    await ctx.db.insert("scientificDataAudit", {
      table: "mixtures",
      documentId: mixtureId,
      operation: "create",
      userId,
      newValue: mixtureData,
      timestamp: now,
    });
    
    return mixtureId;
  }
});

/**
 * Get a mixture by ID
 */
export const getMixture = query({
  args: {
    mixtureId: v.id("mixtures"),
    includeComponents: v.optional(v.boolean())
  },
  handler: async (ctx, args) => {
    // Validate mixture exists
    validateMixtureExists(args.mixtureId);
    
    // Get the mixture
    const mixture = await ctx.db.get(args.mixtureId);
    if (!mixture) {
      return null;
    }
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Check access
    validateMixtureAccess(
      args.mixtureId, 
      userId, 
      mixture.public, 
      mixture.createdBy
    );
    
    // Include components if requested
    if (args.includeComponents) {
      return expandMixtureWithComponents(ctx, mixture);
    }
    
    return mixture;
  }
});

/**
 * Update an existing mixture
 */
export const updateMixture = mutation({
  args: {
    mixtureId: v.id("mixtures"),
    update: v.object({
      name: v.optional(v.string()),
      description: v.optional(v.string()),
      type: v.optional(v.string()),
      projectId: v.optional(v.id("projects")),
      public: v.optional(v.boolean())
    })
  },
  handler: async (ctx, args) => {
    // Validate input
    validateUpdateMixtureInput(args.update);
    validateMixtureExists(args.mixtureId);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Get existing mixture
    const existingMixture = await ctx.db.get(args.mixtureId);
    if (!existingMixture) {
      throw new Error(`Mixture with ID ${args.mixtureId} not found`);
    }
    
    // Check access
    validateMixtureAccess(
      args.mixtureId, 
      userId, 
      existingMixture.public, 
      existingMixture.createdBy
    );
    
    // Prepare update data
    const now = Date.now();
    const updateData = {
      ...args.update,
      updatedAt: now
    };
    
    // Update mixture
    await ctx.db.patch(args.mixtureId, updateData);
    
    // Create audit log entry
    await ctx.db.insert("scientificDataAudit", {
      table: "mixtures",
      documentId: args.mixtureId,
      operation: "update",
      userId,
      previousValue: existingMixture,
      newValue: { ...existingMixture, ...updateData },
      timestamp: now
    });
    
    return args.mixtureId;
  }
});

/**
 * Delete a mixture and its components
 */
export const deleteMixture = mutation({
  args: {
    mixtureId: v.id("mixtures")
  },
  handler: async (ctx, args) => {
    // Validate mixture exists
    validateMixtureExists(args.mixtureId);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Get existing mixture
    const existingMixture = await ctx.db.get(args.mixtureId);
    if (!existingMixture) {
      throw new Error(`Mixture with ID ${args.mixtureId} not found`);
    }
    
    // Check access - only creator can delete
    if (!userId || !existingMixture.createdBy || !userId.equals(existingMixture.createdBy)) {
      throw new Error("Only the creator can delete a mixture");
    }
    
    // Delete all components first
    const components = await ctx.db
      .query("mixtureComponents")
      .withIndex("by_mixture", q => q.eq("mixtureId", args.mixtureId))
      .collect();
    
    // Record time for audit logs
    const now = Date.now();
    
    // Delete each component with audit log
    for (const component of components) {
      await ctx.db.insert("scientificDataAudit", {
        table: "mixtureComponents",
        documentId: component._id,
        operation: "delete",
        userId,
        previousValue: component,
        timestamp: now
      });
      
      await ctx.db.delete(component._id);
    }
    
    // Create audit log entry for mixture deletion
    await ctx.db.insert("scientificDataAudit", {
      table: "mixtures",
      documentId: args.mixtureId,
      operation: "delete",
      userId,
      previousValue: existingMixture,
      timestamp: now
    });
    
    // Delete the mixture
    await ctx.db.delete(args.mixtureId);
    
    return true;
  }
});

/**
 * List mixtures with optional filtering
 */
export const listMixtures = query({
  args: {
    filter: v.optional(v.object({
      name: v.optional(v.string()),
      type: v.optional(v.string()),
      projectId: v.optional(v.id("projects")),
      public: v.optional(v.boolean()),
      includeComponents: v.optional(v.boolean())
    })),
    options: v.optional(v.object({
      limit: v.optional(v.number()),
      cursor: v.optional(v.string())
    }))
  },
  handler: async (ctx, args) => {
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Set up the query
    let query = ctx.db.query("mixtures");
    
    // Apply filters
    if (args.filter) {
      if (args.filter.name) {
        // For a real implementation, we would use a full-text search or similar
        // For now, we'll do a simple filter based on name
        query = query.filter(q => 
          q.includes("name", args.filter!.name!)
        );
      }
      
      if (args.filter.type) {
        query = query.withIndex("by_type", q => 
          q.eq("type", args.filter!.type!)
        );
      }
      
      if (args.filter.projectId) {
        query = query.withIndex("by_project", q => 
          q.eq("projectId", args.filter!.projectId!)
        );
      }
    }
    
    // Filter based on access permissions
    if (!userId) {
      // If no user is logged in, only show public mixtures
      query = query.withIndex("by_public", q => q.eq("public", true));
    } else {
      // If a user is logged in, show public mixtures and their own
      query = query.filter(q => 
        q.or(
          q.eq(q.field("public"), true),
          q.eq(q.field("createdBy"), userId)
        )
      );
    }
    
    // Apply pagination
    if (args.options?.cursor) {
      query = query.withCursor(args.options.cursor);
    }
    
    if (args.options?.limit) {
      query = query.take(args.options.limit);
    } else {
      query = query.take(50); // Default limit
    }
    
    // Execute query
    const mixtures = await query.collect();
    
    // Include components if requested
    if (args.filter?.includeComponents) {
      const mixturesWithComponents: MixtureWithComponents[] = [];
      
      for (const mixture of mixtures) {
        mixturesWithComponents.push(await expandMixtureWithComponents(ctx, mixture));
      }
      
      return mixturesWithComponents;
    }
    
    return mixtures;
  }
});

/**
 * Search mixtures by name
 */
export const searchMixtures = query({
  args: {
    query: v.string(),
    limit: v.optional(v.number())
  },
  handler: async (ctx, args) => {
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Perform a simple search based on name
    // For a real implementation, we would use a more sophisticated search
    let query = ctx.db.query("mixtures")
      .filter(q => 
        q.includes("name", args.query)
      );
    
    // Filter based on access permissions
    if (!userId) {
      // If no user is logged in, only show public mixtures
      query = query.filter(q => q.eq(q.field("public"), true));
    } else {
      // If a user is logged in, show public mixtures and their own
      query = query.filter(q => 
        q.or(
          q.eq(q.field("public"), true),
          q.eq(q.field("createdBy"), userId)
        )
      );
    }
    
    // Apply limit
    const limit = args.limit || 20;
    query = query.take(limit);
    
    // Execute query
    return query.collect();
  }
});
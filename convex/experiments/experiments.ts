/**
 * CRUD operations for experiments
 */

import { mutation, query } from "../_generated/server";
import { v } from "convex/values";
import { Id } from "../_generated/dataModel";
import { 
  Experiment, 
  CreateExperimentInput, 
  UpdateExperimentInput,
  ExperimentFilter,
  ExperimentQueryOptions
} from "./types";
import { 
  validateCreateExperimentInput, 
  validateUpdateExperimentInput,
  validateExperimentExists,
  validateExperimentAccess
} from "./validation";
import { 
  expandExperimentWithResults, 
  expandExperimentWithMixture,
  expandExperimentComplete,
  experimentHasResults
} from "./helpers";
import { getCurrentUser } from "../auth/users";

/**
 * Create a new experiment
 */
export const createExperiment = mutation({
  args: {
    experiment: v.object({
      name: v.string(),
      description: v.optional(v.string()),
      mixtureId: v.optional(v.id("mixtures")),
      protocol: v.optional(v.string()),
      projectId: v.optional(v.id("projects")),
      date: v.optional(v.number()),
      status: v.optional(v.string()),
      public: v.optional(v.boolean())
    })
  },
  handler: async (ctx, args) => {
    // Validate input
    validateCreateExperimentInput(args.experiment);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // If mixtureId is provided, verify it exists
    if (args.experiment.mixtureId) {
      const mixture = await ctx.db.get(args.experiment.mixtureId);
      if (!mixture) {
        throw new Error(`Mixture with ID ${args.experiment.mixtureId} not found`);
      }
    }
    
    // Prepare data for insertion
    const now = Date.now();
    const experimentData = {
      ...args.experiment,
      status: args.experiment.status || "planned",
      conductedBy: userId,
      createdAt: now,
      updatedAt: now,
      public: args.experiment.public ?? false
    };
    
    // Insert experiment
    const experimentId = await ctx.db.insert("experiments", experimentData);
    
    // Create audit log entry
    await ctx.db.insert("scientificDataAudit", {
      table: "experiments",
      documentId: experimentId,
      operation: "create",
      userId,
      newValue: experimentData,
      timestamp: now,
    });
    
    return experimentId;
  }
});

/**
 * Get an experiment by ID
 */
export const getExperiment = query({
  args: {
    experimentId: v.id("experiments"),
    includeResults: v.optional(v.boolean()),
    includeMixture: v.optional(v.boolean())
  },
  handler: async (ctx, args) => {
    // Validate experiment exists
    validateExperimentExists(args.experimentId);
    
    // Get the experiment
    const experiment = await ctx.db.get(args.experimentId);
    if (!experiment) {
      return null;
    }
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Check access
    validateExperimentAccess(
      args.experimentId, 
      userId, 
      experiment.public, 
      experiment.conductedBy
    );
    
    // Include complete details if both flags are true
    if (args.includeResults && args.includeMixture) {
      return expandExperimentComplete(ctx, experiment);
    }
    
    // Include results if requested
    if (args.includeResults) {
      return expandExperimentWithResults(ctx, experiment);
    }
    
    // Include mixture if requested
    if (args.includeMixture) {
      return expandExperimentWithMixture(ctx, experiment);
    }
    
    return experiment;
  }
});

/**
 * Update an existing experiment
 */
export const updateExperiment = mutation({
  args: {
    experimentId: v.id("experiments"),
    update: v.object({
      name: v.optional(v.string()),
      description: v.optional(v.string()),
      mixtureId: v.optional(v.id("mixtures")),
      protocol: v.optional(v.string()),
      projectId: v.optional(v.id("projects")),
      date: v.optional(v.number()),
      status: v.optional(v.string()),
      public: v.optional(v.boolean())
    })
  },
  handler: async (ctx, args) => {
    // Validate input
    validateUpdateExperimentInput(args.update);
    validateExperimentExists(args.experimentId);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Get existing experiment
    const existingExperiment = await ctx.db.get(args.experimentId);
    if (!existingExperiment) {
      throw new Error(`Experiment with ID ${args.experimentId} not found`);
    }
    
    // Check access
    validateExperimentAccess(
      args.experimentId, 
      userId, 
      existingExperiment.public, 
      existingExperiment.conductedBy
    );
    
    // If mixtureId is provided, verify it exists
    if (args.update.mixtureId) {
      const mixture = await ctx.db.get(args.update.mixtureId);
      if (!mixture) {
        throw new Error(`Mixture with ID ${args.update.mixtureId} not found`);
      }
    }
    
    // Prepare update data
    const now = Date.now();
    const updateData = {
      ...args.update,
      updatedAt: now
    };
    
    // Update experiment
    await ctx.db.patch(args.experimentId, updateData);
    
    // Create audit log entry
    await ctx.db.insert("scientificDataAudit", {
      table: "experiments",
      documentId: args.experimentId,
      operation: "update",
      userId,
      previousValue: existingExperiment,
      newValue: { ...existingExperiment, ...updateData },
      timestamp: now
    });
    
    return args.experimentId;
  }
});

/**
 * Delete an experiment and its results
 */
export const deleteExperiment = mutation({
  args: {
    experimentId: v.id("experiments")
  },
  handler: async (ctx, args) => {
    // Validate experiment exists
    validateExperimentExists(args.experimentId);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Get existing experiment
    const existingExperiment = await ctx.db.get(args.experimentId);
    if (!existingExperiment) {
      throw new Error(`Experiment with ID ${args.experimentId} not found`);
    }
    
    // Check access - only creator can delete
    if (!userId || 
        !existingExperiment.conductedBy || 
        !userId.equals(existingExperiment.conductedBy)) {
      throw new Error("Only the creator can delete an experiment");
    }
    
    // Delete all results first
    const results = await ctx.db
      .query("experimentResults")
      .withIndex("by_experiment", q => q.eq("experimentId", args.experimentId))
      .collect();
    
    // Record time for audit logs
    const now = Date.now();
    
    // Delete each result with audit log
    for (const result of results) {
      await ctx.db.insert("scientificDataAudit", {
        table: "experimentResults",
        documentId: result._id,
        operation: "delete",
        userId,
        previousValue: result,
        timestamp: now
      });
      
      await ctx.db.delete(result._id);
    }
    
    // Create audit log entry for experiment deletion
    await ctx.db.insert("scientificDataAudit", {
      table: "experiments",
      documentId: args.experimentId,
      operation: "delete",
      userId,
      previousValue: existingExperiment,
      timestamp: now
    });
    
    // Delete the experiment
    await ctx.db.delete(args.experimentId);
    
    return true;
  }
});

/**
 * List experiments with optional filtering
 */
export const listExperiments = query({
  args: {
    filter: v.optional(v.object({
      name: v.optional(v.string()),
      mixtureId: v.optional(v.id("mixtures")),
      projectId: v.optional(v.id("projects")),
      status: v.optional(v.string()),
      dateRange: v.optional(v.object({
        start: v.optional(v.number()),
        end: v.optional(v.number())
      })),
      public: v.optional(v.boolean()),
      includeResults: v.optional(v.boolean()),
      includeMixture: v.optional(v.boolean())
    })),
    options: v.optional(v.object({
      limit: v.optional(v.number()),
      cursor: v.optional(v.string()),
      sortBy: v.optional(v.string()),
      sortDirection: v.optional(v.string())
    }))
  },
  handler: async (ctx, args) => {
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Set up the query
    let query = ctx.db.query("experiments");
    
    // Apply filters
    if (args.filter) {
      if (args.filter.name) {
        // For a real implementation, we would use a full-text search or similar
        // For now, we'll do a simple filter based on name
        query = query.filter(q => 
          q.includes("name", args.filter!.name!)
        );
      }
      
      if (args.filter.mixtureId) {
        query = query.withIndex("by_mixture", q => 
          q.eq("mixtureId", args.filter!.mixtureId!)
        );
      }
      
      if (args.filter.projectId) {
        query = query.withIndex("by_project", q => 
          q.eq("projectId", args.filter!.projectId!)
        );
      }
      
      if (args.filter.status) {
        query = query.withIndex("by_status", q => 
          q.eq("status", args.filter!.status!)
        );
      }
      
      if (args.filter.dateRange) {
        if (args.filter.dateRange.start !== undefined) {
          query = query.filter(q => 
            q.gte(q.field("date"), args.filter!.dateRange!.start!)
          );
        }
        
        if (args.filter.dateRange.end !== undefined) {
          query = query.filter(q => 
            q.lte(q.field("date"), args.filter!.dateRange!.end!)
          );
        }
      }
    }
    
    // Filter based on access permissions
    if (!userId) {
      // If no user is logged in, only show public experiments
      query = query.withIndex("by_public", q => q.eq("public", true));
    } else {
      // If a user is logged in, show public experiments and their own
      query = query.filter(q => 
        q.or(
          q.eq(q.field("public"), true),
          q.eq(q.field("conductedBy"), userId)
        )
      );
    }
    
    // Apply sorting
    if (args.options?.sortBy) {
      const sortDirection = args.options.sortDirection === "desc" ? "desc" : "asc";
      
      switch (args.options.sortBy) {
        case "name":
          query = query.order(sortDirection);
          break;
        case "date":
          query = query.order("date", sortDirection);
          break;
        case "status":
          query = query.order("status", sortDirection);
          break;
        case "updatedAt":
          query = query.order("updatedAt", sortDirection);
          break;
        default:
          query = query.order("updatedAt", "desc"); // default sort
      }
    } else {
      // Default sort by updatedAt descending
      query = query.order("updatedAt", "desc");
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
    const experiments = await query.collect();
    
    // Handle including additional data
    if (args.filter?.includeResults && args.filter?.includeMixture) {
      const experimentsWithDetails = [];
      
      for (const experiment of experiments) {
        experimentsWithDetails.push(await expandExperimentComplete(ctx, experiment));
      }
      
      return experimentsWithDetails;
    } else if (args.filter?.includeResults) {
      const experimentsWithResults = [];
      
      for (const experiment of experiments) {
        experimentsWithResults.push(await expandExperimentWithResults(ctx, experiment));
      }
      
      return experimentsWithResults;
    } else if (args.filter?.includeMixture) {
      const experimentsWithMixtures = [];
      
      for (const experiment of experiments) {
        experimentsWithMixtures.push(await expandExperimentWithMixture(ctx, experiment));
      }
      
      return experimentsWithMixtures;
    }
    
    return experiments;
  }
});

/**
 * Search experiments by name
 */
export const searchExperiments = query({
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
    let query = ctx.db.query("experiments")
      .filter(q => 
        q.includes("name", args.query)
      );
    
    // Filter based on access permissions
    if (!userId) {
      // If no user is logged in, only show public experiments
      query = query.filter(q => q.eq(q.field("public"), true));
    } else {
      // If a user is logged in, show public experiments and their own
      query = query.filter(q => 
        q.or(
          q.eq(q.field("public"), true),
          q.eq(q.field("conductedBy"), userId)
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

/**
 * Update experiment status
 */
export const updateExperimentStatus = mutation({
  args: {
    experimentId: v.id("experiments"),
    status: v.string()
  },
  handler: async (ctx, args) => {
    // Validate experiment exists
    validateExperimentExists(args.experimentId);
    
    // Validate status
    const validStatuses = ["planned", "in-progress", "completed", "failed"];
    if (!validStatuses.includes(args.status)) {
      throw new Error(
        `Status must be one of: ${validStatuses.join(", ")}`
      );
    }
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Get existing experiment
    const existingExperiment = await ctx.db.get(args.experimentId);
    if (!existingExperiment) {
      throw new Error(`Experiment with ID ${args.experimentId} not found`);
    }
    
    // Check access
    validateExperimentAccess(
      args.experimentId, 
      userId, 
      existingExperiment.public, 
      existingExperiment.conductedBy
    );
    
    // Prepare update data
    const now = Date.now();
    const updateData = {
      status: args.status,
      updatedAt: now
    };
    
    // Update experiment
    await ctx.db.patch(args.experimentId, updateData);
    
    // Create audit log entry
    await ctx.db.insert("scientificDataAudit", {
      table: "experiments",
      documentId: args.experimentId,
      operation: "update",
      userId,
      previousValue: existingExperiment,
      newValue: { ...existingExperiment, ...updateData },
      timestamp: now
    });
    
    return args.experimentId;
  }
});
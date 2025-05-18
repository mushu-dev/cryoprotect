/**
 * CRUD operations for enhanced experiments
 */

import { mutation, query } from "../_generated/server";
import { v } from "convex/values";
import { Id } from "../_generated/dataModel";
import { 
  EnhancedExperiment,
  CreateEnhancedExperimentInput,
  UpdateEnhancedExperimentInput,
  EnhancedExperimentFilter,
  EnhancedExperimentQueryOptions,
  EnhancedExperimentWithDetails
} from "./enhanced_types";
import { 
  validateCreateEnhancedExperimentInput,
  validateUpdateEnhancedExperimentInput,
  validateEnhancedExperimentExists,
  validateEnhancedExperimentAccess,
  validateProtocolExists,
  validateTissueTypeExists
} from "./enhanced_validation";
import { 
  expandEnhancedExperimentWithDetails,
  expandEnhancedExperimentWithResults,
  expandEnhancedExperimentWithProtocol,
  expandEnhancedExperimentWithMixture,
  expandEnhancedExperimentWithTissueTypes,
  expandEnhancedExperimentWithEquipment,
  expandEnhancedExperimentWithTimeSeries
} from "./enhanced_helpers";
import { getCurrentUser } from "../auth/users";

/**
 * Create a new enhanced experiment
 */
export const createEnhancedExperiment = mutation({
  args: {
    experiment: v.object({
      name: v.string(),
      description: v.optional(v.string()),
      experimentTypeId: v.optional(v.string()),
      protocolId: v.optional(v.id("protocols")),
      mixtureId: v.optional(v.id("mixtures")),
      temperature: v.optional(v.number()),
      temperatureUnit: v.optional(v.string()),
      coolingRate: v.optional(v.number()),
      coolingRateUnit: v.optional(v.string()),
      thawingRate: v.optional(v.number()),
      thawingRateUnit: v.optional(v.string()),
      pressure: v.optional(v.number()),
      pressureUnit: v.optional(v.string()),
      parameters: v.optional(v.map(v.string(), v.any())),
      version: v.optional(v.number()),
      provenance: v.optional(v.map(v.string(), v.any())),
      projectId: v.optional(v.id("projects")),
      date: v.optional(v.number()),
      status: v.optional(v.string()),
      tags: v.optional(v.array(v.string())),
      public: v.optional(v.boolean())
    })
  },
  handler: async (ctx, args) => {
    // Validate input
    validateCreateEnhancedExperimentInput(args.experiment);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // If protocolId is provided, verify it exists
    if (args.experiment.protocolId) {
      await validateProtocolExists(ctx.db, args.experiment.protocolId);
    }
    
    // If mixtureId is provided, verify it exists
    if (args.experiment.mixtureId) {
      const mixture = await ctx.db.get(args.experiment.mixtureId);
      if (!mixture) {
        throw new Error(`Mixture with ID ${args.experiment.mixtureId} not found`);
      }
    }
    
    // If projectId is provided, verify it exists
    if (args.experiment.projectId) {
      const project = await ctx.db.get(args.experiment.projectId);
      if (!project) {
        throw new Error(`Project with ID ${args.experiment.projectId} not found`);
      }
    }
    
    // Prepare data for insertion
    const now = Date.now();
    const experimentData = {
      ...args.experiment,
      status: args.experiment.status || "planned",
      version: args.experiment.version || 1,
      conductedBy: userId,
      createdBy: userId,
      createdAt: now,
      updatedAt: now,
      public: args.experiment.public ?? false
    };
    
    // Insert experiment
    const experimentId = await ctx.db.insert("enhancedExperiments", experimentData);
    
    // Create audit log entry
    await ctx.db.insert("scientificDataAudit", {
      table: "enhancedExperiments",
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
 * Get an enhanced experiment by ID
 */
export const getEnhancedExperiment = query({
  args: {
    experimentId: v.id("enhancedExperiments"),
    options: v.optional(v.object({
      includeResults: v.optional(v.boolean()),
      includeProtocol: v.optional(v.boolean()),
      includeMixture: v.optional(v.boolean()),
      includeTissueTypes: v.optional(v.boolean()),
      includeEquipment: v.optional(v.boolean()),
      includeTimeSeries: v.optional(v.boolean())
    }))
  },
  handler: async (ctx, args) => {
    // Get the experiment
    const experiment = await ctx.db.get(args.experimentId);
    if (!experiment) {
      return null;
    }
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Check access
    await validateEnhancedExperimentAccess(
      ctx.db,
      args.experimentId, 
      userId, 
      experiment.public, 
      experiment.conductedBy
    );
    
    // Check if complete details are requested
    if (args.options?.includeResults && 
        args.options?.includeProtocol && 
        args.options?.includeMixture &&
        args.options?.includeTissueTypes &&
        args.options?.includeEquipment &&
        args.options?.includeTimeSeries) {
      return expandEnhancedExperimentWithDetails(ctx, experiment);
    }
    
    // Include individual components as requested
    let result: any = { ...experiment };
    
    if (args.options?.includeResults) {
      result = await expandEnhancedExperimentWithResults(ctx, result);
    }
    
    if (args.options?.includeProtocol && experiment.protocolId) {
      result = await expandEnhancedExperimentWithProtocol(ctx, result);
    }
    
    if (args.options?.includeMixture && experiment.mixtureId) {
      result = await expandEnhancedExperimentWithMixture(ctx, result);
    }
    
    if (args.options?.includeTissueTypes) {
      result = await expandEnhancedExperimentWithTissueTypes(ctx, result);
    }
    
    if (args.options?.includeEquipment) {
      result = await expandEnhancedExperimentWithEquipment(ctx, result);
    }
    
    if (args.options?.includeTimeSeries) {
      result = await expandEnhancedExperimentWithTimeSeries(ctx, result);
    }
    
    return result;
  }
});

/**
 * Update an existing enhanced experiment
 */
export const updateEnhancedExperiment = mutation({
  args: {
    experimentId: v.id("enhancedExperiments"),
    update: v.object({
      name: v.optional(v.string()),
      description: v.optional(v.string()),
      experimentTypeId: v.optional(v.string()),
      protocolId: v.optional(v.id("protocols")),
      mixtureId: v.optional(v.id("mixtures")),
      temperature: v.optional(v.number()),
      temperatureUnit: v.optional(v.string()),
      coolingRate: v.optional(v.number()),
      coolingRateUnit: v.optional(v.string()),
      thawingRate: v.optional(v.number()),
      thawingRateUnit: v.optional(v.string()),
      pressure: v.optional(v.number()),
      pressureUnit: v.optional(v.string()),
      parameters: v.optional(v.map(v.string(), v.any())),
      version: v.optional(v.number()),
      provenance: v.optional(v.map(v.string(), v.any())),
      projectId: v.optional(v.id("projects")),
      date: v.optional(v.number()),
      status: v.optional(v.string()),
      tags: v.optional(v.array(v.string())),
      public: v.optional(v.boolean())
    })
  },
  handler: async (ctx, args) => {
    // Validate input
    validateUpdateEnhancedExperimentInput(args.update);
    await validateEnhancedExperimentExists(ctx.db, args.experimentId);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Get existing experiment
    const existingExperiment = await ctx.db.get(args.experimentId);
    if (!existingExperiment) {
      throw new Error(`Experiment with ID ${args.experimentId} not found`);
    }
    
    // Check access
    await validateEnhancedExperimentAccess(
      ctx.db,
      args.experimentId, 
      userId, 
      existingExperiment.public, 
      existingExperiment.conductedBy
    );
    
    // If protocolId is provided, verify it exists
    if (args.update.protocolId) {
      await validateProtocolExists(ctx.db, args.update.protocolId);
    }
    
    // If mixtureId is provided, verify it exists
    if (args.update.mixtureId) {
      const mixture = await ctx.db.get(args.update.mixtureId);
      if (!mixture) {
        throw new Error(`Mixture with ID ${args.update.mixtureId} not found`);
      }
    }
    
    // If projectId is provided, verify it exists
    if (args.update.projectId) {
      const project = await ctx.db.get(args.update.projectId);
      if (!project) {
        throw new Error(`Project with ID ${args.update.projectId} not found`);
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
      table: "enhancedExperiments",
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
 * Delete an enhanced experiment and its related data
 */
export const deleteEnhancedExperiment = mutation({
  args: {
    experimentId: v.id("enhancedExperiments")
  },
  handler: async (ctx, args) => {
    // Validate experiment exists
    await validateEnhancedExperimentExists(ctx.db, args.experimentId);
    
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
    
    // Record time for audit logs
    const now = Date.now();
    
    // Delete all related data with audit logs
    
    // 1. Delete experiment results
    const results = await ctx.db
      .query("enhancedExperimentResults")
      .withIndex("by_experiment", q => q.eq("experimentId", args.experimentId))
      .collect();
    
    for (const result of results) {
      await ctx.db.insert("scientificDataAudit", {
        table: "enhancedExperimentResults",
        documentId: result._id,
        operation: "delete",
        userId,
        previousValue: result,
        timestamp: now
      });
      
      await ctx.db.delete(result._id);
    }
    
    // 2. Delete experiment equipment links
    const equipmentLinks = await ctx.db
      .query("experimentEquipment")
      .withIndex("by_experiment", q => q.eq("experimentId", args.experimentId))
      .collect();
    
    for (const link of equipmentLinks) {
      await ctx.db.insert("scientificDataAudit", {
        table: "experimentEquipment",
        documentId: link._id,
        operation: "delete",
        userId,
        previousValue: link,
        timestamp: now
      });
      
      await ctx.db.delete(link._id);
    }
    
    // 3. Delete time series and their data points
    const timeSeries = await ctx.db
      .query("timeSeries")
      .withIndex("by_experiment", q => q.eq("experimentId", args.experimentId))
      .collect();
    
    for (const series of timeSeries) {
      // Delete data points for this series
      const dataPoints = await ctx.db
        .query("timeSeriesData")
        .withIndex("by_series", q => q.eq("timeSeriesId", series._id))
        .collect();
      
      for (const point of dataPoints) {
        await ctx.db.insert("scientificDataAudit", {
          table: "timeSeriesData",
          documentId: point._id,
          operation: "delete",
          userId,
          previousValue: point,
          timestamp: now
        });
        
        await ctx.db.delete(point._id);
      }
      
      // Delete the time series itself
      await ctx.db.insert("scientificDataAudit", {
        table: "timeSeries",
        documentId: series._id,
        operation: "delete",
        userId,
        previousValue: series,
        timestamp: now
      });
      
      await ctx.db.delete(series._id);
    }
    
    // 4. Create audit log entry for experiment deletion
    await ctx.db.insert("scientificDataAudit", {
      table: "enhancedExperiments",
      documentId: args.experimentId,
      operation: "delete",
      userId,
      previousValue: existingExperiment,
      timestamp: now
    });
    
    // 5. Delete the experiment
    await ctx.db.delete(args.experimentId);
    
    return true;
  }
});

/**
 * List enhanced experiments with optional filtering
 */
export const listEnhancedExperiments = query({
  args: {
    filter: v.optional(v.object({
      name: v.optional(v.string()),
      experimentTypeId: v.optional(v.string()),
      protocolId: v.optional(v.id("protocols")),
      mixtureId: v.optional(v.id("mixtures")),
      conductedBy: v.optional(v.id("users")),
      projectId: v.optional(v.id("projects")),
      status: v.optional(v.string()),
      dateRange: v.optional(v.object({
        start: v.optional(v.number()),
        end: v.optional(v.number())
      })),
      tags: v.optional(v.array(v.string())),
      public: v.optional(v.boolean()),
      tissueTypeId: v.optional(v.id("tissueTypes"))
    })),
    options: v.optional(v.object({
      limit: v.optional(v.number()),
      cursor: v.optional(v.string()),
      includeResults: v.optional(v.boolean()),
      includeProtocol: v.optional(v.boolean()),
      includeMixture: v.optional(v.boolean()),
      includeTissueTypes: v.optional(v.boolean()),
      includeEquipment: v.optional(v.boolean()),
      includeTimeSeries: v.optional(v.boolean()),
      sortBy: v.optional(v.string()),
      sortDirection: v.optional(v.string())
    }))
  },
  handler: async (ctx, args) => {
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Set up the query
    let query = ctx.db.query("enhancedExperiments");
    
    // Apply filters
    if (args.filter) {
      if (args.filter.name) {
        // For a real implementation, we would use a full-text search or similar
        // For now, we'll do a simple filter based on name
        query = query.filter(q => 
          q.includes("name", args.filter!.name!)
        );
      }
      
      if (args.filter.experimentTypeId) {
        query = query.filter(q => 
          q.eq(q.field("experimentTypeId"), args.filter!.experimentTypeId!)
        );
      }
      
      if (args.filter.protocolId) {
        query = query.withIndex("by_protocol", q => 
          q.eq("protocolId", args.filter!.protocolId!)
        );
      }
      
      if (args.filter.mixtureId) {
        query = query.withIndex("by_mixture", q => 
          q.eq("mixtureId", args.filter!.mixtureId!)
        );
      }
      
      if (args.filter.conductedBy) {
        query = query.withIndex("by_user", q => 
          q.eq("conductedBy", args.filter!.conductedBy!)
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
        query = query.withIndex("by_date");
        
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
      
      if (args.filter.tags && args.filter.tags.length > 0) {
        // Filter experiments that have at least one of the provided tags
        query = query.filter(q => 
          q.includes(q.field("tags"), args.filter!.tags![0])
        );
        
        // Add additional tag filters if needed
        for (let i = 1; i < args.filter.tags.length; i++) {
          query = query.filter(q => 
            q.includes(q.field("tags"), args.filter!.tags![i])
          );
        }
      }
      
      // Filter by tissue type ID requires a join with experiment results
      if (args.filter.tissueTypeId) {
        // Get experiments that have results with this tissue type
        const experimentsWithTissueType = await ctx.db
          .query("enhancedExperimentResults")
          .withIndex("by_tissue", q => 
            q.eq("tissueTypeId", args.filter!.tissueTypeId!)
          )
          .collect();
        
        // Extract unique experiment IDs
        const experimentIds = [...new Set(
          experimentsWithTissueType.map(result => result.experimentId)
        )];
        
        if (experimentIds.length > 0) {
          // Filter to only include these experiments
          query = query.filter(q => 
            q.or(...experimentIds.map(id => q.eq(q.field("_id"), id)))
          );
        } else {
          // No experiments with this tissue type, return empty result
          return [];
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
          query = query.order("name", sortDirection);
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
    
    // If no additional details are requested, return the experiments as is
    if (!args.options?.includeResults && 
        !args.options?.includeProtocol && 
        !args.options?.includeMixture &&
        !args.options?.includeTissueTypes &&
        !args.options?.includeEquipment &&
        !args.options?.includeTimeSeries) {
      return experiments;
    }
    
    // Handle including additional data
    const experimentsWithDetails = [];
    
    for (const experiment of experiments) {
      let result: any = { ...experiment };
      
      if (args.options?.includeResults) {
        result = await expandEnhancedExperimentWithResults(ctx, result);
      }
      
      if (args.options?.includeProtocol && experiment.protocolId) {
        result = await expandEnhancedExperimentWithProtocol(ctx, result);
      }
      
      if (args.options?.includeMixture && experiment.mixtureId) {
        result = await expandEnhancedExperimentWithMixture(ctx, result);
      }
      
      if (args.options?.includeTissueTypes) {
        result = await expandEnhancedExperimentWithTissueTypes(ctx, result);
      }
      
      if (args.options?.includeEquipment) {
        result = await expandEnhancedExperimentWithEquipment(ctx, result);
      }
      
      if (args.options?.includeTimeSeries) {
        result = await expandEnhancedExperimentWithTimeSeries(ctx, result);
      }
      
      experimentsWithDetails.push(result);
    }
    
    return experimentsWithDetails;
  }
});

/**
 * Search enhanced experiments by name and tags
 */
export const searchEnhancedExperiments = query({
  args: {
    query: v.string(),
    limit: v.optional(v.number())
  },
  handler: async (ctx, args) => {
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Perform a simple search based on name or tags
    let query = ctx.db.query("enhancedExperiments")
      .filter(q => 
        q.includes("name", args.query)
        // Note: In a real system, we would also search in tags, descriptions, etc.
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
 * Update enhanced experiment status
 */
export const updateEnhancedExperimentStatus = mutation({
  args: {
    experimentId: v.id("enhancedExperiments"),
    status: v.string()
  },
  handler: async (ctx, args) => {
    // Validate experiment exists
    await validateEnhancedExperimentExists(ctx.db, args.experimentId);
    
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
    await validateEnhancedExperimentAccess(
      ctx.db,
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
      table: "enhancedExperiments",
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
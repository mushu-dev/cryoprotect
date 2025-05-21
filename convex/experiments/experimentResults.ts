/**
 * CRUD operations for experiment results
 */

import { mutation, query } from "../_generated/server";
import { v } from "convex/values";
import { Id } from "../_generated/dataModel";
import { 
  ExperimentResult, 
  CreateExperimentResultInput, 
  UpdateExperimentResultInput,
  ExperimentResultFilter,
  ExperimentResultQueryOptions 
} from "./types";
import { 
  validateCreateExperimentResultInput, 
  validateUpdateExperimentResultInput,
  validateExperimentResultExists,
  validateExperimentExists,
  validateExperimentAccess
} from "./validation";
import { getResultStatistics, groupResultsByParameter } from "./helpers";
import { getCurrentUser } from "../auth/users";

/**
 * Create a new experiment result
 */
export const createExperimentResult = mutation({
  args: {
    result: v.object({
      experimentId: v.id("experiments"),
      parameterName: v.string(),
      value: v.union(v.string(), v.number(), v.boolean(), v.null()),
      numericValue: v.optional(v.number()),
      units: v.optional(v.string()),
      notes: v.optional(v.string())
    })
  },
  handler: async (ctx, args) => {
    // Validate input
    validateCreateExperimentResultInput(args.result);
    validateExperimentExists(args.result.experimentId);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Get the experiment to check access
    const experiment = await ctx.db.get(args.result.experimentId);
    if (!experiment) {
      throw new Error(`Experiment with ID ${args.result.experimentId} not found`);
    }
    
    // Check access - only the conductor can add results
    if (!userId || 
        !experiment.conductedBy || 
        !userId.equals(experiment.conductedBy)) {
      throw new Error("Only the experiment conductor can add results");
    }
    
    // Ensure numeric values are properly stored
    let numericValue = args.result.numericValue;
    
    // If value is a number and numericValue not provided, use value as numericValue
    if (typeof args.result.value === "number" && numericValue === undefined) {
      numericValue = args.result.value;
    }
    
    // Prepare data for insertion
    const now = Date.now();
    const resultData = {
      ...args.result,
      numericValue,
      createdAt: now,
      updatedAt: now
    };
    
    // Insert result
    const resultId = await ctx.db.insert("experimentResults", resultData);
    
    // Create audit log entry
    await ctx.db.insert("scientificDataAudit", {
      table: "experimentResults",
      documentId: resultId,
      operation: "create",
      userId,
      newValue: resultData,
      timestamp: now
    });
    
    // Update the experiment's updatedAt field
    await ctx.db.patch(args.result.experimentId, { updatedAt: now });
    
    // If experiment is in "planned" status, update to "in-progress"
    if (experiment.status === "planned") {
      await ctx.db.patch(args.result.experimentId, { 
        status: "in-progress",
        updatedAt: now
      });
    }
    
    return resultId;
  }
});

/**
 * Get an experiment result by ID
 */
export const getExperimentResult = query({
  args: {
    resultId: v.id("experimentResults")
  },
  handler: async (ctx, args) => {
    // Validate result exists
    validateExperimentResultExists(args.resultId);
    
    // Get the result
    const result = await ctx.db.get(args.resultId);
    if (!result) {
      return null;
    }
    
    // Get the experiment to check access
    const experiment = await ctx.db.get(result.experimentId);
    if (!experiment) {
      return null;
    }
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Check access
    validateExperimentAccess(
      result.experimentId, 
      userId, 
      experiment.public, 
      experiment.conductedBy
    );
    
    return result;
  }
});

/**
 * Update an existing experiment result
 */
export const updateExperimentResult = mutation({
  args: {
    resultId: v.id("experimentResults"),
    update: v.object({
      parameterName: v.optional(v.string()),
      value: v.optional(v.union(v.string(), v.number(), v.boolean(), v.null())),
      numericValue: v.optional(v.number()),
      units: v.optional(v.string()),
      notes: v.optional(v.string())
    })
  },
  handler: async (ctx, args) => {
    // Validate input
    validateUpdateExperimentResultInput(args.update);
    validateExperimentResultExists(args.resultId);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Get existing result
    const existingResult = await ctx.db.get(args.resultId);
    if (!existingResult) {
      throw new Error(`Experiment result with ID ${args.resultId} not found`);
    }
    
    // Get the experiment to check access
    const experiment = await ctx.db.get(existingResult.experimentId);
    if (!experiment) {
      throw new Error(`Experiment with ID ${existingResult.experimentId} not found`);
    }
    
    // Check access - only the conductor can update results
    if (!userId || 
        !experiment.conductedBy || 
        !userId.equals(experiment.conductedBy)) {
      throw new Error("Only the experiment conductor can update results");
    }
    
    // Ensure numeric values are properly stored
    let updateData = { ...args.update };
    
    // If value is a number and numericValue not provided, use value as numericValue
    if (typeof args.update.value === "number" && args.update.numericValue === undefined) {
      updateData.numericValue = args.update.value;
    }
    
    // Prepare update data
    const now = Date.now();
    updateData = {
      ...updateData,
      updatedAt: now
    };
    
    // Update result
    await ctx.db.patch(args.resultId, updateData);
    
    // Create audit log entry
    await ctx.db.insert("scientificDataAudit", {
      table: "experimentResults",
      documentId: args.resultId,
      operation: "update",
      userId,
      previousValue: existingResult,
      newValue: { ...existingResult, ...updateData },
      timestamp: now
    });
    
    // Update the experiment's updatedAt field
    await ctx.db.patch(existingResult.experimentId, { updatedAt: now });
    
    return args.resultId;
  }
});

/**
 * Delete an experiment result
 */
export const deleteExperimentResult = mutation({
  args: {
    resultId: v.id("experimentResults")
  },
  handler: async (ctx, args) => {
    // Validate result exists
    validateExperimentResultExists(args.resultId);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Get existing result
    const existingResult = await ctx.db.get(args.resultId);
    if (!existingResult) {
      throw new Error(`Experiment result with ID ${args.resultId} not found`);
    }
    
    // Get the experiment to check access
    const experiment = await ctx.db.get(existingResult.experimentId);
    if (!experiment) {
      throw new Error(`Experiment with ID ${existingResult.experimentId} not found`);
    }
    
    // Check access - only the conductor can delete results
    if (!userId || 
        !experiment.conductedBy || 
        !userId.equals(experiment.conductedBy)) {
      throw new Error("Only the experiment conductor can delete results");
    }
    
    // Create audit log entry
    const now = Date.now();
    await ctx.db.insert("scientificDataAudit", {
      table: "experimentResults",
      documentId: args.resultId,
      operation: "delete",
      userId,
      previousValue: existingResult,
      timestamp: now
    });
    
    // Delete the result
    await ctx.db.delete(args.resultId);
    
    // Update the experiment's updatedAt field
    await ctx.db.patch(existingResult.experimentId, { updatedAt: now });
    
    return true;
  }
});

/**
 * List results for an experiment
 */
export const listExperimentResults = query({
  args: {
    experimentId: v.id("experiments"),
    options: v.optional(v.object({
      parameterName: v.optional(v.string()),
      limit: v.optional(v.number()),
      cursor: v.optional(v.string())
    }))
  },
  handler: async (ctx, args) => {
    // Validate experiment exists
    validateExperimentExists(args.experimentId);
    
    // Get the experiment to check access
    const experiment = await ctx.db.get(args.experimentId);
    if (!experiment) {
      throw new Error(`Experiment with ID ${args.experimentId} not found`);
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
    
    // Set up the query
    let query = ctx.db.query("experimentResults")
      .withIndex("by_experiment", q => q.eq("experimentId", args.experimentId));
    
    // Filter by parameter name if provided
    if (args.options?.parameterName) {
      query = query.filter(q => 
        q.eq(q.field("parameterName"), args.options!.parameterName!)
      );
    }
    
    // Apply pagination
    if (args.options?.cursor) {
      query = query.withCursor(args.options.cursor);
    }
    
    if (args.options?.limit) {
      query = query.take(args.options.limit);
    } else {
      query = query.take(100); // Higher default limit for results
    }
    
    // Execute query
    return query.collect();
  }
});

/**
 * Get statistics for experiment results
 */
export const getExperimentResultStats = query({
  args: {
    experimentId: v.id("experiments"),
    parameterName: v.optional(v.string())
  },
  handler: async (ctx, args) => {
    // Validate experiment exists
    validateExperimentExists(args.experimentId);
    
    // Get the experiment to check access
    const experiment = await ctx.db.get(args.experimentId);
    if (!experiment) {
      throw new Error(`Experiment with ID ${args.experimentId} not found`);
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
    
    // Get results
    let query = ctx.db.query("experimentResults")
      .withIndex("by_experiment", q => q.eq("experimentId", args.experimentId));
    
    // Filter by parameter name if provided
    if (args.parameterName) {
      query = query.filter(q => 
        q.eq(q.field("parameterName"), args.parameterName)
      );
    }
    
    const results = await query.collect();
    
    // If no parameter name provided, group results by parameter
    if (!args.parameterName) {
      const groupedResults = groupResultsByParameter(results);
      const stats: Record<string, any> = {};
      
      // Calculate statistics for each parameter group
      for (const [paramName, paramResults] of Object.entries(groupedResults)) {
        const paramStats = getResultStatistics(paramResults);
        if (paramStats) {
          stats[paramName] = paramStats;
        }
      }
      
      return stats;
    }
    
    // Calculate statistics for the specific parameter
    return getResultStatistics(results) || null;
  }
});

/**
 * Batch add multiple results to an experiment
 */
export const batchAddExperimentResults = mutation({
  args: {
    experimentId: v.id("experiments"),
    results: v.array(
      v.object({
        parameterName: v.string(),
        value: v.union(v.string(), v.number(), v.boolean(), v.null()),
        numericValue: v.optional(v.number()),
        units: v.optional(v.string()),
        notes: v.optional(v.string())
      })
    )
  },
  handler: async (ctx, args) => {
    // Validate experiment exists
    validateExperimentExists(args.experimentId);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Get the experiment to check access
    const experiment = await ctx.db.get(args.experimentId);
    if (!experiment) {
      throw new Error(`Experiment with ID ${args.experimentId} not found`);
    }
    
    // Check access - only the conductor can add results
    if (!userId || 
        !experiment.conductedBy || 
        !userId.equals(experiment.conductedBy)) {
      throw new Error("Only the experiment conductor can add results");
    }
    
    // Record time for all operations
    const now = Date.now();
    
    // Validate and insert each result
    const resultIds: Id<"experimentResults">[] = [];
    
    for (const result of args.results) {
      // Validate result input
      validateCreateExperimentResultInput({
        ...result,
        experimentId: args.experimentId
      });
      
      // Ensure numeric values are properly stored
      let numericValue = result.numericValue;
      
      // If value is a number and numericValue not provided, use value as numericValue
      if (typeof result.value === "number" && numericValue === undefined) {
        numericValue = result.value;
      }
      
      // Prepare result data
      const resultData = {
        ...result,
        experimentId: args.experimentId,
        numericValue,
        createdAt: now,
        updatedAt: now
      };
      
      // Insert result
      const resultId = await ctx.db.insert("experimentResults", resultData);
      
      // Create audit log entry
      await ctx.db.insert("scientificDataAudit", {
        table: "experimentResults",
        documentId: resultId,
        operation: "create",
        userId,
        newValue: resultData,
        timestamp: now
      });
      
      resultIds.push(resultId);
    }
    
    // Update the experiment's updatedAt field
    await ctx.db.patch(args.experimentId, { updatedAt: now });
    
    // If experiment is in "planned" status, update to "in-progress"
    if (experiment.status === "planned") {
      await ctx.db.patch(args.experimentId, { 
        status: "in-progress",
        updatedAt: now
      });
    }
    
    return resultIds;
  }
});
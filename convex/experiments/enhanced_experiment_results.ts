/**
 * CRUD operations for enhanced experiment results
 */

import { mutation, query } from "../_generated/server";
import { v } from "convex/values";
import { Id } from "../_generated/dataModel";
import { 
  EnhancedExperimentResult,
  CreateEnhancedExperimentResultInput,
  UpdateEnhancedExperimentResultInput
} from "./enhanced_types";
import { 
  validateCreateEnhancedExperimentResultInput,
  validateUpdateEnhancedExperimentResultInput,
  validateEnhancedExperimentExists,
  validateEnhancedExperimentAccess,
  validateTissueTypeExists
} from "./enhanced_validation";
import { getCurrentUser } from "../auth/users";

/**
 * Create a new enhanced experiment result
 */
export const createEnhancedExperimentResult = mutation({
  args: {
    result: v.object({
      experimentId: v.id("enhancedExperiments"),
      moleculeId: v.optional(v.id("molecules")),
      mixtureId: v.optional(v.id("mixtures")),
      tissueTypeId: v.id("tissueTypes"),
      parameterName: v.string(),
      value: v.union(v.string(), v.number(), v.boolean(), v.null()),
      numericValue: v.optional(v.number()),
      units: v.optional(v.string()),
      uncertainty: v.optional(v.object({
        type: v.string(),
        value: v.union(v.number(), v.array(v.number())),
        confidence: v.optional(v.number())
      })),
      provenance: v.optional(v.object({
        method: v.string(),
        reference: v.optional(v.string()),
        timestamp: v.number(),
        operator: v.optional(v.string())
      })),
      notes: v.optional(v.string())
    })
  },
  handler: async (ctx, args) => {
    // Validate input
    validateCreateEnhancedExperimentResultInput(args.result);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Verify experiment exists
    await validateEnhancedExperimentExists(ctx.db, args.result.experimentId);
    
    // Verify tissue type exists
    await validateTissueTypeExists(ctx.db, args.result.tissueTypeId);
    
    // Get experiment
    const experiment = await ctx.db.get(args.result.experimentId);
    if (!experiment) {
      throw new Error(`Experiment with ID ${args.result.experimentId} not found`);
    }
    
    // Check access to experiment
    await validateEnhancedExperimentAccess(
      ctx.db,
      args.result.experimentId,
      userId,
      experiment.public,
      experiment.conductedBy
    );
    
    // If molecule ID is provided, verify it exists
    if (args.result.moleculeId) {
      const molecule = await ctx.db.get(args.result.moleculeId);
      if (!molecule) {
        throw new Error(`Molecule with ID ${args.result.moleculeId} not found`);
      }
    }
    
    // If mixture ID is provided, verify it exists
    if (args.result.mixtureId) {
      const mixture = await ctx.db.get(args.result.mixtureId);
      if (!mixture) {
        throw new Error(`Mixture with ID ${args.result.mixtureId} not found`);
      }
    }
    
    // Prepare data for insertion
    const now = Date.now();
    const resultData = {
      ...args.result,
      createdBy: userId,
      createdAt: now,
      updatedAt: now
    };
    
    // Insert result
    const resultId = await ctx.db.insert("enhancedExperimentResults", resultData);
    
    // Create audit log entry
    await ctx.db.insert("scientificDataAudit", {
      table: "enhancedExperimentResults",
      documentId: resultId,
      operation: "create",
      userId,
      newValue: resultData,
      timestamp: now,
    });
    
    return resultId;
  }
});

/**
 * Get an enhanced experiment result by ID
 */
export const getEnhancedExperimentResult = query({
  args: {
    resultId: v.id("enhancedExperimentResults")
  },
  handler: async (ctx, args) => {
    // Get the result
    const result = await ctx.db.get(args.resultId);
    if (!result) {
      return null;
    }
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Check access to the experiment
    const experiment = await ctx.db.get(result.experimentId);
    if (!experiment) {
      throw new Error(`Experiment with ID ${result.experimentId} not found`);
    }
    
    await validateEnhancedExperimentAccess(
      ctx.db,
      result.experimentId,
      userId,
      experiment.public,
      experiment.conductedBy
    );
    
    return result;
  }
});

/**
 * Update an existing enhanced experiment result
 */
export const updateEnhancedExperimentResult = mutation({
  args: {
    resultId: v.id("enhancedExperimentResults"),
    update: v.object({
      parameterName: v.optional(v.string()),
      value: v.optional(v.union(v.string(), v.number(), v.boolean(), v.null())),
      numericValue: v.optional(v.number()),
      units: v.optional(v.string()),
      uncertainty: v.optional(v.object({
        type: v.string(),
        value: v.union(v.number(), v.array(v.number())),
        confidence: v.optional(v.number())
      })),
      provenance: v.optional(v.object({
        method: v.string(),
        reference: v.optional(v.string()),
        timestamp: v.number(),
        operator: v.optional(v.string())
      })),
      notes: v.optional(v.string())
    })
  },
  handler: async (ctx, args) => {
    // Validate input
    validateUpdateEnhancedExperimentResultInput(args.update);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Get existing result
    const existingResult = await ctx.db.get(args.resultId);
    if (!existingResult) {
      throw new Error(`Result with ID ${args.resultId} not found`);
    }
    
    // Get experiment
    const experiment = await ctx.db.get(existingResult.experimentId);
    if (!experiment) {
      throw new Error(`Experiment with ID ${existingResult.experimentId} not found`);
    }
    
    // Check access to experiment
    await validateEnhancedExperimentAccess(
      ctx.db,
      existingResult.experimentId,
      userId,
      experiment.public,
      experiment.conductedBy
    );
    
    // Prepare update data
    const now = Date.now();
    const updateData = {
      ...args.update,
      updatedAt: now
    };
    
    // Update result
    await ctx.db.patch(args.resultId, updateData);
    
    // Create audit log entry
    await ctx.db.insert("scientificDataAudit", {
      table: "enhancedExperimentResults",
      documentId: args.resultId,
      operation: "update",
      userId,
      previousValue: existingResult,
      newValue: { ...existingResult, ...updateData },
      timestamp: now
    });
    
    return args.resultId;
  }
});

/**
 * Delete an enhanced experiment result
 */
export const deleteEnhancedExperimentResult = mutation({
  args: {
    resultId: v.id("enhancedExperimentResults")
  },
  handler: async (ctx, args) => {
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Get existing result
    const existingResult = await ctx.db.get(args.resultId);
    if (!existingResult) {
      throw new Error(`Result with ID ${args.resultId} not found`);
    }
    
    // Get experiment
    const experiment = await ctx.db.get(existingResult.experimentId);
    if (!experiment) {
      throw new Error(`Experiment with ID ${existingResult.experimentId} not found`);
    }
    
    // Check access to experiment
    await validateEnhancedExperimentAccess(
      ctx.db,
      existingResult.experimentId,
      userId,
      experiment.public,
      experiment.conductedBy
    );
    
    // Record time for audit log
    const now = Date.now();
    
    // Create audit log entry
    await ctx.db.insert("scientificDataAudit", {
      table: "enhancedExperimentResults",
      documentId: args.resultId,
      operation: "delete",
      userId,
      previousValue: existingResult,
      timestamp: now
    });
    
    // Delete the result
    await ctx.db.delete(args.resultId);
    
    return true;
  }
});

/**
 * List results for an experiment
 */
export const listEnhancedExperimentResults = query({
  args: {
    experimentId: v.id("enhancedExperiments"),
    options: v.optional(v.object({
      parameterName: v.optional(v.string()),
      moleculeId: v.optional(v.id("molecules")),
      mixtureId: v.optional(v.id("mixtures")),
      tissueTypeId: v.optional(v.id("tissueTypes"))
    }))
  },
  handler: async (ctx, args) => {
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Get experiment
    const experiment = await ctx.db.get(args.experimentId);
    if (!experiment) {
      throw new Error(`Experiment with ID ${args.experimentId} not found`);
    }
    
    // Check access to experiment
    await validateEnhancedExperimentAccess(
      ctx.db,
      args.experimentId,
      userId,
      experiment.public,
      experiment.conductedBy
    );
    
    // Set up query
    let query = ctx.db
      .query("enhancedExperimentResults")
      .withIndex("by_experiment", q => q.eq("experimentId", args.experimentId));
    
    // Apply filters
    if (args.options) {
      if (args.options.parameterName) {
        query = query.filter(q => 
          q.eq(q.field("parameterName"), args.options!.parameterName!)
        );
      }
      
      if (args.options.moleculeId) {
        query = query.filter(q => 
          q.eq(q.field("moleculeId"), args.options!.moleculeId!)
        );
      }
      
      if (args.options.mixtureId) {
        query = query.filter(q => 
          q.eq(q.field("mixtureId"), args.options!.mixtureId!)
        );
      }
      
      if (args.options.tissueTypeId) {
        query = query.filter(q => 
          q.eq(q.field("tissueTypeId"), args.options!.tissueTypeId!)
        );
      }
    }
    
    // Execute query
    return query.collect();
  }
});
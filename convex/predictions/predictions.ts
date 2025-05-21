/**
 * CRUD operations for predictions
 */

import { mutation, query } from "../_generated/server";
import { v } from "convex/values";
import { Id } from "../_generated/dataModel";
import { 
  Prediction, 
  CreatePredictionInput, 
  UpdatePredictionInput,
  PredictionFilter,
  PredictionQueryOptions
} from "./types";
import { 
  validateCreatePredictionInput, 
  validateUpdatePredictionInput,
  validatePredictionExists,
  validateScientificModelExists,
  validateScientificModelAccess,
  validatePredictionTarget
} from "./validation";
import { 
  expandPredictionWithModel,
  expandPredictionWithTarget,
  calculatePredictionAccuracy,
  groupPredictionsByParameter
} from "./helpers";
import { getCurrentUser } from "../auth/users";

/**
 * Create a new prediction
 */
export const createPrediction = mutation({
  args: {
    prediction: v.object({
      modelId: v.id("scientificModels"),
      moleculeId: v.optional(v.id("molecules")),
      mixtureId: v.optional(v.id("mixtures")),
      parameterName: v.string(),
      value: v.union(v.string(), v.number(), v.boolean(), v.null()),
      numericValue: v.optional(v.number()),
      units: v.optional(v.string()),
      confidence: v.optional(v.number())
    })
  },
  handler: async (ctx, args) => {
    // Validate input
    validateCreatePredictionInput(args.prediction);
    validateScientificModelExists(args.prediction.modelId);
    validatePredictionTarget(args.prediction.moleculeId, args.prediction.mixtureId);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Get the model to check access
    const model = await ctx.db.get(args.prediction.modelId);
    if (!model) {
      throw new Error(`Scientific model with ID ${args.prediction.modelId} not found`);
    }
    
    // Check access to the model - only creator can add predictions
    if (!userId || !model.createdBy || !userId.equals(model.createdBy)) {
      throw new Error("Only the model creator can add predictions");
    }
    
    // Verify molecule exists if moleculeId is provided
    if (args.prediction.moleculeId) {
      const molecule = await ctx.db.get(args.prediction.moleculeId);
      if (!molecule) {
        throw new Error(`Molecule with ID ${args.prediction.moleculeId} not found`);
      }
    }
    
    // Verify mixture exists if mixtureId is provided
    if (args.prediction.mixtureId) {
      const mixture = await ctx.db.get(args.prediction.mixtureId);
      if (!mixture) {
        throw new Error(`Mixture with ID ${args.prediction.mixtureId} not found`);
      }
    }
    
    // Ensure numeric values are properly stored
    let numericValue = args.prediction.numericValue;
    
    // If value is a number and numericValue not provided, use value as numericValue
    if (typeof args.prediction.value === "number" && numericValue === undefined) {
      numericValue = args.prediction.value;
    }
    
    // Prepare data for insertion
    const now = Date.now();
    const predictionData = {
      ...args.prediction,
      numericValue,
      createdAt: now,
      updatedAt: now
    };
    
    // Insert prediction
    const predictionId = await ctx.db.insert("predictions", predictionData);
    
    // Create audit log entry
    await ctx.db.insert("scientificDataAudit", {
      table: "predictions",
      documentId: predictionId,
      operation: "create",
      userId,
      newValue: predictionData,
      timestamp: now
    });
    
    return predictionId;
  }
});

/**
 * Get a prediction by ID
 */
export const getPrediction = query({
  args: {
    predictionId: v.id("predictions"),
    includeModel: v.optional(v.boolean()),
    includeTarget: v.optional(v.boolean())
  },
  handler: async (ctx, args) => {
    // Validate prediction exists
    validatePredictionExists(args.predictionId);
    
    // Get the prediction
    const prediction = await ctx.db.get(args.predictionId);
    if (!prediction) {
      return null;
    }
    
    // Get the model to check access
    const model = await ctx.db.get(prediction.modelId);
    if (!model) {
      return null;
    }
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Check access to the model
    validateScientificModelAccess(
      prediction.modelId, 
      userId, 
      model.public, 
      model.createdBy
    );
    
    // Include target (molecule or mixture) data if requested
    if (args.includeTarget || args.includeModel) {
      return expandPredictionWithTarget(ctx, prediction);
    }
    
    // Include model data if requested
    if (args.includeModel) {
      return expandPredictionWithModel(ctx, prediction);
    }
    
    return prediction;
  }
});

/**
 * Update an existing prediction
 */
export const updatePrediction = mutation({
  args: {
    predictionId: v.id("predictions"),
    update: v.object({
      parameterName: v.optional(v.string()),
      value: v.optional(v.union(v.string(), v.number(), v.boolean(), v.null())),
      numericValue: v.optional(v.number()),
      units: v.optional(v.string()),
      confidence: v.optional(v.number())
    })
  },
  handler: async (ctx, args) => {
    // Validate input
    validateUpdatePredictionInput(args.update);
    validatePredictionExists(args.predictionId);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Get existing prediction
    const existingPrediction = await ctx.db.get(args.predictionId);
    if (!existingPrediction) {
      throw new Error(`Prediction with ID ${args.predictionId} not found`);
    }
    
    // Get the model to check access
    const model = await ctx.db.get(existingPrediction.modelId);
    if (!model) {
      throw new Error(`Scientific model with ID ${existingPrediction.modelId} not found`);
    }
    
    // Check access - only model creator can update predictions
    if (!userId || !model.createdBy || !userId.equals(model.createdBy)) {
      throw new Error("Only the model creator can update predictions");
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
    
    // Update prediction
    await ctx.db.patch(args.predictionId, updateData);
    
    // Create audit log entry
    await ctx.db.insert("scientificDataAudit", {
      table: "predictions",
      documentId: args.predictionId,
      operation: "update",
      userId,
      previousValue: existingPrediction,
      newValue: { ...existingPrediction, ...updateData },
      timestamp: now
    });
    
    return args.predictionId;
  }
});

/**
 * Delete a prediction
 */
export const deletePrediction = mutation({
  args: {
    predictionId: v.id("predictions")
  },
  handler: async (ctx, args) => {
    // Validate prediction exists
    validatePredictionExists(args.predictionId);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Get existing prediction
    const existingPrediction = await ctx.db.get(args.predictionId);
    if (!existingPrediction) {
      throw new Error(`Prediction with ID ${args.predictionId} not found`);
    }
    
    // Get the model to check access
    const model = await ctx.db.get(existingPrediction.modelId);
    if (!model) {
      throw new Error(`Scientific model with ID ${existingPrediction.modelId} not found`);
    }
    
    // Check access - only model creator can delete predictions
    if (!userId || !model.createdBy || !userId.equals(model.createdBy)) {
      throw new Error("Only the model creator can delete predictions");
    }
    
    // Create audit log entry
    const now = Date.now();
    await ctx.db.insert("scientificDataAudit", {
      table: "predictions",
      documentId: args.predictionId,
      operation: "delete",
      userId,
      previousValue: existingPrediction,
      timestamp: now
    });
    
    // Delete the prediction
    await ctx.db.delete(args.predictionId);
    
    return true;
  }
});

/**
 * List predictions with optional filtering
 */
export const listPredictions = query({
  args: {
    filter: v.optional(v.object({
      modelId: v.optional(v.id("scientificModels")),
      moleculeId: v.optional(v.id("molecules")),
      mixtureId: v.optional(v.id("mixtures")),
      parameterName: v.optional(v.string()),
      confidenceThreshold: v.optional(v.number()),
      includeModel: v.optional(v.boolean()),
      includeTarget: v.optional(v.boolean())
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
    let query = ctx.db.query("predictions");
    
    // Apply filters
    if (args.filter) {
      if (args.filter.modelId) {
        query = query.withIndex("by_model", q => 
          q.eq("modelId", args.filter!.modelId!)
        );
        
        // Get model to check access
        const model = await ctx.db.get(args.filter.modelId);
        if (!model) {
          throw new Error(`Scientific model with ID ${args.filter.modelId} not found`);
        }
        
        // Verify access to the model
        validateScientificModelAccess(
          args.filter.modelId, 
          userId, 
          model.public, 
          model.createdBy
        );
      } else {
        // If no specific model is filtered, only show predictions from models user can access
        // This is complex to implement without joins, so we'll handle this in post-processing
        // For a real implementation, we would use a more sophisticated approach
      }
      
      if (args.filter.moleculeId) {
        query = query.withIndex("by_molecule", q => 
          q.eq("moleculeId", args.filter!.moleculeId!)
        );
      }
      
      if (args.filter.mixtureId) {
        query = query.withIndex("by_mixture", q => 
          q.eq("mixtureId", args.filter!.mixtureId!)
        );
      }
      
      if (args.filter.parameterName) {
        query = query.filter(q => 
          q.eq(q.field("parameterName"), args.filter!.parameterName!)
        );
      }
      
      if (args.filter.confidenceThreshold !== undefined) {
        query = query.filter(q => 
          q.gte(q.field("confidence"), args.filter!.confidenceThreshold!)
        );
      }
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
    let predictions = await query.collect();
    
    // Post-process to filter for access rights
    // For each prediction, check if user has access to its model
    if (!args.filter?.modelId) {
      const accessiblePredictions = [];
      
      for (const prediction of predictions) {
        const model = await ctx.db.get(prediction.modelId);
        if (model) {
          try {
            // This will throw if access is denied
            validateScientificModelAccess(
              prediction.modelId, 
              userId, 
              model.public, 
              model.createdBy
            );
            
            accessiblePredictions.push(prediction);
          } catch (error) {
            // Skip predictions user doesn't have access to
            continue;
          }
        }
      }
      
      predictions = accessiblePredictions;
    }
    
    // Include additional data if requested
    if (args.filter?.includeTarget) {
      const predictionsWithTargets = [];
      
      for (const prediction of predictions) {
        predictionsWithTargets.push(await expandPredictionWithTarget(ctx, prediction));
      }
      
      return predictionsWithTargets;
    } else if (args.filter?.includeModel) {
      const predictionsWithModels = [];
      
      for (const prediction of predictions) {
        predictionsWithModels.push(await expandPredictionWithModel(ctx, prediction));
      }
      
      return predictionsWithModels;
    }
    
    return predictions;
  }
});

/**
 * Get predictions for a molecule
 */
export const getMoleculePredictions = query({
  args: {
    moleculeId: v.id("molecules"),
    parameterName: v.optional(v.string()),
    modelId: v.optional(v.id("scientificModels")),
    includeModel: v.optional(v.boolean())
  },
  handler: async (ctx, args) => {
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Set up the query
    let query = ctx.db.query("predictions")
      .withIndex("by_molecule", q => q.eq("moleculeId", args.moleculeId));
    
    // Apply additional filters
    if (args.parameterName) {
      query = query.filter(q => 
        q.eq(q.field("parameterName"), args.parameterName)
      );
    }
    
    if (args.modelId) {
      query = query.filter(q => 
        q.eq(q.field("modelId"), args.modelId)
      );
      
      // Get model to check access
      const model = await ctx.db.get(args.modelId);
      if (!model) {
        throw new Error(`Scientific model with ID ${args.modelId} not found`);
      }
      
      // Verify access to the model
      validateScientificModelAccess(
        args.modelId, 
        userId, 
        model.public, 
        model.createdBy
      );
    }
    
    // Execute query
    let predictions = await query.collect();
    
    // Post-process to filter for access rights if no specific model provided
    if (!args.modelId) {
      const accessiblePredictions = [];
      
      for (const prediction of predictions) {
        const model = await ctx.db.get(prediction.modelId);
        if (model) {
          try {
            // This will throw if access is denied
            validateScientificModelAccess(
              prediction.modelId, 
              userId, 
              model.public, 
              model.createdBy
            );
            
            accessiblePredictions.push(prediction);
          } catch (error) {
            // Skip predictions user doesn't have access to
            continue;
          }
        }
      }
      
      predictions = accessiblePredictions;
    }
    
    // Include model data if requested
    if (args.includeModel) {
      const predictionsWithModels = [];
      
      for (const prediction of predictions) {
        predictionsWithModels.push(await expandPredictionWithModel(ctx, prediction));
      }
      
      return predictionsWithModels;
    }
    
    return predictions;
  }
});

/**
 * Get predictions for a mixture
 */
export const getMixturePredictions = query({
  args: {
    mixtureId: v.id("mixtures"),
    parameterName: v.optional(v.string()),
    modelId: v.optional(v.id("scientificModels")),
    includeModel: v.optional(v.boolean())
  },
  handler: async (ctx, args) => {
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Set up the query
    let query = ctx.db.query("predictions")
      .withIndex("by_mixture", q => q.eq("mixtureId", args.mixtureId));
    
    // Apply additional filters
    if (args.parameterName) {
      query = query.filter(q => 
        q.eq(q.field("parameterName"), args.parameterName)
      );
    }
    
    if (args.modelId) {
      query = query.filter(q => 
        q.eq(q.field("modelId"), args.modelId)
      );
      
      // Get model to check access
      const model = await ctx.db.get(args.modelId);
      if (!model) {
        throw new Error(`Scientific model with ID ${args.modelId} not found`);
      }
      
      // Verify access to the model
      validateScientificModelAccess(
        args.modelId, 
        userId, 
        model.public, 
        model.createdBy
      );
    }
    
    // Execute query
    let predictions = await query.collect();
    
    // Post-process to filter for access rights if no specific model provided
    if (!args.modelId) {
      const accessiblePredictions = [];
      
      for (const prediction of predictions) {
        const model = await ctx.db.get(prediction.modelId);
        if (model) {
          try {
            // This will throw if access is denied
            validateScientificModelAccess(
              prediction.modelId, 
              userId, 
              model.public, 
              model.createdBy
            );
            
            accessiblePredictions.push(prediction);
          } catch (error) {
            // Skip predictions user doesn't have access to
            continue;
          }
        }
      }
      
      predictions = accessiblePredictions;
    }
    
    // Include model data if requested
    if (args.includeModel) {
      const predictionsWithModels = [];
      
      for (const prediction of predictions) {
        predictionsWithModels.push(await expandPredictionWithModel(ctx, prediction));
      }
      
      return predictionsWithModels;
    }
    
    return predictions;
  }
});

/**
 * Batch add multiple predictions
 */
export const batchAddPredictions = mutation({
  args: {
    modelId: v.id("scientificModels"),
    predictions: v.array(
      v.object({
        moleculeId: v.optional(v.id("molecules")),
        mixtureId: v.optional(v.id("mixtures")),
        parameterName: v.string(),
        value: v.union(v.string(), v.number(), v.boolean(), v.null()),
        numericValue: v.optional(v.number()),
        units: v.optional(v.string()),
        confidence: v.optional(v.number())
      })
    )
  },
  handler: async (ctx, args) => {
    // Validate model exists
    validateScientificModelExists(args.modelId);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Get the model to check access
    const model = await ctx.db.get(args.modelId);
    if (!model) {
      throw new Error(`Scientific model with ID ${args.modelId} not found`);
    }
    
    // Check access - only model creator can add predictions
    if (!userId || !model.createdBy || !userId.equals(model.createdBy)) {
      throw new Error("Only the model creator can add predictions");
    }
    
    // Record time for all operations
    const now = Date.now();
    
    // Validate and insert each prediction
    const predictionIds: Id<"predictions">[] = [];
    
    for (const prediction of args.predictions) {
      // Validate prediction input
      validatePredictionTarget(prediction.moleculeId, prediction.mixtureId);
      validateCreatePredictionInput({
        ...prediction,
        modelId: args.modelId
      });
      
      // Verify target exists
      if (prediction.moleculeId) {
        const molecule = await ctx.db.get(prediction.moleculeId);
        if (!molecule) {
          throw new Error(`Molecule with ID ${prediction.moleculeId} not found`);
        }
      }
      
      if (prediction.mixtureId) {
        const mixture = await ctx.db.get(prediction.mixtureId);
        if (!mixture) {
          throw new Error(`Mixture with ID ${prediction.mixtureId} not found`);
        }
      }
      
      // Ensure numeric values are properly stored
      let numericValue = prediction.numericValue;
      
      // If value is a number and numericValue not provided, use value as numericValue
      if (typeof prediction.value === "number" && numericValue === undefined) {
        numericValue = prediction.value;
      }
      
      // Prepare prediction data
      const predictionData = {
        ...prediction,
        modelId: args.modelId,
        numericValue,
        createdAt: now,
        updatedAt: now
      };
      
      // Insert prediction
      const predictionId = await ctx.db.insert("predictions", predictionData);
      
      // Create audit log entry
      await ctx.db.insert("scientificDataAudit", {
        table: "predictions",
        documentId: predictionId,
        operation: "create",
        userId,
        newValue: predictionData,
        timestamp: now
      });
      
      predictionIds.push(predictionId);
    }
    
    return predictionIds;
  }
});
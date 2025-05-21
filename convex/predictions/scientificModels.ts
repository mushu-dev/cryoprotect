/**
 * CRUD operations for scientific models
 */

import { mutation, query } from "../_generated/server";
import { v } from "convex/values";
import { Id } from "../_generated/dataModel";
import { 
  ScientificModel, 
  CreateScientificModelInput, 
  UpdateScientificModelInput,
  ScientificModelFilter,
  ScientificModelQueryOptions
} from "./types";
import { 
  validateCreateScientificModelInput, 
  validateUpdateScientificModelInput,
  validateScientificModelExists,
  validateScientificModelAccess
} from "./validation";
import { 
  formatModelType,
  modelHasPredictions,
  parseModelParameters,
  stringifyModelParameters
} from "./helpers";
import { getCurrentUser } from "../auth/users";

/**
 * Create a new scientific model
 */
export const createScientificModel = mutation({
  args: {
    model: v.object({
      name: v.string(),
      description: v.optional(v.string()),
      version: v.string(),
      type: v.string(),
      parameters: v.optional(v.string()),
      public: v.optional(v.boolean())
    })
  },
  handler: async (ctx, args) => {
    // Validate input
    validateCreateScientificModelInput(args.model);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // User must be logged in to create a model
    if (!userId) {
      throw new Error("Authentication required to create a scientific model");
    }
    
    // Prepare data for insertion
    const now = Date.now();
    const modelData = {
      ...args.model,
      createdBy: userId,
      createdAt: now,
      updatedAt: now,
      public: args.model.public ?? false
    };
    
    // Insert model
    const modelId = await ctx.db.insert("scientificModels", modelData);
    
    // Create audit log entry
    await ctx.db.insert("scientificDataAudit", {
      table: "scientificModels",
      documentId: modelId,
      operation: "create",
      userId,
      newValue: modelData,
      timestamp: now,
    });
    
    return modelId;
  }
});

/**
 * Get a scientific model by ID
 */
export const getScientificModel = query({
  args: {
    modelId: v.id("scientificModels")
  },
  handler: async (ctx, args) => {
    // Validate model exists
    validateScientificModelExists(args.modelId);
    
    // Get the model
    const model = await ctx.db.get(args.modelId);
    if (!model) {
      return null;
    }
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Check access
    validateScientificModelAccess(
      args.modelId, 
      userId, 
      model.public, 
      model.createdBy
    );
    
    return model;
  }
});

/**
 * Update an existing scientific model
 */
export const updateScientificModel = mutation({
  args: {
    modelId: v.id("scientificModels"),
    update: v.object({
      name: v.optional(v.string()),
      description: v.optional(v.string()),
      version: v.optional(v.string()),
      type: v.optional(v.string()),
      parameters: v.optional(v.string()),
      public: v.optional(v.boolean())
    })
  },
  handler: async (ctx, args) => {
    // Validate input
    validateUpdateScientificModelInput(args.update);
    validateScientificModelExists(args.modelId);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Get existing model
    const existingModel = await ctx.db.get(args.modelId);
    if (!existingModel) {
      throw new Error(`Scientific model with ID ${args.modelId} not found`);
    }
    
    // Check access - only creator can update
    if (!userId || !existingModel.createdBy || !userId.equals(existingModel.createdBy)) {
      throw new Error("Only the creator can update a scientific model");
    }
    
    // Prepare update data
    const now = Date.now();
    const updateData = {
      ...args.update,
      updatedAt: now
    };
    
    // Update model
    await ctx.db.patch(args.modelId, updateData);
    
    // Create audit log entry
    await ctx.db.insert("scientificDataAudit", {
      table: "scientificModels",
      documentId: args.modelId,
      operation: "update",
      userId,
      previousValue: existingModel,
      newValue: { ...existingModel, ...updateData },
      timestamp: now
    });
    
    return args.modelId;
  }
});

/**
 * Delete a scientific model
 */
export const deleteScientificModel = mutation({
  args: {
    modelId: v.id("scientificModels")
  },
  handler: async (ctx, args) => {
    // Validate model exists
    validateScientificModelExists(args.modelId);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Get existing model
    const existingModel = await ctx.db.get(args.modelId);
    if (!existingModel) {
      throw new Error(`Scientific model with ID ${args.modelId} not found`);
    }
    
    // Check access - only creator can delete
    if (!userId || !existingModel.createdBy || !userId.equals(existingModel.createdBy)) {
      throw new Error("Only the creator can delete a scientific model");
    }
    
    // Check if the model has predictions
    const hasPredictions = await modelHasPredictions(ctx, args.modelId);
    if (hasPredictions) {
      throw new Error("Cannot delete a model that has predictions. Delete the predictions first.");
    }
    
    // Create audit log entry for model deletion
    const now = Date.now();
    await ctx.db.insert("scientificDataAudit", {
      table: "scientificModels",
      documentId: args.modelId,
      operation: "delete",
      userId,
      previousValue: existingModel,
      timestamp: now
    });
    
    // Delete the model
    await ctx.db.delete(args.modelId);
    
    return true;
  }
});

/**
 * List scientific models with optional filtering
 */
export const listScientificModels = query({
  args: {
    filter: v.optional(v.object({
      name: v.optional(v.string()),
      type: v.optional(v.string()),
      version: v.optional(v.string()),
      createdBy: v.optional(v.id("users")),
      public: v.optional(v.boolean())
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
    let query = ctx.db.query("scientificModels");
    
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
      
      if (args.filter.createdBy) {
        query = query.withIndex("by_creator", q => 
          q.eq("createdBy", args.filter!.createdBy!)
        );
      }
      
      if (args.filter.version) {
        // We'd typically use the by_version index, but need to include name for it
        // For simplicity, using filter here
        query = query.filter(q => 
          q.eq(q.field("version"), args.filter!.version!)
        );
      }
    }
    
    // Filter based on access permissions
    if (!userId) {
      // If no user is logged in, only show public models
      query = query.withIndex("by_public", q => q.eq("public", true));
    } else {
      // If a user is logged in, show public models and their own
      query = query.filter(q => 
        q.or(
          q.eq(q.field("public"), true),
          q.eq(q.field("createdBy"), userId)
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
        case "type":
          query = query.order("type", sortDirection);
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
    return query.collect();
  }
});

/**
 * Search scientific models by name
 */
export const searchScientificModels = query({
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
    let query = ctx.db.query("scientificModels")
      .filter(q => 
        q.includes("name", args.query)
      );
    
    // Filter based on access permissions
    if (!userId) {
      // If no user is logged in, only show public models
      query = query.filter(q => q.eq(q.field("public"), true));
    } else {
      // If a user is logged in, show public models and their own
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

/**
 * Update scientific model parameters
 */
export const updateModelParameters = mutation({
  args: {
    modelId: v.id("scientificModels"),
    parameters: v.any()
  },
  handler: async (ctx, args) => {
    // Validate model exists
    validateScientificModelExists(args.modelId);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Get existing model
    const existingModel = await ctx.db.get(args.modelId);
    if (!existingModel) {
      throw new Error(`Scientific model with ID ${args.modelId} not found`);
    }
    
    // Check access - only creator can update
    if (!userId || !existingModel.createdBy || !userId.equals(existingModel.createdBy)) {
      throw new Error("Only the creator can update model parameters");
    }
    
    // Convert parameters to string
    const parametersString = stringifyModelParameters(args.parameters);
    
    // Prepare update data
    const now = Date.now();
    const updateData = {
      parameters: parametersString,
      updatedAt: now
    };
    
    // Update model
    await ctx.db.patch(args.modelId, updateData);
    
    // Create audit log entry
    await ctx.db.insert("scientificDataAudit", {
      table: "scientificModels",
      documentId: args.modelId,
      operation: "update",
      userId,
      previousValue: existingModel,
      newValue: { ...existingModel, ...updateData },
      timestamp: now
    });
    
    return args.modelId;
  }
});

/**
 * Get scientific model parameters as an object
 */
export const getModelParameters = query({
  args: {
    modelId: v.id("scientificModels")
  },
  handler: async (ctx, args) => {
    // Validate model exists
    validateScientificModelExists(args.modelId);
    
    // Get the model
    const model = await ctx.db.get(args.modelId);
    if (!model) {
      return null;
    }
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Check access
    validateScientificModelAccess(
      args.modelId, 
      userId, 
      model.public, 
      model.createdBy
    );
    
    // Parse parameters
    return parseModelParameters(model.parameters);
  }
});
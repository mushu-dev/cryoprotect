/**
 * CRUD operations for data sources
 */

import { mutation, query } from "../_generated/server";
import { v } from "convex/values";
import { Id } from "../_generated/dataModel";
import { 
  DataSource,
  CreateDataSourceInput, 
  UpdateDataSourceInput,
  DataSourceFilter,
  DataSourceQueryOptions
} from "./types";
import { 
  validateCreateDataSourceInput, 
  validateUpdateDataSourceInput,
  validateDataSourceExists,
  validateDataSourceNotReferenced
} from "./validation";
import { 
  formatDataSourceType,
  getDataSourceUsageStats,
  expandDataSourceWithUsage,
  isDataSourceInUse
} from "./helpers";
import { getCurrentUser } from "../auth/users";

/**
 * Create a new data source
 */
export const createDataSource = mutation({
  args: {
    dataSource: v.object({
      name: v.string(),
      description: v.optional(v.string()),
      url: v.optional(v.string()),
      type: v.string(),
      version: v.optional(v.string())
    })
  },
  handler: async (ctx, args) => {
    // Validate input
    validateCreateDataSourceInput(args.dataSource);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // User must be logged in to create a data source
    if (!userId) {
      throw new Error("Authentication required to create a data source");
    }
    
    // Prepare data for insertion
    const now = Date.now();
    const dataSourceData = {
      ...args.dataSource,
      createdAt: now,
      updatedAt: now
    };
    
    // Insert data source
    const dataSourceId = await ctx.db.insert("dataSources", dataSourceData);
    
    // Create audit log entry
    await ctx.db.insert("scientificDataAudit", {
      table: "dataSources",
      documentId: dataSourceId,
      operation: "create",
      userId,
      newValue: dataSourceData,
      timestamp: now,
    });
    
    return dataSourceId;
  }
});

/**
 * Get a data source by ID
 */
export const getDataSource = query({
  args: {
    dataSourceId: v.id("dataSources"),
    includeUsageStats: v.optional(v.boolean())
  },
  handler: async (ctx, args) => {
    // Validate data source exists
    validateDataSourceExists(args.dataSourceId);
    
    // Get the data source
    const dataSource = await ctx.db.get(args.dataSourceId);
    if (!dataSource) {
      return null;
    }
    
    // Include usage statistics if requested
    if (args.includeUsageStats) {
      return expandDataSourceWithUsage(ctx, dataSource);
    }
    
    return dataSource;
  }
});

/**
 * Update an existing data source
 */
export const updateDataSource = mutation({
  args: {
    dataSourceId: v.id("dataSources"),
    update: v.object({
      name: v.optional(v.string()),
      description: v.optional(v.string()),
      url: v.optional(v.string()),
      type: v.optional(v.string()),
      version: v.optional(v.string())
    })
  },
  handler: async (ctx, args) => {
    // Validate input
    validateUpdateDataSourceInput(args.update);
    validateDataSourceExists(args.dataSourceId);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // User must be logged in to update a data source
    if (!userId) {
      throw new Error("Authentication required to update a data source");
    }
    
    // Get existing data source
    const existingDataSource = await ctx.db.get(args.dataSourceId);
    if (!existingDataSource) {
      throw new Error(`Data source with ID ${args.dataSourceId} not found`);
    }
    
    // Prepare update data
    const now = Date.now();
    const updateData = {
      ...args.update,
      updatedAt: now
    };
    
    // Update data source
    await ctx.db.patch(args.dataSourceId, updateData);
    
    // Create audit log entry
    await ctx.db.insert("scientificDataAudit", {
      table: "dataSources",
      documentId: args.dataSourceId,
      operation: "update",
      userId,
      previousValue: existingDataSource,
      newValue: { ...existingDataSource, ...updateData },
      timestamp: now
    });
    
    return args.dataSourceId;
  }
});

/**
 * Delete a data source
 */
export const deleteDataSource = mutation({
  args: {
    dataSourceId: v.id("dataSources")
  },
  handler: async (ctx, args) => {
    // Validate data source exists
    validateDataSourceExists(args.dataSourceId);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // User must be logged in to delete a data source
    if (!userId) {
      throw new Error("Authentication required to delete a data source");
    }
    
    // Get existing data source
    const existingDataSource = await ctx.db.get(args.dataSourceId);
    if (!existingDataSource) {
      throw new Error(`Data source with ID ${args.dataSourceId} not found`);
    }
    
    // Check if the data source is in use
    const usageStats = await getDataSourceUsageStats(ctx, args.dataSourceId);
    validateDataSourceNotReferenced(
      args.dataSourceId,
      usageStats.moleculeCount,
      usageStats.propertyCount,
      usageStats.experimentCount
    );
    
    // Create audit log entry
    const now = Date.now();
    await ctx.db.insert("scientificDataAudit", {
      table: "dataSources",
      documentId: args.dataSourceId,
      operation: "delete",
      userId,
      previousValue: existingDataSource,
      timestamp: now
    });
    
    // Delete the data source
    await ctx.db.delete(args.dataSourceId);
    
    return true;
  }
});

/**
 * List data sources with optional filtering
 */
export const listDataSources = query({
  args: {
    filter: v.optional(v.object({
      name: v.optional(v.string()),
      type: v.optional(v.string()),
      includeUsageStats: v.optional(v.boolean())
    })),
    options: v.optional(v.object({
      limit: v.optional(v.number()),
      cursor: v.optional(v.string()),
      sortBy: v.optional(v.string()),
      sortDirection: v.optional(v.string())
    }))
  },
  handler: async (ctx, args) => {
    // Set up the query
    let query = ctx.db.query("dataSources");
    
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
          query = query.order("name", "asc"); // default sort
      }
    } else {
      // Default sort by name ascending
      query = query.order("name", "asc");
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
    const dataSources = await query.collect();
    
    // Include usage statistics if requested
    if (args.filter?.includeUsageStats) {
      const dataSourcesWithUsage = [];
      
      for (const dataSource of dataSources) {
        dataSourcesWithUsage.push(await expandDataSourceWithUsage(ctx, dataSource));
      }
      
      return dataSourcesWithUsage;
    }
    
    return dataSources;
  }
});

/**
 * Search data sources by name
 */
export const searchDataSources = query({
  args: {
    query: v.string(),
    type: v.optional(v.string()),
    limit: v.optional(v.number())
  },
  handler: async (ctx, args) => {
    // Set up the query
    let query = ctx.db.query("dataSources")
      .filter(q => 
        q.includes("name", args.query)
      );
    
    // Filter by type if provided
    if (args.type) {
      query = query.filter(q => 
        q.eq(q.field("type"), args.type)
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
 * Check if a data source can be safely deleted
 */
export const canDeleteDataSource = query({
  args: {
    dataSourceId: v.id("dataSources")
  },
  handler: async (ctx, args) => {
    // Validate data source exists
    validateDataSourceExists(args.dataSourceId);
    
    // Check if the data source is in use
    const inUse = await isDataSourceInUse(ctx, args.dataSourceId);
    
    return !inUse;
  }
});
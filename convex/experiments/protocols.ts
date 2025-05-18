/**
 * CRUD operations for experiment protocols
 */

import { mutation, query } from "../_generated/server";
import { v } from "convex/values";
import { Id } from "../_generated/dataModel";
import { 
  Protocol,
  CreateProtocolInput,
  UpdateProtocolInput
} from "./enhanced_types";
import { 
  validateCreateProtocolInput,
  validateUpdateProtocolInput,
  validateProtocolExists
} from "./enhanced_validation";
import { getCurrentUser } from "../auth/users";

/**
 * Create a new protocol
 */
export const createProtocol = mutation({
  args: {
    protocol: v.object({
      name: v.string(),
      description: v.optional(v.string()),
      steps: v.array(v.object({
        name: v.string(),
        description: v.optional(v.string()),
        parameters: v.optional(v.map(v.string(), v.any())),
        duration: v.optional(v.number()),
        durationUnit: v.optional(v.string()),
        temperature: v.optional(v.number()),
        temperatureUnit: v.optional(v.string())
      })),
      parameters: v.optional(v.map(v.string(), v.any())),
      version: v.optional(v.number()),
      parentId: v.optional(v.id("protocols")),
      isTemplate: v.optional(v.boolean()),
      category: v.optional(v.string()),
      public: v.optional(v.boolean())
    })
  },
  handler: async (ctx, args) => {
    // Validate input
    validateCreateProtocolInput(args.protocol);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // If parentId is provided, verify it exists
    if (args.protocol.parentId) {
      await validateProtocolExists(ctx.db, args.protocol.parentId);
    }
    
    // Prepare data for insertion
    const now = Date.now();
    const protocolData = {
      ...args.protocol,
      version: args.protocol.version || 1,
      isTemplate: args.protocol.isTemplate ?? false,
      createdBy: userId,
      createdAt: now,
      updatedAt: now,
      public: args.protocol.public ?? false
    };
    
    // Insert protocol
    const protocolId = await ctx.db.insert("protocols", protocolData);
    
    // Create audit log entry
    await ctx.db.insert("scientificDataAudit", {
      table: "protocols",
      documentId: protocolId,
      operation: "create",
      userId,
      newValue: protocolData,
      timestamp: now,
    });
    
    return protocolId;
  }
});

/**
 * Get a protocol by ID
 */
export const getProtocol = query({
  args: {
    protocolId: v.id("protocols")
  },
  handler: async (ctx, args) => {
    // Get the protocol
    const protocol = await ctx.db.get(args.protocolId);
    if (!protocol) {
      return null;
    }
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Check access
    if (!protocol.public) {
      if (!userId) {
        throw new Error("Authentication required to access this protocol");
      }
      
      if (!protocol.createdBy || !userId.equals(protocol.createdBy)) {
        throw new Error("You do not have permission to access this protocol");
      }
    }
    
    return protocol;
  }
});

/**
 * Update an existing protocol
 */
export const updateProtocol = mutation({
  args: {
    protocolId: v.id("protocols"),
    update: v.object({
      name: v.optional(v.string()),
      description: v.optional(v.string()),
      steps: v.optional(v.array(v.object({
        name: v.string(),
        description: v.optional(v.string()),
        parameters: v.optional(v.map(v.string(), v.any())),
        duration: v.optional(v.number()),
        durationUnit: v.optional(v.string()),
        temperature: v.optional(v.number()),
        temperatureUnit: v.optional(v.string())
      }))),
      parameters: v.optional(v.map(v.string(), v.any())),
      version: v.optional(v.number()),
      category: v.optional(v.string()),
      public: v.optional(v.boolean())
    })
  },
  handler: async (ctx, args) => {
    // Validate input
    validateUpdateProtocolInput(args.update);
    await validateProtocolExists(ctx.db, args.protocolId);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Get existing protocol
    const existingProtocol = await ctx.db.get(args.protocolId);
    if (!existingProtocol) {
      throw new Error(`Protocol with ID ${args.protocolId} not found`);
    }
    
    // Check access - only creator can update
    if (!userId || 
        !existingProtocol.createdBy || 
        !userId.equals(existingProtocol.createdBy)) {
      throw new Error("Only the creator can update a protocol");
    }
    
    // Prepare update data
    const now = Date.now();
    const updateData = {
      ...args.update,
      updatedAt: now
    };
    
    // Update protocol
    await ctx.db.patch(args.protocolId, updateData);
    
    // Create audit log entry
    await ctx.db.insert("scientificDataAudit", {
      table: "protocols",
      documentId: args.protocolId,
      operation: "update",
      userId,
      previousValue: existingProtocol,
      newValue: { ...existingProtocol, ...updateData },
      timestamp: now
    });
    
    return args.protocolId;
  }
});

/**
 * Create a new protocol version
 */
export const createProtocolVersion = mutation({
  args: {
    protocol: v.object({
      name: v.string(),
      description: v.optional(v.string()),
      steps: v.array(v.object({
        name: v.string(),
        description: v.optional(v.string()),
        parameters: v.optional(v.map(v.string(), v.any())),
        duration: v.optional(v.number()),
        durationUnit: v.optional(v.string()),
        temperature: v.optional(v.number()),
        temperatureUnit: v.optional(v.string())
      })),
      parameters: v.optional(v.map(v.string(), v.any())),
      parentId: v.id("protocols"),
      category: v.optional(v.string()),
      public: v.optional(v.boolean())
    })
  },
  handler: async (ctx, args) => {
    // Validate input
    validateCreateProtocolInput(args.protocol);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Verify parent exists
    await validateProtocolExists(ctx.db, args.protocol.parentId);
    
    // Get parent protocol to calculate new version number
    const parentProtocol = await ctx.db.get(args.protocol.parentId);
    if (!parentProtocol) {
      throw new Error(`Parent protocol with ID ${args.protocol.parentId} not found`);
    }
    
    // Ensure user has access to parent
    if (!parentProtocol.public) {
      if (!userId) {
        throw new Error("Authentication required to create a new version");
      }
      
      if (!parentProtocol.createdBy || !userId.equals(parentProtocol.createdBy)) {
        throw new Error("You do not have permission to create a new version");
      }
    }
    
    // Get the highest version number among protocols with the same name
    const highestVersionProtocol = await ctx.db
      .query("protocols")
      .filter(q => q.eq(q.field("name"), args.protocol.name))
      .order("version", "desc")
      .first();
    
    const newVersion = highestVersionProtocol 
      ? (highestVersionProtocol.version || 1) + 1 
      : 1;
    
    // Prepare data for insertion
    const now = Date.now();
    const protocolData = {
      ...args.protocol,
      version: newVersion,
      isTemplate: parentProtocol.isTemplate,
      createdBy: userId,
      createdAt: now,
      updatedAt: now,
      public: args.protocol.public ?? parentProtocol.public
    };
    
    // Insert protocol
    const protocolId = await ctx.db.insert("protocols", protocolData);
    
    // Create audit log entry
    await ctx.db.insert("scientificDataAudit", {
      table: "protocols",
      documentId: protocolId,
      operation: "create",
      userId,
      newValue: protocolData,
      timestamp: now,
    });
    
    return protocolId;
  }
});

/**
 * Delete a protocol
 */
export const deleteProtocol = mutation({
  args: {
    protocolId: v.id("protocols")
  },
  handler: async (ctx, args) => {
    // Validate protocol exists
    await validateProtocolExists(ctx.db, args.protocolId);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Get existing protocol
    const existingProtocol = await ctx.db.get(args.protocolId);
    if (!existingProtocol) {
      throw new Error(`Protocol with ID ${args.protocolId} not found`);
    }
    
    // Check access - only creator can delete
    if (!userId || 
        !existingProtocol.createdBy || 
        !userId.equals(existingProtocol.createdBy)) {
      throw new Error("Only the creator can delete a protocol");
    }
    
    // Check if protocol is in use by any experiments
    const experimentUsingProtocol = await ctx.db
      .query("enhancedExperiments")
      .withIndex("by_protocol", q => q.eq("protocolId", args.protocolId))
      .first();
    
    if (experimentUsingProtocol) {
      throw new Error("Cannot delete protocol that is in use by experiments");
    }
    
    // Record time for audit log
    const now = Date.now();
    
    // Create audit log entry
    await ctx.db.insert("scientificDataAudit", {
      table: "protocols",
      documentId: args.protocolId,
      operation: "delete",
      userId,
      previousValue: existingProtocol,
      timestamp: now
    });
    
    // Delete the protocol
    await ctx.db.delete(args.protocolId);
    
    return true;
  }
});

/**
 * List protocols with optional filtering
 */
export const listProtocols = query({
  args: {
    filter: v.optional(v.object({
      name: v.optional(v.string()),
      category: v.optional(v.string()),
      isTemplate: v.optional(v.boolean()),
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
    let query = ctx.db.query("protocols");
    
    // Apply filters
    if (args.filter) {
      if (args.filter.name) {
        query = query.filter(q => 
          q.includes("name", args.filter!.name!)
        );
      }
      
      if (args.filter.category) {
        query = query.withIndex("by_category", q => 
          q.eq("category", args.filter!.category!)
        );
      }
      
      if (args.filter.isTemplate !== undefined) {
        query = query.withIndex("by_template", q => 
          q.eq("isTemplate", args.filter!.isTemplate!)
        );
      }
    }
    
    // Filter based on access permissions
    if (!userId) {
      // If no user is logged in, only show public protocols
      query = query.withIndex("by_public", q => q.eq("public", true));
    } else {
      // If a user is logged in, show public protocols and their own
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
        case "category":
          query = query.order("category", sortDirection);
          break;
        case "version":
          query = query.order("version", sortDirection);
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
 * Get protocol version history
 */
export const getProtocolVersionHistory = query({
  args: {
    name: v.string()
  },
  handler: async (ctx, args) => {
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Get all versions of the protocol with this name
    let query = ctx.db
      .query("protocols")
      .filter(q => q.eq(q.field("name"), args.name))
      .order("version", "desc");
    
    // Filter based on access permissions
    if (!userId) {
      // If no user is logged in, only show public protocols
      query = query.filter(q => q.eq(q.field("public"), true));
    } else {
      // If a user is logged in, show public protocols and their own
      query = query.filter(q => 
        q.or(
          q.eq(q.field("public"), true),
          q.eq(q.field("createdBy"), userId)
        )
      );
    }
    
    // Execute query
    return query.collect();
  }
});
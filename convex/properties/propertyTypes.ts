/**
 * Functions for managing property types
 */

import { mutation, query } from "../_generated/server";
import { v } from "convex/values";
import { Id } from "../_generated/dataModel";
import { PropertyType, CreatePropertyTypeInput, UpdatePropertyTypeInput, PropertyTypeFilter, PropertyTypeQueryOptions } from "./types";
import { validateCreatePropertyTypeInput, validateUpdatePropertyTypeInput } from "./validation";
import { getPropertyTypeByName } from "./helpers";
import { createAuditLog, getCurrentTimestamp, pagination, validateAuthenticated, validateRole } from "../utils/common";

/**
 * Create a new property type
 */
export const createPropertyType = mutation({
  args: {
    input: v.object({
      name: v.string(),
      displayName: v.string(),
      description: v.optional(v.string()),
      dataType: v.union(v.literal("string"), v.literal("number"), v.literal("boolean")),
      units: v.optional(v.string()),
      defaultUnits: v.optional(v.string()),
      category: v.optional(v.string()),
      isCalculated: v.optional(v.boolean()),
    }),
  },
  handler: async (ctx, args): Promise<PropertyType> => {
    // Validate user authentication
    validateAuthenticated(ctx);
    
    // Only admins and scientists can create property types
    await validateRole(ctx, ["admin", "scientist"]);
    
    // Validate input
    const validationResult = validateCreatePropertyTypeInput(args.input);
    if (!validationResult.valid) {
      throw new Error(validationResult.error);
    }
    
    // Check if property type with this name already exists
    const existingPropertyType = await getPropertyTypeByName(ctx.db, args.input.name);
    if (existingPropertyType) {
      throw new Error(`Property type with name ${args.input.name} already exists`);
    }
    
    // Create property type document
    const now = getCurrentTimestamp();
    const propertyTypeId = await ctx.db.insert("propertyTypes", {
      name: args.input.name,
      displayName: args.input.displayName,
      description: args.input.description,
      dataType: args.input.dataType,
      units: args.input.units,
      defaultUnits: args.input.defaultUnits,
      category: args.input.category,
      isCalculated: args.input.isCalculated ?? false,
      createdAt: now,
      updatedAt: now,
    });
    
    // Fetch the inserted property type
    const propertyType = await ctx.db.get(propertyTypeId) as PropertyType;
    
    // Create audit log
    await createAuditLog(
      ctx,
      "propertyTypes",
      propertyTypeId.toString(),
      "create",
      undefined,
      propertyType,
      "Initial creation"
    );
    
    return propertyType;
  },
});

/**
 * Update an existing property type
 */
export const updatePropertyType = mutation({
  args: {
    id: v.id("propertyTypes"),
    input: v.object({
      displayName: v.optional(v.string()),
      description: v.optional(v.string()),
      units: v.optional(v.string()),
      defaultUnits: v.optional(v.string()),
      category: v.optional(v.string()),
      isCalculated: v.optional(v.boolean()),
    }),
  },
  handler: async (ctx, args): Promise<PropertyType> => {
    // Validate user authentication
    validateAuthenticated(ctx);
    
    // Only admins and scientists can update property types
    await validateRole(ctx, ["admin", "scientist"]);
    
    // Validate input
    const validationResult = validateUpdatePropertyTypeInput(args.input);
    if (!validationResult.valid) {
      throw new Error(validationResult.error);
    }
    
    // Get the current property type
    const currentPropertyType = await ctx.db.get(args.id);
    if (!currentPropertyType) {
      throw new Error("Property type not found");
    }
    
    // Create update object
    const updateObject = {
      ...args.input,
      updatedAt: getCurrentTimestamp(),
    };
    
    // Update the property type
    const updatedPropertyType = await ctx.db.patch(args.id, updateObject) as PropertyType;
    
    // Create audit log
    await createAuditLog(
      ctx,
      "propertyTypes",
      args.id.toString(),
      "update",
      currentPropertyType,
      updatedPropertyType,
      "User update"
    );
    
    return updatedPropertyType;
  },
});

/**
 * Get a property type by ID
 */
export const getPropertyType = query({
  args: {
    id: v.id("propertyTypes"),
  },
  handler: async (ctx, args): Promise<PropertyType | null> => {
    return await ctx.db.get(args.id) as PropertyType | null;
  },
});

/**
 * Search for property types based on filters
 */
export const searchPropertyTypes = query({
  args: {
    filter: v.optional(
      v.object({
        name: v.optional(v.string()),
        dataType: v.optional(v.union(
          v.literal("string"),
          v.literal("number"),
          v.literal("boolean")
        )),
        category: v.optional(v.string()),
        isCalculated: v.optional(v.boolean()),
      })
    ),
    options: v.optional(
      v.object({
        limit: v.optional(v.number()),
        cursor: v.optional(v.string()),
      })
    ),
  },
  handler: async (ctx, args): Promise<{
    propertyTypes: PropertyType[];
    cursor: string | null;
    totalCount: number;
  }> => {
    const filter = args.filter || {};
    const options = args.options || {};
    const limit = options.limit || pagination.defaultPageSize;
    
    // Start building the query
    let query = ctx.db.query("propertyTypes");
    
    // Apply filters
    if (filter.name) {
      query = query.withSearchIndex("by_name", (q) => 
        q.search("name", filter.name)
      );
    }
    
    if (filter.dataType) {
      query = query.filter(q => q.eq(q.field("dataType"), filter.dataType));
    }
    
    if (filter.category) {
      query = query.withIndex("by_category", (q) => 
        q.eq("category", filter.category)
      );
    }
    
    if (filter.isCalculated !== undefined) {
      query = query.filter(q => q.eq(q.field("isCalculated"), filter.isCalculated));
    }
    
    // Apply pagination
    const results = await query
      .paginate({ cursor: options.cursor, numItems: limit });
    
    // Get total count (approximate)
    const totalCount = await query.collect();
    
    return {
      propertyTypes: results.page as PropertyType[],
      cursor: results.continueCursor,
      totalCount: totalCount.length,
    };
  },
});

/**
 * Delete a property type (admin only)
 * 
 * Note: This operation will fail if the property type is in use
 */
export const deletePropertyType = mutation({
  args: {
    id: v.id("propertyTypes"),
    force: v.optional(v.boolean()),
  },
  handler: async (ctx, args): Promise<{ success: boolean; message: string }> => {
    // Validate user authentication and role
    await validateRole(ctx, ["admin"]);
    
    // Get the current property type
    const propertyType = await ctx.db.get(args.id);
    if (!propertyType) {
      throw new Error("Property type not found");
    }
    
    // Check if property type is in use if not force deleting
    if (!args.force) {
      const propertiesUsingType = await ctx.db
        .query("molecularProperties")
        .withIndex("by_property_type", q => q.eq("propertyTypeId", args.id))
        .first();
      
      if (propertiesUsingType) {
        throw new Error(
          "Cannot delete property type that is in use. Use force=true to delete anyway."
        );
      }
    }
    
    // If force is true, delete all dependent records
    if (args.force) {
      const properties = await ctx.db
        .query("molecularProperties")
        .withIndex("by_property_type", q => q.eq("propertyTypeId", args.id))
        .collect();
      
      for (const property of properties) {
        await ctx.db.delete(property._id);
      }
    }
    
    // Create audit log before deletion
    await createAuditLog(
      ctx,
      "propertyTypes",
      args.id.toString(),
      "delete",
      propertyType,
      undefined,
      `Deleted property type ${propertyType.name}`
    );
    
    // Delete the property type
    await ctx.db.delete(args.id);
    
    return {
      success: true,
      message: `Successfully deleted property type ${propertyType.name}`,
    };
  },
});
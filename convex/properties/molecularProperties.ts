/**
 * Functions for managing molecular properties
 */

import { mutation, query } from "../_generated/server";
import { v } from "convex/values";
import { Id } from "../_generated/dataModel";
import { MolecularProperty, CreateMolecularPropertyInput, UpdateMolecularPropertyInput, MolecularPropertyFilter, MolecularPropertyQueryOptions } from "./types";
import { validateCreateMolecularPropertyInput, validateUpdateMolecularPropertyInput } from "./validation";
import { getPropertyTypeById, getMolecularPropertyById, moleculePropertyExists, extractNumericValue, getMoleculeProperties } from "./helpers";
import { createAuditLog, getCurrentTimestamp, pagination, validateAuthenticated, validateRole } from "../utils/common";

/**
 * Create a new molecular property
 */
export const createMolecularProperty = mutation({
  args: {
    input: v.object({
      moleculeId: v.id("molecules"),
      propertyTypeId: v.id("propertyTypes"),
      value: v.union(v.string(), v.number(), v.boolean(), v.null()),
      numericValue: v.optional(v.number()),
      units: v.optional(v.string()),
      source: v.optional(v.id("dataSources")),
      calculationMethod: v.optional(v.string()),
      confidence: v.optional(v.number()),
    }),
  },
  handler: async (ctx, args): Promise<MolecularProperty> => {
    // Validate user authentication
    validateAuthenticated(ctx);
    
    // Get property type for validation
    const propertyType = await getPropertyTypeById(ctx.db, args.input.propertyTypeId);
    if (!propertyType) {
      throw new Error("Property type not found");
    }
    
    // Make sure molecule exists
    const molecule = await ctx.db.get(args.input.moleculeId);
    if (!molecule) {
      throw new Error("Molecule not found");
    }
    
    // Validate input based on property type
    const validationResult = validateCreateMolecularPropertyInput(
      args.input,
      propertyType.dataType as "string" | "number" | "boolean"
    );
    if (!validationResult.valid) {
      throw new Error(validationResult.error);
    }
    
    // Check if this property already exists for this molecule
    const propertyExists = await moleculePropertyExists(
      ctx.db,
      args.input.moleculeId,
      args.input.propertyTypeId
    );
    
    if (propertyExists) {
      throw new Error(
        `Property ${propertyType.name} already exists for this molecule`
      );
    }
    
    // Calculate numeric value if not provided
    let numericValue = args.input.numericValue;
    if (propertyType.dataType === "number" && numericValue === undefined) {
      numericValue = extractNumericValue(args.input.value);
    }
    
    // Create property document
    const now = getCurrentTimestamp();
    const propertyId = await ctx.db.insert("molecularProperties", {
      moleculeId: args.input.moleculeId,
      propertyTypeId: args.input.propertyTypeId,
      value: args.input.value,
      numericValue,
      units: args.input.units || propertyType.defaultUnits,
      source: args.input.source,
      calculationMethod: args.input.calculationMethod,
      confidence: args.input.confidence,
      createdAt: now,
      updatedAt: now,
    });
    
    // Fetch the inserted property
    const property = await ctx.db.get(propertyId) as MolecularProperty;
    
    // Create audit log
    await createAuditLog(
      ctx,
      "molecularProperties",
      propertyId.toString(),
      "create",
      undefined,
      property,
      `Created property ${propertyType.name} for molecule`
    );
    
    return property;
  },
});

/**
 * Update an existing molecular property
 */
export const updateMolecularProperty = mutation({
  args: {
    id: v.id("molecularProperties"),
    input: v.object({
      value: v.optional(v.union(v.string(), v.number(), v.boolean(), v.null())),
      numericValue: v.optional(v.number()),
      units: v.optional(v.string()),
      source: v.optional(v.id("dataSources")),
      calculationMethod: v.optional(v.string()),
      confidence: v.optional(v.number()),
    }),
  },
  handler: async (ctx, args): Promise<MolecularProperty> => {
    // Validate user authentication
    validateAuthenticated(ctx);
    
    // Get the current property
    const currentProperty = await getMolecularPropertyById(ctx.db, args.id);
    if (!currentProperty) {
      throw new Error("Molecular property not found");
    }
    
    // Get property type for validation
    const propertyType = await getPropertyTypeById(ctx.db, currentProperty.propertyTypeId);
    if (!propertyType) {
      throw new Error("Property type not found");
    }
    
    // Validate input based on property type
    const validationResult = validateUpdateMolecularPropertyInput(
      args.input,
      propertyType.dataType as "string" | "number" | "boolean"
    );
    if (!validationResult.valid) {
      throw new Error(validationResult.error);
    }
    
    // Calculate numeric value if value is updated but numeric value isn't
    let updateObject: any = { ...args.input };
    
    if (
      propertyType.dataType === "number" &&
      args.input.value !== undefined &&
      args.input.numericValue === undefined
    ) {
      updateObject.numericValue = extractNumericValue(args.input.value);
    }
    
    // Add update timestamp
    updateObject.updatedAt = getCurrentTimestamp();
    
    // Update the property
    const updatedProperty = await ctx.db.patch(args.id, updateObject) as MolecularProperty;
    
    // Create audit log
    await createAuditLog(
      ctx,
      "molecularProperties",
      args.id.toString(),
      "update",
      currentProperty,
      updatedProperty,
      `Updated property ${propertyType.name} for molecule`
    );
    
    return updatedProperty;
  },
});

/**
 * Get a molecular property by ID
 */
export const getMolecularProperty = query({
  args: {
    id: v.id("molecularProperties"),
  },
  handler: async (ctx, args): Promise<MolecularProperty | null> => {
    return await ctx.db.get(args.id) as MolecularProperty | null;
  },
});

/**
 * Get properties for a molecule
 */
export const getMoleculePropertiesById = query({
  args: {
    moleculeId: v.id("molecules"),
  },
  handler: async (ctx, args): Promise<Array<MolecularProperty & { propertyType: any }>> => {
    return await getMoleculeProperties(ctx.db, args.moleculeId);
  },
});

/**
 * Search for molecular properties based on filters
 */
export const searchMolecularProperties = query({
  args: {
    filter: v.optional(
      v.object({
        moleculeId: v.optional(v.id("molecules")),
        propertyTypeId: v.optional(v.id("propertyTypes")),
        valueEquals: v.optional(v.union(v.string(), v.number(), v.boolean())),
        valueRange: v.optional(
          v.object({
            min: v.optional(v.number()),
            max: v.optional(v.number()),
          })
        ),
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
    properties: MolecularProperty[];
    cursor: string | null;
    totalCount: number;
  }> => {
    const filter = args.filter || {};
    const options = args.options || {};
    const limit = options.limit || pagination.defaultPageSize;
    
    // Start building the query
    let query = ctx.db.query("molecularProperties");
    
    // Apply filters
    if (filter.moleculeId) {
      query = query.withIndex("by_molecule", (q) => 
        q.eq("moleculeId", filter.moleculeId)
      );
    }
    
    if (filter.propertyTypeId) {
      query = query.withIndex("by_property_type", (q) => 
        q.eq("propertyTypeId", filter.propertyTypeId)
      );
    }
    
    if (filter.moleculeId && filter.propertyTypeId) {
      query = query.withIndex("by_molecule_property", (q) => 
        q.eq("moleculeId", filter.moleculeId).eq("propertyTypeId", filter.propertyTypeId)
      );
    }
    
    if (filter.valueEquals !== undefined) {
      query = query.filter(q => q.eq(q.field("value"), filter.valueEquals));
    }
    
    if (filter.valueRange) {
      if (filter.valueRange.min !== undefined) {
        query = query.filter(q => q.gte(q.field("numericValue"), filter.valueRange!.min!));
      }
      
      if (filter.valueRange.max !== undefined) {
        query = query.filter(q => q.lte(q.field("numericValue"), filter.valueRange!.max!));
      }
    }
    
    // Apply pagination
    const results = await query
      .paginate({ cursor: options.cursor, numItems: limit });
    
    // Get total count (approximate)
    const totalCount = await query.collect();
    
    return {
      properties: results.page as MolecularProperty[],
      cursor: results.continueCursor,
      totalCount: totalCount.length,
    };
  },
});

/**
 * Delete a molecular property
 */
export const deleteMolecularProperty = mutation({
  args: {
    id: v.id("molecularProperties"),
  },
  handler: async (ctx, args): Promise<{ success: boolean; message: string }> => {
    // Validate user authentication
    validateAuthenticated(ctx);
    
    // Get the current property
    const property = await getMolecularPropertyById(ctx.db, args.id);
    if (!property) {
      throw new Error("Molecular property not found");
    }
    
    // Get property type
    const propertyType = await getPropertyTypeById(ctx.db, property.propertyTypeId);
    
    // Create audit log before deletion
    await createAuditLog(
      ctx,
      "molecularProperties",
      args.id.toString(),
      "delete",
      property,
      undefined,
      `Deleted property ${propertyType?.name || "unknown"} for molecule`
    );
    
    // Delete the property
    await ctx.db.delete(args.id);
    
    return {
      success: true,
      message: `Successfully deleted property`,
    };
  },
});

/**
 * Create or update multiple properties for a molecule in batch
 */
export const batchUpdateMoleculeProperties = mutation({
  args: {
    moleculeId: v.id("molecules"),
    properties: v.array(
      v.object({
        propertyTypeId: v.id("propertyTypes"),
        value: v.union(v.string(), v.number(), v.boolean(), v.null()),
        units: v.optional(v.string()),
        source: v.optional(v.id("dataSources")),
        calculationMethod: v.optional(v.string()),
        confidence: v.optional(v.number()),
      })
    ),
  },
  handler: async (ctx, args): Promise<{
    created: number;
    updated: number;
    failed: number;
    errors: Array<{ propertyTypeId: Id<"propertyTypes">; error: string }>;
  }> => {
    // Validate user authentication
    validateAuthenticated(ctx);
    
    // Make sure molecule exists
    const molecule = await ctx.db.get(args.moleculeId);
    if (!molecule) {
      throw new Error("Molecule not found");
    }
    
    const results = {
      created: 0,
      updated: 0,
      failed: 0,
      errors: [] as Array<{ propertyTypeId: Id<"propertyTypes">; error: string }>,
    };
    
    // Process each property in the batch
    for (const propertyData of args.properties) {
      try {
        // Get property type for validation
        const propertyType = await getPropertyTypeById(ctx.db, propertyData.propertyTypeId);
        if (!propertyType) {
          throw new Error("Property type not found");
        }
        
        // Check if property already exists
        const existingProperty = await ctx.db
          .query("molecularProperties")
          .withIndex("by_molecule_property", (q) => 
            q.eq("moleculeId", args.moleculeId).eq("propertyTypeId", propertyData.propertyTypeId)
          )
          .first();
        
        if (existingProperty) {
          // Update existing property
          const input: UpdateMolecularPropertyInput = {
            value: propertyData.value,
            units: propertyData.units,
            source: propertyData.source,
            calculationMethod: propertyData.calculationMethod,
            confidence: propertyData.confidence,
          };
          
          // Validate update
          const validationResult = validateUpdateMolecularPropertyInput(
            input,
            propertyType.dataType as "string" | "number" | "boolean"
          );
          
          if (!validationResult.valid) {
            throw new Error(validationResult.error);
          }
          
          // Calculate numeric value if needed
          if (
            propertyType.dataType === "number" &&
            input.numericValue === undefined
          ) {
            input.numericValue = extractNumericValue(input.value);
          }
          
          // Update the property
          await ctx.db.patch(existingProperty._id, {
            ...input,
            updatedAt: getCurrentTimestamp(),
          });
          
          results.updated++;
        } else {
          // Create new property
          const input: CreateMolecularPropertyInput = {
            moleculeId: args.moleculeId,
            propertyTypeId: propertyData.propertyTypeId,
            value: propertyData.value,
            units: propertyData.units,
            source: propertyData.source,
            calculationMethod: propertyData.calculationMethod,
            confidence: propertyData.confidence,
          };
          
          // Validate create
          const validationResult = validateCreateMolecularPropertyInput(
            input,
            propertyType.dataType as "string" | "number" | "boolean"
          );
          
          if (!validationResult.valid) {
            throw new Error(validationResult.error);
          }
          
          // Calculate numeric value if needed
          if (
            propertyType.dataType === "number" &&
            input.numericValue === undefined
          ) {
            input.numericValue = extractNumericValue(input.value);
          }
          
          // Create the property
          const now = getCurrentTimestamp();
          await ctx.db.insert("molecularProperties", {
            ...input,
            createdAt: now,
            updatedAt: now,
          });
          
          results.created++;
        }
      } catch (error) {
        results.failed++;
        results.errors.push({
          propertyTypeId: propertyData.propertyTypeId,
          error: error instanceof Error ? error.message : "Unknown error",
        });
      }
    }
    
    return results;
  },
});
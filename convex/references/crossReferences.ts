/**
 * CRUD operations for molecule cross references
 */

import { mutation, query } from "../_generated/server";
import { v } from "convex/values";
import { Id } from "../_generated/dataModel";
import { 
  MoleculeCrossReference, 
  CreateCrossReferenceInput, 
  UpdateCrossReferenceInput,
  CrossReferenceFilter,
  CrossReferenceQueryOptions
} from "./types";
import { 
  validateCreateCrossReferenceInput, 
  validateUpdateCrossReferenceInput,
  validateCrossReferenceExists,
  validateMoleculeExists
} from "./validation";
import { 
  formatCrossReferenceUrl,
  expandCrossReferenceWithMolecule
} from "./helpers";
import { getCurrentUser } from "../auth/users";

/**
 * Create a new molecule cross reference
 */
export const createCrossReference = mutation({
  args: {
    crossRef: v.object({
      moleculeId: v.id("molecules"),
      databaseName: v.string(),
      identifier: v.string(),
      url: v.optional(v.string())
    })
  },
  handler: async (ctx, args) => {
    // Validate input
    validateCreateCrossReferenceInput(args.crossRef);
    validateMoleculeExists(args.crossRef.moleculeId);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // User must be logged in to create a cross reference
    if (!userId) {
      throw new Error("Authentication required to create a cross reference");
    }
    
    // Verify molecule exists
    const molecule = await ctx.db.get(args.crossRef.moleculeId);
    if (!molecule) {
      throw new Error(`Molecule with ID ${args.crossRef.moleculeId} not found`);
    }
    
    // Generate URL if not provided
    let url = args.crossRef.url;
    if (!url) {
      url = formatCrossReferenceUrl(args.crossRef.databaseName, args.crossRef.identifier);
    }
    
    // Prepare data for insertion
    const now = Date.now();
    const crossRefData = {
      ...args.crossRef,
      url,
      createdAt: now,
      updatedAt: now
    };
    
    // Check for duplicates
    const existingRefs = await ctx.db
      .query("moleculeCrossReferences")
      .withIndex("by_molecule", q => q.eq("moleculeId", args.crossRef.moleculeId))
      .filter(q => 
        q.and(
          q.eq(q.field("databaseName"), args.crossRef.databaseName),
          q.eq(q.field("identifier"), args.crossRef.identifier)
        )
      )
      .collect();
    
    if (existingRefs.length > 0) {
      throw new Error(
        `Cross reference for ${args.crossRef.databaseName} with identifier ${args.crossRef.identifier} already exists for this molecule`
      );
    }
    
    // Insert cross reference
    const crossRefId = await ctx.db.insert("moleculeCrossReferences", crossRefData);
    
    // Create audit log entry
    await ctx.db.insert("scientificDataAudit", {
      table: "moleculeCrossReferences",
      documentId: crossRefId,
      operation: "create",
      userId,
      newValue: crossRefData,
      timestamp: now,
    });
    
    return crossRefId;
  }
});

/**
 * Get a cross reference by ID
 */
export const getCrossReference = query({
  args: {
    crossRefId: v.id("moleculeCrossReferences"),
    includeMolecule: v.optional(v.boolean())
  },
  handler: async (ctx, args) => {
    // Validate cross reference exists
    validateCrossReferenceExists(args.crossRefId);
    
    // Get the cross reference
    const crossRef = await ctx.db.get(args.crossRefId);
    if (!crossRef) {
      return null;
    }
    
    // Include molecule data if requested
    if (args.includeMolecule) {
      return expandCrossReferenceWithMolecule(ctx, crossRef);
    }
    
    return crossRef;
  }
});

/**
 * Update an existing cross reference
 */
export const updateCrossReference = mutation({
  args: {
    crossRefId: v.id("moleculeCrossReferences"),
    update: v.object({
      databaseName: v.optional(v.string()),
      identifier: v.optional(v.string()),
      url: v.optional(v.string())
    })
  },
  handler: async (ctx, args) => {
    // Validate input
    validateUpdateCrossReferenceInput(args.update);
    validateCrossReferenceExists(args.crossRefId);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // User must be logged in to update a cross reference
    if (!userId) {
      throw new Error("Authentication required to update a cross reference");
    }
    
    // Get existing cross reference
    const existingCrossRef = await ctx.db.get(args.crossRefId);
    if (!existingCrossRef) {
      throw new Error(`Cross reference with ID ${args.crossRefId} not found`);
    }
    
    // Prepare update data
    let updateData = { ...args.update };
    
    // If database name or identifier is changing, update the URL if not explicitly provided
    if ((args.update.databaseName || args.update.identifier) && !args.update.url) {
      const databaseName = args.update.databaseName || existingCrossRef.databaseName;
      const identifier = args.update.identifier || existingCrossRef.identifier;
      
      updateData.url = formatCrossReferenceUrl(
        databaseName, 
        identifier, 
        existingCrossRef.url
      );
    }
    
    // Add updatedAt timestamp
    const now = Date.now();
    updateData = {
      ...updateData,
      updatedAt: now
    };
    
    // Check for duplicates if changing database or identifier
    if (args.update.databaseName || args.update.identifier) {
      const databaseName = args.update.databaseName || existingCrossRef.databaseName;
      const identifier = args.update.identifier || existingCrossRef.identifier;
      
      const existingRefs = await ctx.db
        .query("moleculeCrossReferences")
        .withIndex("by_molecule", q => q.eq("moleculeId", existingCrossRef.moleculeId))
        .filter(q => 
          q.and(
            q.eq(q.field("databaseName"), databaseName),
            q.eq(q.field("identifier"), identifier),
            q.neq(q.field("_id"), args.crossRefId)
          )
        )
        .collect();
      
      if (existingRefs.length > 0) {
        throw new Error(
          `Cross reference for ${databaseName} with identifier ${identifier} already exists for this molecule`
        );
      }
    }
    
    // Update cross reference
    await ctx.db.patch(args.crossRefId, updateData);
    
    // Create audit log entry
    await ctx.db.insert("scientificDataAudit", {
      table: "moleculeCrossReferences",
      documentId: args.crossRefId,
      operation: "update",
      userId,
      previousValue: existingCrossRef,
      newValue: { ...existingCrossRef, ...updateData },
      timestamp: now
    });
    
    return args.crossRefId;
  }
});

/**
 * Delete a cross reference
 */
export const deleteCrossReference = mutation({
  args: {
    crossRefId: v.id("moleculeCrossReferences")
  },
  handler: async (ctx, args) => {
    // Validate cross reference exists
    validateCrossReferenceExists(args.crossRefId);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // User must be logged in to delete a cross reference
    if (!userId) {
      throw new Error("Authentication required to delete a cross reference");
    }
    
    // Get existing cross reference
    const existingCrossRef = await ctx.db.get(args.crossRefId);
    if (!existingCrossRef) {
      throw new Error(`Cross reference with ID ${args.crossRefId} not found`);
    }
    
    // Create audit log entry
    const now = Date.now();
    await ctx.db.insert("scientificDataAudit", {
      table: "moleculeCrossReferences",
      documentId: args.crossRefId,
      operation: "delete",
      userId,
      previousValue: existingCrossRef,
      timestamp: now
    });
    
    // Delete the cross reference
    await ctx.db.delete(args.crossRefId);
    
    return true;
  }
});

/**
 * List cross references with optional filtering
 */
export const listCrossReferences = query({
  args: {
    filter: v.optional(v.object({
      moleculeId: v.optional(v.id("molecules")),
      databaseName: v.optional(v.string()),
      includeMolecule: v.optional(v.boolean())
    })),
    options: v.optional(v.object({
      limit: v.optional(v.number()),
      cursor: v.optional(v.string())
    }))
  },
  handler: async (ctx, args) => {
    // Set up the query
    let query = ctx.db.query("moleculeCrossReferences");
    
    // Apply filters
    if (args.filter) {
      if (args.filter.moleculeId) {
        query = query.withIndex("by_molecule", q => 
          q.eq("moleculeId", args.filter!.moleculeId!)
        );
      }
      
      if (args.filter.databaseName) {
        query = query.withIndex("by_database", q => 
          q.eq("databaseName", args.filter!.databaseName!)
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
    const crossRefs = await query.collect();
    
    // Include molecule data if requested
    if (args.filter?.includeMolecule) {
      const crossRefsWithMolecules = [];
      
      for (const crossRef of crossRefs) {
        crossRefsWithMolecules.push(await expandCrossReferenceWithMolecule(ctx, crossRef));
      }
      
      return crossRefsWithMolecules;
    }
    
    return crossRefs;
  }
});

/**
 * Get all cross references for a molecule
 */
export const getMoleculeCrossReferences = query({
  args: {
    moleculeId: v.id("molecules")
  },
  handler: async (ctx, args) => {
    // Validate molecule exists
    validateMoleculeExists(args.moleculeId);
    
    // Query for cross references
    const crossRefs = await ctx.db
      .query("moleculeCrossReferences")
      .withIndex("by_molecule", q => q.eq("moleculeId", args.moleculeId))
      .collect();
    
    return crossRefs;
  }
});

/**
 * Batch add multiple cross references for a molecule
 */
export const batchAddCrossReferences = mutation({
  args: {
    moleculeId: v.id("molecules"),
    crossRefs: v.array(
      v.object({
        databaseName: v.string(),
        identifier: v.string(),
        url: v.optional(v.string())
      })
    )
  },
  handler: async (ctx, args) => {
    // Validate molecule exists
    validateMoleculeExists(args.moleculeId);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // User must be logged in to add cross references
    if (!userId) {
      throw new Error("Authentication required to add cross references");
    }
    
    // Verify molecule exists
    const molecule = await ctx.db.get(args.moleculeId);
    if (!molecule) {
      throw new Error(`Molecule with ID ${args.moleculeId} not found`);
    }
    
    // Record time for all operations
    const now = Date.now();
    
    // Get existing cross references to avoid duplicates
    const existingRefs = await ctx.db
      .query("moleculeCrossReferences")
      .withIndex("by_molecule", q => q.eq("moleculeId", args.moleculeId))
      .collect();
    
    // Create a set for checking duplicates
    const existingRefSet = new Set(
      existingRefs.map(ref => `${ref.databaseName}:${ref.identifier}`)
    );
    
    // Validate and insert each cross reference
    const crossRefIds: Id<"moleculeCrossReferences">[] = [];
    
    for (const crossRef of args.crossRefs) {
      // Validate cross reference input
      validateCreateCrossReferenceInput({
        ...crossRef,
        moleculeId: args.moleculeId
      });
      
      // Check for duplicates
      const refKey = `${crossRef.databaseName}:${crossRef.identifier}`;
      if (existingRefSet.has(refKey)) {
        continue; // Skip duplicates
      }
      
      // Generate URL if not provided
      let url = crossRef.url;
      if (!url) {
        url = formatCrossReferenceUrl(crossRef.databaseName, crossRef.identifier);
      }
      
      // Prepare cross reference data
      const crossRefData = {
        ...crossRef,
        moleculeId: args.moleculeId,
        url,
        createdAt: now,
        updatedAt: now
      };
      
      // Insert cross reference
      const crossRefId = await ctx.db.insert("moleculeCrossReferences", crossRefData);
      
      // Create audit log entry
      await ctx.db.insert("scientificDataAudit", {
        table: "moleculeCrossReferences",
        documentId: crossRefId,
        operation: "create",
        userId,
        newValue: crossRefData,
        timestamp: now
      });
      
      crossRefIds.push(crossRefId);
      existingRefSet.add(refKey); // Add to set to prevent duplicates in batch
    }
    
    return crossRefIds;
  }
});
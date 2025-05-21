/**
 * Functions for updating molecules
 */

import { mutation } from "../_generated/server";
import { v } from "convex/values";
import { Id } from "../_generated/dataModel";
import { Molecule, MoleculeConsolidationResult } from "./types";
import { validateUpdateMoleculeInput } from "./validation";
import { getMoleculeById, markMoleculeAsConsolidated } from "./helpers";
import { createAuditLog, getCurrentTimestamp, validateAuthenticated, validateRole } from "../utils/common";

/**
 * Update a molecule
 */
export const updateMolecule = mutation({
  args: {
    id: v.id("molecules"),
    input: v.object({
      name: v.optional(v.string()),
      pubchemCid: v.optional(v.string()),
      canonicalSmiles: v.optional(v.string()),
      inchiKey: v.optional(v.string()),
      formula: v.optional(v.string()),
      status: v.optional(v.union(
        v.literal("active"),
        v.literal("deprecated"),
        v.literal("consolidated")
      )),
      dataSource: v.optional(v.id("dataSources")),
      sourceId: v.optional(v.string()),
    }),
  },
  handler: async (ctx, args): Promise<Molecule> => {
    // Validate user authentication
    validateAuthenticated(ctx);
    
    // Validate input
    const validationResult = validateUpdateMoleculeInput(args.input);
    if (!validationResult.valid) {
      throw new Error(validationResult.error);
    }
    
    // Get the current molecule
    const currentMolecule = await getMoleculeById(ctx.db, args.id);
    if (!currentMolecule) {
      throw new Error("Molecule not found");
    }
    
    // If molecule is consolidated, prevent updates
    if (currentMolecule.consolidated) {
      throw new Error("Cannot update a consolidated molecule");
    }
    
    // Create update object
    const updateObject = {
      ...args.input,
      updatedAt: getCurrentTimestamp(),
    };
    
    // Update the molecule
    const updatedMolecule = await ctx.db.patch(args.id, updateObject) as Molecule;
    
    // Create audit log
    await createAuditLog(
      ctx,
      "molecules",
      args.id.toString(),
      "update",
      currentMolecule,
      updatedMolecule,
      "User update"
    );
    
    return updatedMolecule;
  },
});

/**
 * Mark a molecule as deprecated
 */
export const deprecateMolecule = mutation({
  args: {
    id: v.id("molecules"),
    reason: v.optional(v.string()),
  },
  handler: async (ctx, args): Promise<Molecule> => {
    // Validate user authentication and role
    await validateRole(ctx, ["admin", "scientist"]);
    
    // Get the current molecule
    const currentMolecule = await getMoleculeById(ctx.db, args.id);
    if (!currentMolecule) {
      throw new Error("Molecule not found");
    }
    
    // If molecule is consolidated, prevent updates
    if (currentMolecule.consolidated) {
      throw new Error("Cannot deprecate a consolidated molecule");
    }
    
    // Update the molecule
    const updatedMolecule = await ctx.db.patch(args.id, {
      status: "deprecated",
      updatedAt: getCurrentTimestamp(),
    }) as Molecule;
    
    // Create audit log
    await createAuditLog(
      ctx,
      "molecules",
      args.id.toString(),
      "update",
      currentMolecule,
      updatedMolecule,
      args.reason || "Marked as deprecated"
    );
    
    return updatedMolecule;
  },
});

/**
 * Consolidate multiple molecules into a primary molecule
 */
export const consolidateMolecules = mutation({
  args: {
    primaryId: v.id("molecules"),
    moleculeIds: v.array(v.id("molecules")),
    reason: v.optional(v.string()),
  },
  handler: async (ctx, args): Promise<MoleculeConsolidationResult> => {
    // Validate user authentication and role
    await validateRole(ctx, ["admin", "scientist"]);
    
    // Get the primary molecule
    const primaryMolecule = await getMoleculeById(ctx.db, args.primaryId);
    if (!primaryMolecule) {
      throw new Error("Primary molecule not found");
    }
    
    // Cannot use a consolidated molecule as primary
    if (primaryMolecule.consolidated) {
      throw new Error("Cannot use a consolidated molecule as primary");
    }
    
    // Remove primary ID from list if present
    const secondaryIds = args.moleculeIds.filter(
      id => id.toString() !== args.primaryId.toString()
    );
    
    // Validate all secondary molecules exist
    const secondaryMolecules = await Promise.all(
      secondaryIds.map(id => getMoleculeById(ctx.db, id))
    );
    
    if (secondaryMolecules.some(m => m === null)) {
      throw new Error("One or more secondary molecules not found");
    }
    
    // Track consolidation results
    const result: MoleculeConsolidationResult = {
      primaryId: args.primaryId,
      consolidatedIds: [],
      totalProperties: 0,
      mergedProperties: 0,
    };
    
    // Process each secondary molecule
    for (const moleculeId of secondaryIds) {
      // Mark molecule as consolidated
      await markMoleculeAsConsolidated(ctx.db, moleculeId, args.primaryId);
      
      // Get all properties for this molecule
      const properties = await ctx.db
        .query("molecularProperties")
        .withIndex("by_molecule", q => q.eq("moleculeId", moleculeId))
        .collect();
      
      result.totalProperties += properties.length;
      
      // For each property, check if primary molecule already has it
      for (const property of properties) {
        const existingProperty = await ctx.db
          .query("molecularProperties")
          .withIndex("by_molecule_property", q => 
            q.eq("moleculeId", args.primaryId).eq("propertyTypeId", property.propertyTypeId)
          )
          .first();
        
        if (!existingProperty) {
          // If no existing property, create one for the primary molecule
          await ctx.db.insert("molecularProperties", {
            moleculeId: args.primaryId,
            propertyTypeId: property.propertyTypeId,
            value: property.value,
            numericValue: property.numericValue,
            units: property.units,
            source: property.source,
            calculationMethod: property.calculationMethod,
            confidence: property.confidence,
            createdAt: getCurrentTimestamp(),
            updatedAt: getCurrentTimestamp(),
          });
          
          result.mergedProperties++;
        }
      }
      
      // Add to consolidated IDs
      result.consolidatedIds.push(moleculeId);
      
      // Create audit log for the consolidation
      await createAuditLog(
        ctx,
        "molecules",
        moleculeId.toString(),
        "update",
        await ctx.db.get(moleculeId),
        { consolidated: true, consolidatedWith: args.primaryId },
        args.reason || `Consolidated with ${primaryMolecule.name}`
      );
    }
    
    return result;
  },
});
/**
 * Functions for deleting molecules
 */

import { mutation } from "../_generated/server";
import { v } from "convex/values";
import { Id } from "../_generated/dataModel";
import { getMoleculeById } from "./helpers";
import { createAuditLog, validateRole } from "../utils/common";

/**
 * Delete a molecule
 * 
 * Note: This operation is restricted to admin users only
 * and will fail if the molecule is referenced in mixtures or has properties
 */
export const deleteMolecule = mutation({
  args: {
    id: v.id("molecules"),
    force: v.optional(v.boolean()),
    reason: v.optional(v.string()),
  },
  handler: async (ctx, args): Promise<{ success: boolean; message: string }> => {
    // Validate user has admin role
    await validateRole(ctx, ["admin"]);
    
    // Get the current molecule
    const currentMolecule = await getMoleculeById(ctx.db, args.id);
    if (!currentMolecule) {
      throw new Error("Molecule not found");
    }
    
    // Check for references to this molecule if not force deleting
    if (!args.force) {
      // Check for properties
      const properties = await ctx.db
        .query("molecularProperties")
        .withIndex("by_molecule", q => q.eq("moleculeId", args.id))
        .first();
      
      if (properties) {
        throw new Error(
          "Cannot delete molecule with properties. Use force=true to delete anyway."
        );
      }
      
      // Check for mixture components
      const mixtureComponents = await ctx.db
        .query("mixtureComponents")
        .withIndex("by_molecule", q => q.eq("moleculeId", args.id))
        .first();
      
      if (mixtureComponents) {
        throw new Error(
          "Cannot delete molecule used in mixtures. Use force=true to delete anyway."
        );
      }
      
      // Check for predictions
      const predictions = await ctx.db
        .query("predictions")
        .withIndex("by_molecule", q => q.eq("moleculeId", args.id))
        .first();
      
      if (predictions) {
        throw new Error(
          "Cannot delete molecule with predictions. Use force=true to delete anyway."
        );
      }
      
      // Check for toxicity data
      const toxicityData = await ctx.db
        .query("toxicityData")
        .withIndex("by_molecule", q => q.eq("moleculeId", args.id))
        .first();
      
      if (toxicityData) {
        throw new Error(
          "Cannot delete molecule with toxicity data. Use force=true to delete anyway."
        );
      }
      
      // Check for cross references
      const crossReferences = await ctx.db
        .query("moleculeCrossReferences")
        .withIndex("by_molecule", q => q.eq("moleculeId", args.id))
        .first();
      
      if (crossReferences) {
        throw new Error(
          "Cannot delete molecule with cross references. Use force=true to delete anyway."
        );
      }
      
      // Check for synonyms
      const synonyms = await ctx.db
        .query("moleculeSynonyms")
        .withIndex("by_molecule", q => q.eq("moleculeId", args.id))
        .first();
      
      if (synonyms) {
        throw new Error(
          "Cannot delete molecule with synonyms. Use force=true to delete anyway."
        );
      }
    }
    
    // If force is true, delete all related records
    if (args.force) {
      // Delete properties
      const properties = await ctx.db
        .query("molecularProperties")
        .withIndex("by_molecule", q => q.eq("moleculeId", args.id))
        .collect();
      
      for (const property of properties) {
        await ctx.db.delete(property._id);
      }
      
      // Delete mixture components
      const mixtureComponents = await ctx.db
        .query("mixtureComponents")
        .withIndex("by_molecule", q => q.eq("moleculeId", args.id))
        .collect();
      
      for (const component of mixtureComponents) {
        await ctx.db.delete(component._id);
      }
      
      // Delete predictions
      const predictions = await ctx.db
        .query("predictions")
        .withIndex("by_molecule", q => q.eq("moleculeId", args.id))
        .collect();
      
      for (const prediction of predictions) {
        await ctx.db.delete(prediction._id);
      }
      
      // Delete toxicity data
      const toxicityData = await ctx.db
        .query("toxicityData")
        .withIndex("by_molecule", q => q.eq("moleculeId", args.id))
        .collect();
      
      for (const data of toxicityData) {
        await ctx.db.delete(data._id);
      }
      
      // Delete cross references
      const crossReferences = await ctx.db
        .query("moleculeCrossReferences")
        .withIndex("by_molecule", q => q.eq("moleculeId", args.id))
        .collect();
      
      for (const ref of crossReferences) {
        await ctx.db.delete(ref._id);
      }
      
      // Delete synonyms
      const synonyms = await ctx.db
        .query("moleculeSynonyms")
        .withIndex("by_molecule", q => q.eq("moleculeId", args.id))
        .collect();
      
      for (const synonym of synonyms) {
        await ctx.db.delete(synonym._id);
      }
    }
    
    // Create audit log before deletion
    await createAuditLog(
      ctx,
      "molecules",
      args.id.toString(),
      "delete",
      currentMolecule,
      undefined,
      args.reason || `Deleted molecule ${currentMolecule.name}`
    );
    
    // Delete the molecule
    await ctx.db.delete(args.id);
    
    return {
      success: true,
      message: `Successfully deleted molecule ${currentMolecule.name}`,
    };
  },
});

/**
 * Delete multiple molecules in a batch
 * 
 * Note: This operation is restricted to admin users only
 */
export const deleteMoleculesBatch = mutation({
  args: {
    ids: v.array(v.id("molecules")),
    force: v.optional(v.boolean()),
    reason: v.optional(v.string()),
  },
  handler: async (ctx, args): Promise<{
    success: number;
    failed: number;
    errors: Array<{ id: Id<"molecules">; error: string }>;
  }> => {
    // Validate user has admin role
    await validateRole(ctx, ["admin"]);
    
    const results = {
      success: 0,
      failed: 0,
      errors: [] as Array<{ id: Id<"molecules">; error: string }>,
    };
    
    // Process each molecule in the batch
    for (const id of args.ids) {
      try {
        // Call the individual delete function
        await deleteMolecule.handler(ctx, {
          id,
          force: args.force,
          reason: args.reason,
        });
        
        results.success++;
      } catch (error) {
        results.failed++;
        results.errors.push({
          id,
          error: error instanceof Error ? error.message : "Unknown error",
        });
      }
    }
    
    return results;
  },
});
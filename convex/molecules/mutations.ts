import { mutation } from "../_generated/server";
import { v } from "convex/values";

/**
 * Create a new molecule
 */
export const create = mutation({
  args: {
    name: v.string(),
    pubchemCid: v.optional(v.string()),
    canonicalSmiles: v.optional(v.string()),
    inchiKey: v.optional(v.string()),
    formula: v.optional(v.string()),
    status: v.string(),
  },
  handler: async (ctx, args) => {
    const timestamp = Date.now();
    
    return ctx.db.insert("molecules", {
      name: args.name,
      pubchemCid: args.pubchemCid,
      canonicalSmiles: args.canonicalSmiles,
      inchiKey: args.inchiKey,
      formula: args.formula,
      status: args.status,
      createdAt: timestamp,
      updatedAt: timestamp,
    });
  },
});

/**
 * Update an existing molecule
 */
export const update = mutation({
  args: {
    id: v.id("molecules"),
    name: v.optional(v.string()),
    pubchemCid: v.optional(v.string()),
    canonicalSmiles: v.optional(v.string()),
    inchiKey: v.optional(v.string()),
    formula: v.optional(v.string()),
    status: v.optional(v.string()),
  },
  handler: async (ctx, args) => {
    const { id, ...updates } = args;
    
    // Get current molecule to ensure it exists
    const molecule = await ctx.db.get(id);
    if (!molecule) {
      throw new Error(`Molecule with ID ${id} not found`);
    }
    
    // Update the molecule with the specified fields
    return ctx.db.patch(id, {
      ...updates,
      updatedAt: Date.now(),
    });
  },
});

/**
 * Remove a molecule (soft delete by marking as deprecated)
 */
export const remove = mutation({
  args: {
    id: v.id("molecules"),
  },
  handler: async (ctx, args) => {
    // Mark the molecule as deprecated instead of hard deleting
    return ctx.db.patch(args.id, { 
      status: "deprecated",
      updatedAt: Date.now() 
    });
  },
});
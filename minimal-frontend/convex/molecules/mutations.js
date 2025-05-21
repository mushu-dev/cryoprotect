import { mutation } from "../_generated/server";
import { v } from "convex/values";

/**
 * Add a new molecule
 */
export const addMolecule = mutation({
  args: {
    name: v.string(),
    formula: v.optional(v.string()),
    pubchemCid: v.optional(v.string()),
    smiles: v.optional(v.string()),
    molecularWeight: v.optional(v.number()),
    isCryoprotectant: v.optional(v.boolean()),
    description: v.optional(v.string()),
  },
  handler: async (ctx, args) => {
    const now = Date.now();
    
    const moleculeId = await ctx.db.insert("molecules", {
      ...args,
      createdAt: now,
      updatedAt: now,
    });
    
    return moleculeId;
  },
});

/**
 * Update an existing molecule
 */
export const updateMolecule = mutation({
  args: {
    id: v.id("molecules"),
    name: v.optional(v.string()),
    formula: v.optional(v.string()),
    pubchemCid: v.optional(v.string()),
    smiles: v.optional(v.string()),
    molecularWeight: v.optional(v.number()),
    isCryoprotectant: v.optional(v.boolean()),
    description: v.optional(v.string()),
  },
  handler: async (ctx, args) => {
    const { id, ...updates } = args;
    
    const existing = await ctx.db.get(id);
    if (!existing) {
      throw new Error(`Molecule with ID ${id} not found`);
    }
    
    await ctx.db.patch(id, {
      ...updates,
      updatedAt: Date.now(),
    });
    
    return id;
  },
});

/**
 * Delete a molecule
 */
export const deleteMolecule = mutation({
  args: { id: v.id("molecules") },
  handler: async (ctx, args) => {
    const existing = await ctx.db.get(args.id);
    if (!existing) {
      throw new Error(`Molecule with ID ${args.id} not found`);
    }
    
    await ctx.db.delete(args.id);
    
    // Also delete any properties associated with this molecule
    const properties = await ctx.db.query("molecularProperties")
      .withIndex("by_molecule", q => q.eq("moleculeId", args.id))
      .collect();
    
    for (const property of properties) {
      await ctx.db.delete(property._id);
    }
    
    return args.id;
  },
});
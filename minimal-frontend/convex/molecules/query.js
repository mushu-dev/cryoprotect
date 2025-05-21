import { query } from "../_generated/server";
import { v } from "convex/values";

/**
 * Get all recent molecules
 */
export const getRecentMolecules = query({
  handler: async (ctx) => {
    return await ctx.db.query("molecules")
      .order("desc")
      .take(20);
  },
});

/**
 * Get a molecule by ID
 */
export const getMolecule = query({
  args: { id: v.id("molecules") },
  handler: async (ctx, args) => {
    const molecule = await ctx.db.get(args.id);
    
    if (!molecule) {
      return null;
    }
    
    // Get properties for this molecule
    const properties = await ctx.db.query("molecularProperties")
      .withIndex("by_molecule", q => q.eq("moleculeId", args.id))
      .collect();
    
    return {
      ...molecule,
      properties
    };
  },
});

/**
 * Get molecules by IDs
 */
export const getMoleculesByIds = query({
  args: { ids: v.array(v.id("molecules")) },
  handler: async (ctx, args) => {
    return await Promise.all(
      args.ids.map(async (id) => {
        return await ctx.db.get(id);
      })
    );
  },
});

/**
 * Search molecules by name
 */
export const searchMolecules = query({
  args: { query: v.optional(v.string()) },
  handler: async (ctx, args) => {
    let molecules = ctx.db.query("molecules");
    
    if (args.query) {
      molecules = molecules.withSearchIndex("by_name", q => 
        q.search("name", args.query)
      );
    }
    
    return await molecules.take(100);
  },
});
import { query } from "../_generated/server";
import { v } from "convex/values";

/**
 * Get all recent mixtures
 */
export const getRecentMixtures = query({
  handler: async (ctx) => {
    return await ctx.db.query("mixtures")
      .order("desc")
      .take(20);
  },
});

/**
 * Get a mixture by ID
 */
export const getMixture = query({
  args: { id: v.id("mixtures") },
  handler: async (ctx, args) => {
    const mixture = await ctx.db.get(args.id);
    
    if (!mixture) {
      return null;
    }
    
    // Get components for this mixture
    const componentDocs = await ctx.db.query("mixtureComponents")
      .withIndex("by_mixture", q => q.eq("mixtureId", args.id))
      .collect();
    
    // Get molecule data for each component
    const components = await Promise.all(
      componentDocs.map(async (comp) => {
        const molecule = await ctx.db.get(comp.moleculeId);
        return {
          ...comp,
          molecule
        };
      })
    );
    
    return {
      ...mixture,
      components
    };
  },
});

/**
 * Get mixtures by IDs
 */
export const getMixturesByIds = query({
  args: { ids: v.array(v.id("mixtures")) },
  handler: async (ctx, args) => {
    return await Promise.all(
      args.ids.map(async (id) => {
        return await ctx.db.get(id);
      })
    );
  },
});

/**
 * Search mixtures by name
 */
export const searchMixtures = query({
  args: { query: v.optional(v.string()) },
  handler: async (ctx, args) => {
    let mixtures = ctx.db.query("mixtures");
    
    if (args.query) {
      mixtures = mixtures.withSearchIndex("by_name", q => 
        q.search("name", args.query)
      );
    }
    
    return await mixtures.take(100);
  },
});
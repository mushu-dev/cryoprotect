import { mutation } from "../_generated/server";
import { v } from "convex/values";

/**
 * Add a new mixture
 */
export const addMixture = mutation({
  args: {
    name: v.string(),
    description: v.optional(v.string()),
    freezingPoint: v.optional(v.number()),
    components: v.optional(v.array(v.object({
      moleculeId: v.id("molecules"),
      concentration: v.number(),
      units: v.string(),
      role: v.optional(v.string()),
    }))),
  },
  handler: async (ctx, args) => {
    const { components, ...mixtureData } = args;
    const now = Date.now();
    
    // Insert the mixture
    const mixtureId = await ctx.db.insert("mixtures", {
      ...mixtureData,
      createdAt: now,
      updatedAt: now,
    });
    
    // Insert components if provided
    if (components && components.length > 0) {
      for (const component of components) {
        await ctx.db.insert("mixtureComponents", {
          mixtureId,
          moleculeId: component.moleculeId,
          concentration: component.concentration,
          units: component.units,
          role: component.role,
          createdAt: now,
        });
      }
    }
    
    return mixtureId;
  },
});

/**
 * Update an existing mixture
 */
export const updateMixture = mutation({
  args: {
    id: v.id("mixtures"),
    name: v.optional(v.string()),
    description: v.optional(v.string()),
    freezingPoint: v.optional(v.number()),
  },
  handler: async (ctx, args) => {
    const { id, ...updates } = args;
    
    const existing = await ctx.db.get(id);
    if (!existing) {
      throw new Error(`Mixture with ID ${id} not found`);
    }
    
    await ctx.db.patch(id, {
      ...updates,
      updatedAt: Date.now(),
    });
    
    return id;
  },
});

/**
 * Delete a mixture
 */
export const deleteMixture = mutation({
  args: { id: v.id("mixtures") },
  handler: async (ctx, args) => {
    const existing = await ctx.db.get(args.id);
    if (!existing) {
      throw new Error(`Mixture with ID ${args.id} not found`);
    }
    
    await ctx.db.delete(args.id);
    
    // Also delete any components associated with this mixture
    const components = await ctx.db.query("mixtureComponents")
      .withIndex("by_mixture", q => q.eq("mixtureId", args.id))
      .collect();
    
    for (const component of components) {
      await ctx.db.delete(component._id);
    }
    
    return args.id;
  },
});
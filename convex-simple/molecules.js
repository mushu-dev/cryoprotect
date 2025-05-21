// Simple CRUD functions for molecules
import { mutation, query } from "./_generated/server";
import { v } from "convex/values";

// Create a molecule
export const createMolecule = mutation({
  args: {
    name: v.string(),
    formula: v.optional(v.string()),
  },
  handler: async (ctx, args) => {
    const id = await ctx.db.insert("molecules", {
      name: args.name,
      formula: args.formula,
      createdAt: Date.now(),
    });
    return id;
  },
});

// List all molecules
export const listMolecules = query({
  args: {},
  handler: async (ctx) => {
    return await ctx.db.query("molecules").collect();
  },
});

// Get a molecule by ID
export const getMolecule = query({
  args: { id: v.id("molecules") },
  handler: async (ctx, args) => {
    return await ctx.db.get(args.id);
  },
});
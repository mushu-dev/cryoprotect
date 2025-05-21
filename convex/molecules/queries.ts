import { query } from "../_generated/server";
import { v } from "convex/values";

/**
 * Get a list of all molecules
 */
export const list = query({
  args: {
    limit: v.optional(v.number()),
    status: v.optional(v.string()),
  },
  handler: async (ctx, args) => {
    const limit = args.limit ?? 50;
    let moleculesQuery = ctx.db.query("molecules");
    
    // Filter by status if provided
    if (args.status) {
      moleculesQuery = moleculesQuery.filter(q => q.eq(q.field("status"), args.status));
    }
    
    // Get the most recently updated molecules first
    const molecules = await moleculesQuery
      .order("desc", q => q.field("updatedAt"))
      .take(limit);
    
    return molecules;
  },
});

/**
 * Get a single molecule by ID
 */
export const getById = query({
  args: {
    id: v.id("molecules"),
  },
  handler: async (ctx, args) => {
    return ctx.db.get(args.id);
  },
});

/**
 * Search for molecules by name
 */
export const search = query({
  args: {
    query: v.string(),
    limit: v.optional(v.number()),
  },
  handler: async (ctx, args) => {
    const limit = args.limit ?? 20;
    const searchQuery = args.query.toLowerCase();
    
    // This is a simple implementation - in a production app you'd want
    // to use more sophisticated search with proper indexing
    const molecules = await ctx.db.query("molecules")
      .filter(q => 
        q.eq(q.field("status"), "active")
      )
      .order("desc", q => q.field("updatedAt"))
      .collect();
    
    // Filter client-side by name (simple search)
    const filtered = molecules.filter(
      molecule => molecule.name.toLowerCase().includes(searchQuery)
    );
    
    return filtered.slice(0, limit);
  },
});
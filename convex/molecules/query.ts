/**
 * Query functions for molecules
 */

import { query } from "../_generated/server";
import { v } from "convex/values";
import { Id } from "../_generated/dataModel";
import { MoleculeFilter, MoleculeQueryOptions, MoleculeWithProperties } from "./types";
import { getMoleculeById } from "./helpers";
import { isAuthenticated, pagination } from "../utils/common";

/**
 * Get a molecule by ID
 */
export const getMolecule = query({
  args: {
    id: v.id("molecules"),
    includeProperties: v.optional(v.boolean()),
  },
  handler: async (ctx, args): Promise<MoleculeWithProperties | null> => {
    return await getMoleculeById(ctx.db, args.id, args.includeProperties || false);
  },
});

/**
 * Search for molecules by various criteria
 */
export const searchMolecules = query({
  args: {
    filter: v.optional(
      v.object({
        name: v.optional(v.string()),
        pubchemCid: v.optional(v.string()),
        inchiKey: v.optional(v.string()),
        formula: v.optional(v.string()),
        status: v.optional(v.union(
          v.literal("active"),
          v.literal("deprecated"),
          v.literal("consolidated")
        )),
        consolidated: v.optional(v.boolean()),
        includeDeprecated: v.optional(v.boolean()),
      })
    ),
    options: v.optional(
      v.object({
        limit: v.optional(v.number()),
        cursor: v.optional(v.string()),
        includeProperties: v.optional(v.boolean()),
      })
    ),
  },
  handler: async (ctx, args): Promise<{
    molecules: MoleculeWithProperties[];
    cursor: string | null;
    totalCount: number;
  }> => {
    const filter = args.filter || {};
    const options = args.options || {};
    const limit = options.limit || pagination.defaultPageSize;
    const includeProperties = options.includeProperties || false;
    
    // Start building the query
    let query = ctx.db.query("molecules");
    
    // Apply filters
    if (filter.name) {
      query = query.withSearchIndex("by_name", (q) => 
        q.search("name", filter.name)
      );
    }
    
    if (filter.pubchemCid) {
      query = query.withIndex("by_pubchemCid", (q) => 
        q.eq("pubchemCid", filter.pubchemCid)
      );
    }
    
    if (filter.inchiKey) {
      query = query.withIndex("by_inchiKey", (q) => 
        q.eq("inchiKey", filter.inchiKey)
      );
    }
    
    if (filter.status) {
      query = query.withIndex("by_status", (q) => 
        q.eq("status", filter.status)
      );
    } else if (!filter.includeDeprecated) {
      // By default, only show active molecules
      query = query.filter(q => q.eq(q.field("status"), "active"));
    }
    
    if (filter.consolidated !== undefined) {
      query = query.withIndex("by_consolidated", (q) => 
        q.eq("consolidated", filter.consolidated)
      );
    }
    
    // Apply pagination
    const results = await query
      .paginate({ cursor: options.cursor, numItems: limit });
    
    // Get molecules with properties if requested
    const molecules = await Promise.all(
      results.page.map(async (molecule) => {
        if (includeProperties) {
          return await getMoleculeById(ctx.db, molecule._id, true);
        }
        return molecule as MoleculeWithProperties;
      })
    );
    
    // Get total count (approximate)
    const totalCount = await query.collect();
    
    return {
      molecules,
      cursor: results.continueCursor,
      totalCount: totalCount.length,
    };
  },
});

/**
 * Get molecules by IDs
 */
export const getMoleculesByIds = query({
  args: {
    ids: v.array(v.id("molecules")),
    includeProperties: v.optional(v.boolean()),
  },
  handler: async (ctx, args): Promise<MoleculeWithProperties[]> => {
    const includeProperties = args.includeProperties || false;
    
    return await Promise.all(
      args.ids.map(async (id) => {
        const molecule = await getMoleculeById(ctx.db, id, includeProperties);
        return molecule || null;
      })
    ).then(results => results.filter(Boolean) as MoleculeWithProperties[]);
  },
});

/**
 * Get recent molecules
 */
export const getRecentMolecules = query({
  args: {
    limit: v.optional(v.number()),
    includeProperties: v.optional(v.boolean()),
  },
  handler: async (ctx, args): Promise<MoleculeWithProperties[]> => {
    const limit = args.limit || 10;
    const includeProperties = args.includeProperties || false;
    
    // Query recent active molecules
    const recentMolecules = await ctx.db
      .query("molecules")
      .withIndex("by_status", (q) => q.eq("status", "active"))
      .order("desc")
      .take(limit);
    
    // Get molecules with properties if requested
    return await Promise.all(
      recentMolecules.map(async (molecule) => {
        if (includeProperties) {
          return await getMoleculeById(ctx.db, molecule._id, true);
        }
        return molecule as MoleculeWithProperties;
      })
    );
  },
});
import { mutation } from "../_generated/server";
import { v } from "convex/values";

/**
 * API endpoint for our backend adapter to delete data from Convex tables
 * This function handles deletes from the Supabase-compatible adapter
 */
export const execute = mutation({
  args: {
    table: v.string(),
    filters: v.optional(v.any()),
  },
  handler: async (ctx, args) => {
    try {
      const { table, filters } = args;
      
      // Find records that match the filters
      let query = ctx.db.query(table);
      
      if (filters && Object.keys(filters).length > 0) {
        query = query.filter(q => {
          // Start with a condition that will always be true
          let conditions = q.eq(1, 1);
          
          // Add conditions for each filter
          for (const [field, value] of Object.entries(filters)) {
            conditions = q.and(conditions, q.eq(q.field(field), value));
          }
          
          return conditions;
        });
      }
      
      const records = await query.collect();
      
      // Delete each record and keep track of deleted IDs
      const deleted = [];
      for (const record of records) {
        await ctx.db.delete(record._id);
        deleted.push({
          id: record._id
        });
      }
      
      return { data: deleted, error: null };
    } catch (error) {
      console.error('Error deleting data:', error);
      return { 
        data: [], 
        error: error instanceof Error ? error.message : 'Unknown error' 
      };
    }
  },
});
import { mutation } from "../_generated/server";
import { v } from "convex/values";

/**
 * API endpoint for our backend adapter to update data in Convex tables
 * This function handles updates from the Supabase-compatible adapter
 */
export const execute = mutation({
  args: {
    table: v.string(),
    data: v.any(),
    filters: v.optional(v.any()),
  },
  handler: async (ctx, args) => {
    try {
      const { table, data, filters } = args;
      
      // Add updated timestamp
      const timestamp = Date.now();
      const dataWithTimestamp = {
        ...data,
        updatedAt: timestamp
      };
      
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
      
      // Update each record
      const updated = [];
      for (const record of records) {
        await ctx.db.patch(record._id, dataWithTimestamp);
        
        // Get the updated record
        const updatedRecord = await ctx.db.get(record._id);
        if (updatedRecord) {
          updated.push({
            id: record._id,
            ...updatedRecord
          });
        }
      }
      
      return { data: updated, error: null };
    } catch (error) {
      console.error('Error updating data:', error);
      return { 
        data: [], 
        error: error instanceof Error ? error.message : 'Unknown error' 
      };
    }
  },
});
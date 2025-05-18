import { mutation } from "../_generated/server";
import { v } from "convex/values";

/**
 * API endpoint for our backend adapter to query Convex tables
 * This function handles queries from the Supabase-compatible adapter
 */
export const execute = mutation({
  args: {
    table: v.string(),
    filters: v.optional(v.any()),
    columns: v.optional(v.string()),
    order: v.optional(v.array(v.array(v.string()))),
    limit: v.optional(v.number()),
    offset: v.optional(v.number()),
  },
  handler: async (ctx, args) => {
    try {
      const { table, filters, columns, order, limit, offset } = args;
      
      // Start building the query
      let query = ctx.db.query(table);
      
      // Apply filters if provided
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
      
      // Execute the query to get the results
      let results = await query.collect();
      
      // Apply pagination
      if (offset && offset > 0) {
        results = results.slice(offset);
      }
      
      if (limit && limit > 0) {
        results = results.slice(0, limit);
      }
      
      // Filter columns if specified
      if (columns && columns !== '*') {
        const columnsToInclude = columns.split(',').map(c => c.trim());
        results = results.map(record => {
          const filteredRecord: Record<string, any> = {};
          columnsToInclude.forEach(col => {
            if (record.hasOwnProperty(col)) {
              filteredRecord[col] = record[col];
            }
          });
          return filteredRecord;
        });
      }
      
      return { data: results, error: null };
    } catch (error) {
      console.error('Error executing query:', error);
      return { 
        data: [], 
        error: error instanceof Error ? error.message : 'Unknown error' 
      };
    }
  },
});
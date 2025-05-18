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
      
      // Get table name in proper case (important for Convex)
      let tableName = table;
      
      // Map Supabase table names to Convex table names
      const tableMap: Record<string, string> = {
        'molecules': 'molecules',
        'mixtures': 'mixtures',
        'mixture_components': 'mixtureComponents',
        'molecular_properties': 'molecularProperties',
        'property_types': 'propertyTypes',
        'experiments': 'experiments',
        'experiment_results': 'experimentResults',
        'predictions': 'predictions',
        'scientific_models': 'scientificModels',
        'data_sources': 'dataSources',
        'cross_references': 'crossReferences',
        'synonyms': 'synonyms',
        'users': 'users',
        'teams': 'teams',
        'team_members': 'teamMembers',
        'projects': 'projects'
      };
      
      // Use mapped table name if available
      if (tableMap[tableName]) {
        tableName = tableMap[tableName];
      }
      
      // Find records that match the filters
      let query = ctx.db.query(tableName);
      
      if (filters && Object.keys(filters).length > 0) {
        query = applyFilters(query, filters);
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

/**
 * Helper function to apply filters to a query
 */
function applyFilters(query: any, filters: Record<string, any>) {
  return query.filter(q => {
    // Start with a condition that will always be true
    let conditions = q.eq(1, 1);
    
    for (const [key, value] of Object.entries(filters)) {
      // Handle special operators (e.g., _neq, _gte, _lte)
      if (key.includes('_')) {
        const [field, operator] = key.split('_');
        
        switch (operator) {
          case 'neq':
            conditions = q.and(conditions, q.neq(q.field(field), value));
            break;
          case 'gt':
            conditions = q.and(conditions, q.gt(q.field(field), value));
            break;
          case 'gte':
            conditions = q.and(conditions, q.gte(q.field(field), value));
            break;
          case 'lt':
            conditions = q.and(conditions, q.lt(q.field(field), value));
            break;
          case 'lte':
            conditions = q.and(conditions, q.lte(q.field(field), value));
            break;
          default:
            // Unknown operator, ignore
            break;
        }
      } else {
        // Simple equality
        conditions = q.and(conditions, q.eq(q.field(key), value));
      }
    }
    
    return conditions;
  });
}
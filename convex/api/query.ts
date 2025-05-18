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
      
      // Get table name in proper case (important for Convex)
      let tableName = table;
      
      // Map Supabase table names to Convex table names
      // This is needed because our existing backend expects certain table names
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
      
      // Start building the query
      let query = ctx.db.query(tableName);
      
      // Apply filters
      if (filters && Object.keys(filters).length > 0) {
        query = applyFilters(query, filters);
      }
      
      // Apply ordering
      if (order && order.length > 0) {
        // In our adapter, order is an array of [column, direction] pairs
        const [column, direction] = order[0];
        const dir = direction === 'desc' ? 'desc' : 'asc';
        
        query = query.order(dir, q => q.field(column));
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
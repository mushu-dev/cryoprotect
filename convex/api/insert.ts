import { mutation } from "../_generated/server";
import { v } from "convex/values";

/**
 * API endpoint for our backend adapter to insert data into Convex tables
 * This function handles inserts from the Supabase-compatible adapter
 */
export const execute = mutation({
  args: {
    table: v.string(),
    data: v.any(),
  },
  handler: async (ctx, args) => {
    try {
      const { table, data } = args;
      
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
      
      // Add timestamps
      const timestamp = Date.now();
      
      // Handle array of records or single record
      if (Array.isArray(data)) {
        // Insert multiple records
        const ids = [];
        for (const record of data) {
          const recordWithTimestamps = {
            ...record,
            createdAt: timestamp,
            updatedAt: timestamp
          };
          
          const id = await ctx.db.insert(tableName, recordWithTimestamps);
          ids.push(id);
          
          // Get the inserted record
          const insertedRecord = await ctx.db.get(id);
          if (insertedRecord) {
            ids.push({
              id,
              ...insertedRecord
            });
          }
        }
        
        return { data: ids, error: null };
      } else {
        // Insert a single record
        const recordWithTimestamps = {
          ...data,
          createdAt: timestamp,
          updatedAt: timestamp
        };
        
        const id = await ctx.db.insert(tableName, recordWithTimestamps);
        
        // Get the inserted record
        const insertedRecord = await ctx.db.get(id);
        
        return { 
          data: insertedRecord ? [{ id, ...insertedRecord }] : [{ id }], 
          error: null 
        };
      }
    } catch (error) {
      console.error('Error inserting data:', error);
      return { 
        data: [], 
        error: error instanceof Error ? error.message : 'Unknown error' 
      };
    }
  },
});
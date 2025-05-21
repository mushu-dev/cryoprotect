/**
 * CRUD operations for molecule synonyms
 */

import { mutation, query } from "../_generated/server";
import { v } from "convex/values";
import { Id } from "../_generated/dataModel";
import { 
  MoleculeSynonym, 
  CreateSynonymInput, 
  UpdateSynonymInput,
  SynonymFilter,
  SynonymQueryOptions
} from "./types";
import { 
  validateCreateSynonymInput, 
  validateUpdateSynonymInput,
  validateSynonymExists,
  validateMoleculeExists
} from "./validation";
import { 
  expandSynonymWithMolecule,
  groupSynonymsByType,
  formatSynonym
} from "./helpers";
import { getCurrentUser } from "../auth/users";

/**
 * Create a new molecule synonym
 */
export const createSynonym = mutation({
  args: {
    synonym: v.object({
      moleculeId: v.id("molecules"),
      name: v.string(),
      type: v.optional(v.string()),
      source: v.optional(v.id("dataSources"))
    })
  },
  handler: async (ctx, args) => {
    // Validate input
    validateCreateSynonymInput(args.synonym);
    validateMoleculeExists(args.synonym.moleculeId);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // User must be logged in to create a synonym
    if (!userId) {
      throw new Error("Authentication required to create a synonym");
    }
    
    // Verify molecule exists
    const molecule = await ctx.db.get(args.synonym.moleculeId);
    if (!molecule) {
      throw new Error(`Molecule with ID ${args.synonym.moleculeId} not found`);
    }
    
    // Verify data source exists if provided
    if (args.synonym.source) {
      const dataSource = await ctx.db.get(args.synonym.source);
      if (!dataSource) {
        throw new Error(`Data source with ID ${args.synonym.source} not found`);
      }
    }
    
    // Prepare data for insertion
    const now = Date.now();
    const synonymData = {
      ...args.synonym,
      createdAt: now,
      updatedAt: now
    };
    
    // Check for duplicates
    const existingSynonyms = await ctx.db
      .query("moleculeSynonyms")
      .withIndex("by_molecule", q => q.eq("moleculeId", args.synonym.moleculeId))
      .filter(q => q.eq(q.field("name"), args.synonym.name))
      .collect();
    
    if (existingSynonyms.length > 0) {
      throw new Error(
        `Synonym "${args.synonym.name}" already exists for this molecule`
      );
    }
    
    // Insert synonym
    const synonymId = await ctx.db.insert("moleculeSynonyms", synonymData);
    
    // Create audit log entry
    await ctx.db.insert("scientificDataAudit", {
      table: "moleculeSynonyms",
      documentId: synonymId,
      operation: "create",
      userId,
      newValue: synonymData,
      timestamp: now,
    });
    
    return synonymId;
  }
});

/**
 * Get a synonym by ID
 */
export const getSynonym = query({
  args: {
    synonymId: v.id("moleculeSynonyms"),
    includeMolecule: v.optional(v.boolean())
  },
  handler: async (ctx, args) => {
    // Validate synonym exists
    validateSynonymExists(args.synonymId);
    
    // Get the synonym
    const synonym = await ctx.db.get(args.synonymId);
    if (!synonym) {
      return null;
    }
    
    // Include molecule data if requested
    if (args.includeMolecule) {
      return expandSynonymWithMolecule(ctx, synonym);
    }
    
    return synonym;
  }
});

/**
 * Update an existing synonym
 */
export const updateSynonym = mutation({
  args: {
    synonymId: v.id("moleculeSynonyms"),
    update: v.object({
      name: v.optional(v.string()),
      type: v.optional(v.string()),
      source: v.optional(v.id("dataSources"))
    })
  },
  handler: async (ctx, args) => {
    // Validate input
    validateUpdateSynonymInput(args.update);
    validateSynonymExists(args.synonymId);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // User must be logged in to update a synonym
    if (!userId) {
      throw new Error("Authentication required to update a synonym");
    }
    
    // Get existing synonym
    const existingSynonym = await ctx.db.get(args.synonymId);
    if (!existingSynonym) {
      throw new Error(`Synonym with ID ${args.synonymId} not found`);
    }
    
    // Verify data source exists if provided
    if (args.update.source) {
      const dataSource = await ctx.db.get(args.update.source);
      if (!dataSource) {
        throw new Error(`Data source with ID ${args.update.source} not found`);
      }
    }
    
    // Prepare update data
    const now = Date.now();
    const updateData = {
      ...args.update,
      updatedAt: now
    };
    
    // Check for duplicates if name is changing
    if (args.update.name && args.update.name !== existingSynonym.name) {
      const existingSynonyms = await ctx.db
        .query("moleculeSynonyms")
        .withIndex("by_molecule", q => q.eq("moleculeId", existingSynonym.moleculeId))
        .filter(q => 
          q.and(
            q.eq(q.field("name"), args.update.name!),
            q.neq(q.field("_id"), args.synonymId)
          )
        )
        .collect();
      
      if (existingSynonyms.length > 0) {
        throw new Error(
          `Synonym "${args.update.name}" already exists for this molecule`
        );
      }
    }
    
    // Update synonym
    await ctx.db.patch(args.synonymId, updateData);
    
    // Create audit log entry
    await ctx.db.insert("scientificDataAudit", {
      table: "moleculeSynonyms",
      documentId: args.synonymId,
      operation: "update",
      userId,
      previousValue: existingSynonym,
      newValue: { ...existingSynonym, ...updateData },
      timestamp: now
    });
    
    return args.synonymId;
  }
});

/**
 * Delete a synonym
 */
export const deleteSynonym = mutation({
  args: {
    synonymId: v.id("moleculeSynonyms")
  },
  handler: async (ctx, args) => {
    // Validate synonym exists
    validateSynonymExists(args.synonymId);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // User must be logged in to delete a synonym
    if (!userId) {
      throw new Error("Authentication required to delete a synonym");
    }
    
    // Get existing synonym
    const existingSynonym = await ctx.db.get(args.synonymId);
    if (!existingSynonym) {
      throw new Error(`Synonym with ID ${args.synonymId} not found`);
    }
    
    // Create audit log entry
    const now = Date.now();
    await ctx.db.insert("scientificDataAudit", {
      table: "moleculeSynonyms",
      documentId: args.synonymId,
      operation: "delete",
      userId,
      previousValue: existingSynonym,
      timestamp: now
    });
    
    // Delete the synonym
    await ctx.db.delete(args.synonymId);
    
    return true;
  }
});

/**
 * List synonyms with optional filtering
 */
export const listSynonyms = query({
  args: {
    filter: v.optional(v.object({
      moleculeId: v.optional(v.id("molecules")),
      name: v.optional(v.string()),
      type: v.optional(v.string()),
      source: v.optional(v.id("dataSources")),
      includeMolecule: v.optional(v.boolean())
    })),
    options: v.optional(v.object({
      limit: v.optional(v.number()),
      cursor: v.optional(v.string()),
      groupByType: v.optional(v.boolean())
    }))
  },
  handler: async (ctx, args) => {
    // Set up the query
    let query = ctx.db.query("moleculeSynonyms");
    
    // Apply filters
    if (args.filter) {
      if (args.filter.moleculeId) {
        query = query.withIndex("by_molecule", q => 
          q.eq("moleculeId", args.filter!.moleculeId!)
        );
      }
      
      if (args.filter.name) {
        query = query.withIndex("by_name", q => 
          q.eq("name", args.filter!.name!)
        );
      }
      
      if (args.filter.type) {
        query = query.filter(q => 
          q.eq(q.field("type"), args.filter!.type!)
        );
      }
      
      if (args.filter.source) {
        query = query.filter(q => 
          q.eq(q.field("source"), args.filter!.source!)
        );
      }
    }
    
    // Apply pagination
    if (args.options?.cursor) {
      query = query.withCursor(args.options.cursor);
    }
    
    if (args.options?.limit) {
      query = query.take(args.options.limit);
    } else {
      query = query.take(50); // Default limit
    }
    
    // Execute query
    const synonyms = await query.collect();
    
    // Group by type if requested
    if (args.options?.groupByType) {
      return groupSynonymsByType(synonyms);
    }
    
    // Include molecule data if requested
    if (args.filter?.includeMolecule) {
      const synonymsWithMolecules = [];
      
      for (const synonym of synonyms) {
        synonymsWithMolecules.push(await expandSynonymWithMolecule(ctx, synonym));
      }
      
      return synonymsWithMolecules;
    }
    
    return synonyms;
  }
});

/**
 * Get all synonyms for a molecule
 */
export const getMoleculeSynonyms = query({
  args: {
    moleculeId: v.id("molecules"),
    groupByType: v.optional(v.boolean())
  },
  handler: async (ctx, args) => {
    // Validate molecule exists
    validateMoleculeExists(args.moleculeId);
    
    // Query for synonyms
    const synonyms = await ctx.db
      .query("moleculeSynonyms")
      .withIndex("by_molecule", q => q.eq("moleculeId", args.moleculeId))
      .collect();
    
    // Group by type if requested
    if (args.groupByType) {
      return groupSynonymsByType(synonyms);
    }
    
    return synonyms;
  }
});

/**
 * Search for molecules by synonym
 */
export const searchBySynonym = query({
  args: {
    query: v.string(),
    limit: v.optional(v.number())
  },
  handler: async (ctx, args) => {
    // Perform a search by synonym name
    let query = ctx.db.query("moleculeSynonyms")
      .filter(q => 
        q.includes("name", args.query)
      );
    
    // Apply limit
    const limit = args.limit || 20;
    query = query.take(limit);
    
    // Execute query
    const synonyms = await query.collect();
    
    // Expand with molecule details
    const results = [];
    
    for (const synonym of synonyms) {
      results.push(await expandSynonymWithMolecule(ctx, synonym));
    }
    
    return results;
  }
});

/**
 * Batch add multiple synonyms for a molecule
 */
export const batchAddSynonyms = mutation({
  args: {
    moleculeId: v.id("molecules"),
    synonyms: v.array(
      v.object({
        name: v.string(),
        type: v.optional(v.string()),
        source: v.optional(v.id("dataSources"))
      })
    )
  },
  handler: async (ctx, args) => {
    // Validate molecule exists
    validateMoleculeExists(args.moleculeId);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // User must be logged in to add synonyms
    if (!userId) {
      throw new Error("Authentication required to add synonyms");
    }
    
    // Verify molecule exists
    const molecule = await ctx.db.get(args.moleculeId);
    if (!molecule) {
      throw new Error(`Molecule with ID ${args.moleculeId} not found`);
    }
    
    // Record time for all operations
    const now = Date.now();
    
    // Get existing synonyms to avoid duplicates
    const existingSynonyms = await ctx.db
      .query("moleculeSynonyms")
      .withIndex("by_molecule", q => q.eq("moleculeId", args.moleculeId))
      .collect();
    
    // Create a set for checking duplicates
    const existingSynonymSet = new Set(
      existingSynonyms.map(synonym => synonym.name.toLowerCase())
    );
    
    // Validate data sources once
    const sourcesToCheck = new Set<string>();
    for (const synonym of args.synonyms) {
      if (synonym.source) {
        sourcesToCheck.add(synonym.source.toString());
      }
    }
    
    for (const sourceId of sourcesToCheck) {
      const source = await ctx.db.get(sourceId as Id<"dataSources">);
      if (!source) {
        throw new Error(`Data source with ID ${sourceId} not found`);
      }
    }
    
    // Validate and insert each synonym
    const synonymIds: Id<"moleculeSynonyms">[] = [];
    
    for (const synonym of args.synonyms) {
      // Validate synonym input
      validateCreateSynonymInput({
        ...synonym,
        moleculeId: args.moleculeId
      });
      
      // Check for duplicates
      const synKey = synonym.name.toLowerCase();
      if (existingSynonymSet.has(synKey)) {
        continue; // Skip duplicates
      }
      
      // Prepare synonym data
      const synonymData = {
        ...synonym,
        moleculeId: args.moleculeId,
        createdAt: now,
        updatedAt: now
      };
      
      // Insert synonym
      const synonymId = await ctx.db.insert("moleculeSynonyms", synonymData);
      
      // Create audit log entry
      await ctx.db.insert("scientificDataAudit", {
        table: "moleculeSynonyms",
        documentId: synonymId,
        operation: "create",
        userId,
        newValue: synonymData,
        timestamp: now
      });
      
      synonymIds.push(synonymId);
      existingSynonymSet.add(synKey); // Add to set to prevent duplicates in batch
    }
    
    return synonymIds;
  }
});
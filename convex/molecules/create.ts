/**
 * Functions for creating new molecules
 */

import { mutation } from "../_generated/server";
import { v } from "convex/values";
import { CreateMoleculeInput, Molecule } from "./types";
import { validateCreateMoleculeInput } from "./validation";
import { findDuplicateMolecules } from "./helpers";
import { createAuditLog, getCurrentTimestamp, validateAuthenticated } from "../utils/common";

/**
 * Create a new molecule
 */
export const createMolecule = mutation({
  args: {
    input: v.object({
      name: v.string(),
      pubchemCid: v.optional(v.string()),
      canonicalSmiles: v.optional(v.string()),
      inchiKey: v.optional(v.string()),
      formula: v.optional(v.string()),
      dataSource: v.optional(v.id("dataSources")),
      sourceId: v.optional(v.string()),
    }),
  },
  handler: async (ctx, args): Promise<Molecule> => {
    // Validate user authentication
    validateAuthenticated(ctx);
    
    // Validate input
    const validationResult = validateCreateMoleculeInput(args.input);
    if (!validationResult.valid) {
      throw new Error(validationResult.error);
    }
    
    // Check for duplicates
    const duplicates = await findDuplicateMolecules(ctx.db, args.input);
    if (duplicates.length > 0) {
      throw new Error(`Duplicate molecule found with identifier: ${duplicates[0].name}`);
    }
    
    // Create molecule document
    const now = getCurrentTimestamp();
    const moleculeId = await ctx.db.insert("molecules", {
      name: args.input.name,
      pubchemCid: args.input.pubchemCid,
      canonicalSmiles: args.input.canonicalSmiles,
      inchiKey: args.input.inchiKey,
      formula: args.input.formula,
      status: "active",
      consolidated: false,
      createdAt: now,
      updatedAt: now,
      dataSource: args.input.dataSource,
      sourceId: args.input.sourceId,
    });
    
    // Fetch the inserted molecule
    const molecule = await ctx.db.get(moleculeId) as Molecule;
    
    // Create audit log
    await createAuditLog(
      ctx,
      "molecules",
      moleculeId.toString(),
      "create",
      undefined,
      molecule,
      "Initial creation"
    );
    
    return molecule;
  },
});

/**
 * Create multiple molecules in a batch
 */
export const createMoleculesBatch = mutation({
  args: {
    molecules: v.array(
      v.object({
        name: v.string(),
        pubchemCid: v.optional(v.string()),
        canonicalSmiles: v.optional(v.string()),
        inchiKey: v.optional(v.string()),
        formula: v.optional(v.string()),
        dataSource: v.optional(v.id("dataSources")),
        sourceId: v.optional(v.string()),
      })
    ),
    skipDuplicates: v.optional(v.boolean()),
  },
  handler: async (ctx, args): Promise<{
    created: number;
    skipped: number;
    failed: number;
    errors: Array<{ input: CreateMoleculeInput; error: string }>;
  }> => {
    // Validate user authentication
    validateAuthenticated(ctx);
    
    const results = {
      created: 0,
      skipped: 0,
      failed: 0,
      errors: [] as Array<{ input: CreateMoleculeInput; error: string }>,
    };
    
    const now = getCurrentTimestamp();
    
    // Process each molecule in the batch
    for (const input of args.molecules) {
      try {
        // Validate input
        const validationResult = validateCreateMoleculeInput(input);
        if (!validationResult.valid) {
          throw new Error(validationResult.error);
        }
        
        // Check for duplicates
        const duplicates = await findDuplicateMolecules(ctx.db, input);
        if (duplicates.length > 0) {
          if (args.skipDuplicates) {
            results.skipped++;
            continue;
          } else {
            throw new Error(`Duplicate molecule found with identifier: ${duplicates[0].name}`);
          }
        }
        
        // Create molecule document
        const moleculeId = await ctx.db.insert("molecules", {
          name: input.name,
          pubchemCid: input.pubchemCid,
          canonicalSmiles: input.canonicalSmiles,
          inchiKey: input.inchiKey,
          formula: input.formula,
          status: "active",
          consolidated: false,
          createdAt: now,
          updatedAt: now,
          dataSource: input.dataSource,
          sourceId: input.sourceId,
        });
        
        // Create audit log
        await createAuditLog(
          ctx,
          "molecules",
          moleculeId.toString(),
          "create",
          undefined,
          await ctx.db.get(moleculeId),
          "Batch creation"
        );
        
        results.created++;
      } catch (error) {
        results.failed++;
        results.errors.push({
          input,
          error: error instanceof Error ? error.message : 'Unknown error',
        });
      }
    }
    
    return results;
  },
});
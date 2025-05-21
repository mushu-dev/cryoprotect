/**
 * Protocol Experiment Links API for CryoProtect
 * 
 * This module provides functions for creating and managing links between
 * protocols and experiments for experimental data tracking.
 */

import { mutation, query } from "./_generated/server";
import { v } from "convex/values";
import { Id } from "./_generated/dataModel";
import { getCurrentUser } from "./auth/users";
import { getProtocol } from "./experiments/protocols";

/**
 * Create a new link between a protocol and experiment
 */
export const createLink = mutation({
  args: {
    protocolId: v.id("protocols"),
    experimentId: v.string(),
    metadata: v.optional(v.object({})),
  },
  handler: async (ctx, args) => {
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;

    // Verify protocol exists
    const protocol = await getProtocol.internal(ctx, { protocolId: args.protocolId });
    if (!protocol) {
      throw new Error(`Protocol with ID ${args.protocolId} not found`);
    }

    // Check if link already exists
    const existingLink = await ctx.db
      .query("protocol_experiment_links")
      .withIndex("by_protocol_experiment", (q) => 
        q.eq("protocolId", args.protocolId).eq("experimentId", args.experimentId)
      )
      .first();

    if (existingLink) {
      return existingLink._id;
    }

    // Create link
    const linkId = await ctx.db.insert("protocol_experiment_links", {
      protocolId: args.protocolId,
      experimentId: args.experimentId,
      createdBy: userId,
      createdAt: Date.now(),
      metadata: args.metadata || {},
    });

    // Create audit log entry
    await ctx.db.insert("scientificDataAudit", {
      table: "protocol_experiment_links",
      documentId: linkId,
      operation: "create",
      userId,
      timestamp: Date.now(),
    });

    return linkId;
  },
});

/**
 * Get links for a protocol or experiment
 */
export const getLinks = query({
  args: {
    protocolId: v.optional(v.id("protocols")),
    experimentId: v.optional(v.string()),
  },
  handler: async (ctx, args) => {
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;

    // Make sure at least one argument is provided
    if (!args.protocolId && !args.experimentId) {
      throw new Error("Either protocolId or experimentId must be provided");
    }

    // Query based on provided arguments
    let query = ctx.db.query("protocol_experiment_links");

    if (args.protocolId && args.experimentId) {
      // Query for a specific protocol-experiment pair
      query = query.withIndex("by_protocol_experiment", (q) => 
        q.eq("protocolId", args.protocolId).eq("experimentId", args.experimentId)
      );
    } else if (args.protocolId) {
      // Query for a specific protocol
      query = query.withIndex("by_protocol", (q) => 
        q.eq("protocolId", args.protocolId)
      );
    } else if (args.experimentId) {
      // Query for a specific experiment
      query = query.withIndex("by_experiment", (q) => 
        q.eq("experimentId", args.experimentId)
      );
    }

    // Execute query
    const links = await query.collect();

    // If not logged in, only return public links
    if (!userId) {
      // In this case, we need to check if the linked protocols are public
      const filteredLinks = [];
      for (const link of links) {
        const protocol = await ctx.db.get(link.protocolId);
        if (protocol && protocol.public) {
          filteredLinks.push(link);
        }
      }
      return filteredLinks;
    }

    return links;
  },
});

/**
 * Delete a link between a protocol and experiment
 */
export const deleteLink = mutation({
  args: {
    linkId: v.id("protocol_experiment_links"),
  },
  handler: async (ctx, args) => {
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;

    // Get link
    const link = await ctx.db.get(args.linkId);
    if (!link) {
      throw new Error(`Link with ID ${args.linkId} not found`);
    }

    // Get protocol to check permissions
    const protocol = await ctx.db.get(link.protocolId);
    if (!protocol) {
      throw new Error(`Protocol with ID ${link.protocolId} not found`);
    }

    // Check permissions - only allow delete if user created the protocol or link
    if (
      !userId || 
      (protocol.createdBy && !userId.equals(protocol.createdBy) && 
       link.createdBy && !userId.equals(link.createdBy))
    ) {
      throw new Error("You do not have permission to delete this link");
    }

    // Create audit log entry
    await ctx.db.insert("scientificDataAudit", {
      table: "protocol_experiment_links",
      documentId: args.linkId,
      operation: "delete",
      userId,
      previousValue: link,
      timestamp: Date.now(),
    });

    // Delete link
    await ctx.db.delete(args.linkId);

    return true;
  },
});

/**
 * Get a single link by ID
 */
export const getLinkById = query({
  args: {
    linkId: v.id("protocol_experiment_links"),
  },
  handler: async (ctx, args) => {
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;

    // Get link
    const link = await ctx.db.get(args.linkId);
    if (!link) {
      return null;
    }

    // Get protocol to check permissions
    const protocol = await ctx.db.get(link.protocolId);
    if (!protocol) {
      return null;
    }

    // Check permissions - only allow view if protocol is public or user created it
    if (
      !protocol.public && 
      (!userId || (protocol.createdBy && !userId.equals(protocol.createdBy)))
    ) {
      throw new Error("You do not have permission to view this link");
    }

    return link;
  },
});

/**
 * Record results for a protocol step in an experiment
 */
export const recordStepResults = mutation({
  args: {
    linkId: v.id("protocol_experiment_links"),
    stepId: v.string(),
    results: v.any(),
  },
  handler: async (ctx, args) => {
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;

    // Get link
    const link = await ctx.db.get(args.linkId);
    if (!link) {
      throw new Error(`Link with ID ${args.linkId} not found`);
    }

    // Check if result already exists
    const existingResult = await ctx.db
      .query("protocol_step_results")
      .withIndex("by_link_step", (q) => 
        q.eq("linkId", args.linkId).eq("stepId", args.stepId)
      )
      .first();

    if (existingResult) {
      // Update existing result
      await ctx.db.patch(existingResult._id, {
        results: args.results,
        updatedBy: userId,
        updatedAt: Date.now(),
      });

      // Create audit log entry
      await ctx.db.insert("scientificDataAudit", {
        table: "protocol_step_results",
        documentId: existingResult._id,
        operation: "update",
        userId,
        previousValue: existingResult,
        timestamp: Date.now(),
      });

      return existingResult._id;
    } else {
      // Create new result
      const resultId = await ctx.db.insert("protocol_step_results", {
        linkId: args.linkId,
        stepId: args.stepId,
        protocolId: link.protocolId,
        experimentId: link.experimentId,
        results: args.results,
        createdBy: userId,
        createdAt: Date.now(),
        updatedAt: Date.now(),
      });

      // Create audit log entry
      await ctx.db.insert("scientificDataAudit", {
        table: "protocol_step_results",
        documentId: resultId,
        operation: "create",
        userId,
        timestamp: Date.now(),
      });

      return resultId;
    }
  },
});

/**
 * Get results for a protocol step in an experiment
 */
export const getStepResults = query({
  args: {
    linkId: v.optional(v.id("protocol_experiment_links")),
    protocolId: v.optional(v.id("protocols")),
    experimentId: v.optional(v.string()),
    stepId: v.string(),
  },
  handler: async (ctx, args) => {
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;

    let link: any = null;

    // If linkId is provided, use it directly
    if (args.linkId) {
      link = await ctx.db.get(args.linkId);
      if (!link) {
        throw new Error(`Link with ID ${args.linkId} not found`);
      }
    } 
    // If protocolId and experimentId are provided, find the link
    else if (args.protocolId && args.experimentId) {
      link = await ctx.db
        .query("protocol_experiment_links")
        .withIndex("by_protocol_experiment", (q) => 
          q.eq("protocolId", args.protocolId).eq("experimentId", args.experimentId)
        )
        .first();

      if (!link) {
        throw new Error(`No link found for protocol ${args.protocolId} and experiment ${args.experimentId}`);
      }
    } else {
      throw new Error("Either linkId or both protocolId and experimentId must be provided");
    }

    // Check permissions
    const protocol = await ctx.db.get(link.protocolId);
    if (!protocol) {
      throw new Error(`Protocol with ID ${link.protocolId} not found`);
    }

    if (
      !protocol.public && 
      (!userId || (protocol.createdBy && !userId.equals(protocol.createdBy)))
    ) {
      throw new Error("You do not have permission to view results for this protocol");
    }

    // Get result
    const result = await ctx.db
      .query("protocol_step_results")
      .withIndex("by_link_step", (q) => 
        q.eq("linkId", link._id).eq("stepId", args.stepId)
      )
      .first();

    return result;
  },
});

/**
 * Get all results for an experiment
 */
export const getAllResultsForExperiment = query({
  args: {
    experimentId: v.string(),
  },
  handler: async (ctx, args) => {
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;

    // Get all links for this experiment
    const links = await ctx.db
      .query("protocol_experiment_links")
      .withIndex("by_experiment", (q) => q.eq("experimentId", args.experimentId))
      .collect();

    if (links.length === 0) {
      return [];
    }

    // For each link, check permissions and get results
    const allResults = [];
    for (const link of links) {
      // Check permissions
      const protocol = await ctx.db.get(link.protocolId);
      if (!protocol) {
        continue;
      }

      if (
        !protocol.public && 
        (!userId || (protocol.createdBy && !userId.equals(protocol.createdBy)))
      ) {
        continue;
      }

      // Get results for this link
      const results = await ctx.db
        .query("protocol_step_results")
        .withIndex("by_link", (q) => q.eq("linkId", link._id))
        .collect();

      // Add protocol info to each result
      const resultsWithProtocol = results.map(result => ({
        ...result,
        protocol: {
          id: protocol._id,
          name: protocol.name,
          version: protocol.version,
        }
      }));

      allResults.push(...resultsWithProtocol);
    }

    return allResults;
  },
});
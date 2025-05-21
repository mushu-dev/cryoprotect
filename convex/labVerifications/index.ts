/**
 * Lab Verification System Implementation
 * 
 * This module provides the API for creating, updating, and reviewing
 * lab verifications to ensure experimental reproducibility.
 */

import { mutation, query } from "../_generated/server";
import { v } from "convex/values";
import { Id } from "../_generated/dataModel";
import { 
  VerificationStatus, 
  CreateVerificationParams, 
  UpdateVerificationParams,
  VerifyExperimentParams,
  LabVerificationFilter,
  LabVerificationQueryOptions
} from "./types";
import {
  validateCreateVerification,
  validateUpdateVerification,
  validateVerifyExperiment,
  validateExperimentExists,
  validateVerificationExists,
  validateCanVerify,
  calculateOverallRating,
  canUpdateVerificationStatus
} from "./validation";
import { getCurrentUser } from "../auth/users";

/**
 * Request verification for an experiment
 */
export const requestVerification = mutation({
  args: {
    verification: v.object({
      experimentId: v.id("enhancedExperiments"),
      equipmentUsed: v.array(v.string()),
      requestNotes: v.optional(v.string()),
      methodologyDescription: v.optional(v.string()),
      controlProcedures: v.optional(v.string()),
      evidenceUrls: v.optional(v.array(v.string())),
      reviewedProtocolSteps: v.optional(v.array(v.string()))
    })
  },
  handler: async (ctx, args) => {
    // Validate input
    validateCreateVerification(args.verification);
    
    // Validate experiment exists
    await validateExperimentExists(ctx.db, args.verification.experimentId);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    if (!userId) {
      throw new Error("Authentication required to request verification");
    }
    
    // Check if verification already exists for this experiment
    const existingVerification = await ctx.db
      .query("labVerifications")
      .withIndex("by_experiment", q => 
        q.eq("experimentId", args.verification.experimentId)
      )
      .first();
    
    if (existingVerification) {
      throw new Error("Verification already exists for this experiment");
    }
    
    // Prepare verification data
    const now = Date.now();
    const verificationData = {
      experimentId: args.verification.experimentId,
      verificationStatus: "pending" as VerificationStatus,
      requesterId: userId,
      equipmentUsed: args.verification.equipmentUsed,
      methodologyDescription: args.verification.methodologyDescription,
      controlProcedures: args.verification.controlProcedures,
      requestNotes: args.verification.requestNotes,
      evidenceUrls: args.verification.evidenceUrls,
      reviewedProtocolSteps: args.verification.reviewedProtocolSteps,
      requestDate: now,
      createdAt: now,
      updatedAt: now
    };
    
    // Insert verification
    const verificationId = await ctx.db.insert("labVerifications", verificationData);
    
    // Update experiment status to reflect verification in progress
    const experiment = await ctx.db.get(args.verification.experimentId);
    if (experiment) {
      await ctx.db.patch(args.verification.experimentId, {
        verificationStatus: "pending_verification",
        updatedAt: now
      });
    }
    
    return verificationId;
  }
});

/**
 * Verify an experiment (approve, reject, or request revisions)
 */
export const verifyExperiment = mutation({
  args: {
    verificationId: v.id("labVerifications"),
    status: v.string(),
    verifierNotes: v.optional(v.string()),
    reproducibilityRating: v.optional(v.number()),
    qualityRating: v.optional(v.number()),
    documentationRating: v.optional(v.number()),
    evidenceUrls: v.optional(v.array(v.string()))
  },
  handler: async (ctx, args) => {
    // Validate input
    validateVerifyExperiment({
      verificationId: args.verificationId,
      verificationStatus: args.status as VerificationStatus,
      verifierNotes: args.verifierNotes,
      reproducibilityRating: args.reproducibilityRating,
      qualityRating: args.qualityRating,
      documentationRating: args.documentationRating,
      evidenceUrls: args.evidenceUrls
    });
    
    // Validate verification exists
    await validateVerificationExists(ctx.db, args.verificationId);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Validate user can verify
    await validateCanVerify(ctx.db, userId);
    
    // Get existing verification
    const verification = await ctx.db.get(args.verificationId);
    if (!verification) {
      throw new Error(`Verification with ID ${args.verificationId} not found`);
    }
    
    // Check if status transition is allowed
    if (!canUpdateVerificationStatus(
      verification.verificationStatus as VerificationStatus,
      args.status as VerificationStatus
    )) {
      throw new Error(`Cannot transition from ${verification.verificationStatus} to ${args.status}`);
    }
    
    // Calculate overall rating if all component ratings are provided
    const overallRating = calculateOverallRating(
      args.reproducibilityRating,
      args.qualityRating,
      args.documentationRating
    );
    
    // Prepare update data
    const now = Date.now();
    const updateData = {
      verificationStatus: args.status,
      verifierId: userId,
      verifierNotes: args.verifierNotes,
      verificationDate: now,
      reproducibilityRating: args.reproducibilityRating,
      qualityRating: args.qualityRating,
      documentationRating: args.documentationRating,
      overallRating,
      updatedAt: now
    };
    
    // Add evidence URLs if provided
    if (args.evidenceUrls) {
      // Merge with existing URLs if any
      const existingUrls = verification.evidenceUrls || [];
      updateData.evidenceUrls = [...new Set([...existingUrls, ...args.evidenceUrls])];
    }
    
    // Update verification
    await ctx.db.patch(args.verificationId, updateData);
    
    // Update experiment status to reflect verification result
    if (verification.experimentId) {
      let experimentStatus: string;
      
      switch (args.status) {
        case "verified":
          experimentStatus = "verified";
          break;
        case "rejected":
          experimentStatus = "verification_failed";
          break;
        case "needs_revision":
          experimentStatus = "verification_revision_needed";
          break;
        default:
          experimentStatus = "pending_verification";
      }
      
      await ctx.db.patch(verification.experimentId, {
        verificationStatus: experimentStatus,
        updatedAt: now
      });
    }
    
    return args.verificationId;
  }
});

/**
 * Update an existing verification
 */
export const updateVerification = mutation({
  args: {
    verificationId: v.id("labVerifications"),
    update: v.object({
      verificationStatus: v.optional(v.string()),
      equipmentUsed: v.optional(v.array(v.string())),
      requestNotes: v.optional(v.string()),
      verifierNotes: v.optional(v.string()),
      methodologyDescription: v.optional(v.string()),
      controlProcedures: v.optional(v.string()),
      evidenceUrls: v.optional(v.array(v.string())),
      reviewedProtocolSteps: v.optional(v.array(v.string())),
      reproducibilityRating: v.optional(v.number()),
      qualityRating: v.optional(v.number()),
      documentationRating: v.optional(v.number())
    })
  },
  handler: async (ctx, args) => {
    // Validate input
    validateUpdateVerification({
      verificationId: args.verificationId,
      ...args.update
    });
    
    // Validate verification exists
    await validateVerificationExists(ctx.db, args.verificationId);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    if (!userId) {
      throw new Error("Authentication required to update verification");
    }
    
    // Get existing verification
    const verification = await ctx.db.get(args.verificationId);
    if (!verification) {
      throw new Error(`Verification with ID ${args.verificationId} not found`);
    }
    
    // Validate permissions
    const isRequester = verification.requesterId && verification.requesterId.equals(userId);
    const isVerifier = verification.verifierId && verification.verifierId.equals(userId);
    
    // Only requesters can update certain fields, and only verifiers can update other fields
    const requesterFields = ["equipmentUsed", "requestNotes", "methodologyDescription", 
                             "controlProcedures", "reviewedProtocolSteps"];
    const verifierFields = ["verifierNotes", "reproducibilityRating", "qualityRating", 
                            "documentationRating", "verificationStatus"];
    
    const requesterUpdating = Object.keys(args.update).some(field => requesterFields.includes(field));
    const verifierUpdating = Object.keys(args.update).some(field => verifierFields.includes(field));
    
    if (requesterUpdating && !isRequester) {
      throw new Error("Only the verification requester can update these fields");
    }
    
    if (verifierUpdating && !isVerifier) {
      throw new Error("Only the verifier can update these fields");
    }
    
    // Check status transition if status is being updated
    if (args.update.verificationStatus && 
        !canUpdateVerificationStatus(
          verification.verificationStatus as VerificationStatus,
          args.update.verificationStatus as VerificationStatus
        )) {
      throw new Error(`Cannot transition from ${verification.verificationStatus} to ${args.update.verificationStatus}`);
    }
    
    // Calculate overall rating if all component ratings are provided
    let overallRating = verification.overallRating;
    
    // If any rating is updated, recalculate overall rating
    if (args.update.reproducibilityRating !== undefined || 
        args.update.qualityRating !== undefined || 
        args.update.documentationRating !== undefined) {
      
      const reproducibilityRating = args.update.reproducibilityRating ?? verification.reproducibilityRating;
      const qualityRating = args.update.qualityRating ?? verification.qualityRating;
      const documentationRating = args.update.documentationRating ?? verification.documentationRating;
      
      overallRating = calculateOverallRating(
        reproducibilityRating,
        qualityRating,
        documentationRating
      );
    }
    
    // Prepare update data
    const now = Date.now();
    const updateData = {
      ...args.update,
      overallRating,
      updatedAt: now
    };
    
    // Handle evidence URLs if provided - merge with existing
    if (args.update.evidenceUrls) {
      const existingUrls = verification.evidenceUrls || [];
      updateData.evidenceUrls = [...new Set([...existingUrls, ...args.update.evidenceUrls])];
    }
    
    // Update verification
    await ctx.db.patch(args.verificationId, updateData);
    
    // If status changed, update experiment status
    if (args.update.verificationStatus && verification.experimentId) {
      let experimentStatus: string;
      
      switch (args.update.verificationStatus) {
        case "verified":
          experimentStatus = "verified";
          break;
        case "rejected":
          experimentStatus = "verification_failed";
          break;
        case "needs_revision":
          experimentStatus = "verification_revision_needed";
          break;
        default:
          experimentStatus = "pending_verification";
      }
      
      await ctx.db.patch(verification.experimentId, {
        verificationStatus: experimentStatus,
        updatedAt: now
      });
    }
    
    return args.verificationId;
  }
});

/**
 * Get a verification by ID
 */
export const getVerification = query({
  args: {
    verificationId: v.id("labVerifications"),
    includeExperiment: v.optional(v.boolean())
  },
  handler: async (ctx, args) => {
    // Get verification
    const verification = await ctx.db.get(args.verificationId);
    if (!verification) {
      return null;
    }
    
    // Include experiment details if requested
    if (args.includeExperiment && verification.experimentId) {
      const experiment = await ctx.db.get(verification.experimentId);
      if (experiment) {
        return {
          ...verification,
          experiment: {
            _id: experiment._id,
            name: experiment.name,
            description: experiment.description,
            status: experiment.status,
            protocolId: experiment.protocolId,
            conductedBy: experiment.conductedBy,
            date: experiment.date
          }
        };
      }
    }
    
    return verification;
  }
});

/**
 * Get a verification by experiment ID
 */
export const getVerificationByExperiment = query({
  args: {
    experimentId: v.id("enhancedExperiments"),
    includeExperiment: v.optional(v.boolean())
  },
  handler: async (ctx, args) => {
    // Get verification for this experiment
    const verification = await ctx.db
      .query("labVerifications")
      .withIndex("by_experiment", q => 
        q.eq("experimentId", args.experimentId)
      )
      .first();
    
    if (!verification) {
      return null;
    }
    
    // Include experiment details if requested
    if (args.includeExperiment) {
      const experiment = await ctx.db.get(args.experimentId);
      if (experiment) {
        return {
          ...verification,
          experiment: {
            _id: experiment._id,
            name: experiment.name,
            description: experiment.description,
            status: experiment.status,
            protocolId: experiment.protocolId,
            conductedBy: experiment.conductedBy,
            date: experiment.date
          }
        };
      }
    }
    
    return verification;
  }
});

/**
 * List verifications with optional filtering
 */
export const listVerifications = query({
  args: {
    filter: v.optional(v.object({
      status: v.optional(v.string()),
      verifierId: v.optional(v.id("users")),
      requesterId: v.optional(v.id("users")),
      experimentId: v.optional(v.id("enhancedExperiments")),
      dateRange: v.optional(v.object({
        start: v.optional(v.number()),
        end: v.optional(v.number())
      }))
    })),
    options: v.optional(v.object({
      limit: v.optional(v.number()),
      cursor: v.optional(v.string()),
      sortBy: v.optional(v.string()),
      sortDirection: v.optional(v.string()),
      includeExperiment: v.optional(v.boolean())
    }))
  },
  handler: async (ctx, args) => {
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Setup query
    let query = ctx.db.query("labVerifications");
    
    // Apply filters
    if (args.filter) {
      // Filter by status
      if (args.filter.status && isValidStatus(args.filter.status)) {
        query = query.filter(q => 
          q.eq(q.field("verificationStatus"), args.filter!.status!)
        );
      }
      
      // Filter by verifier
      if (args.filter.verifierId) {
        query = query.filter(q => 
          q.eq(q.field("verifierId"), args.filter!.verifierId!)
        );
      }
      
      // Filter by requester
      if (args.filter.requesterId) {
        query = query.filter(q => 
          q.eq(q.field("requesterId"), args.filter!.requesterId!)
        );
      }
      
      // Filter by experiment
      if (args.filter.experimentId) {
        query = query.withIndex("by_experiment", q => 
          q.eq("experimentId", args.filter!.experimentId!)
        );
      }
      
      // Filter by date range
      if (args.filter.dateRange) {
        if (args.filter.dateRange.start !== undefined) {
          query = query.filter(q => 
            q.gte(q.field("requestDate"), args.filter!.dateRange!.start!)
          );
        }
        
        if (args.filter.dateRange.end !== undefined) {
          query = query.filter(q => 
            q.lte(q.field("requestDate"), args.filter!.dateRange!.end!)
          );
        }
      }
    }
    
    // Apply sorting
    if (args.options?.sortBy) {
      const direction = args.options.sortDirection === "desc" ? "desc" : "asc";
      
      const validSortFields = ["requestDate", "verificationDate", "verificationStatus", "overallRating"];
      const sortField = validSortFields.includes(args.options.sortBy) 
        ? args.options.sortBy 
        : "requestDate";
      
      query = query.order(sortField as any, direction);
    } else {
      // Default sort by request date, newest first
      query = query.order("requestDate", "desc");
    }
    
    // Apply pagination
    if (args.options?.cursor) {
      query = query.withCursor(args.options.cursor);
    }
    
    const limit = args.options?.limit || 50;
    query = query.take(limit);
    
    // Execute query
    let verifications = await query.collect();
    
    // Include experiment details if requested
    if (args.options?.includeExperiment) {
      const experimentIds = verifications
        .map(v => v.experimentId)
        .filter((id): id is Id<"enhancedExperiments"> => id !== undefined);
      
      if (experimentIds.length > 0) {
        const experimentsMap = new Map();
        
        // Fetch all experiments in one query
        for (const id of experimentIds) {
          const experiment = await ctx.db.get(id);
          if (experiment) {
            experimentsMap.set(id.toString(), {
              _id: experiment._id,
              name: experiment.name,
              description: experiment.description,
              status: experiment.status,
              protocolId: experiment.protocolId,
              conductedBy: experiment.conductedBy,
              date: experiment.date
            });
          }
        }
        
        // Add experiment details to verifications
        verifications = verifications.map(verification => {
          if (verification.experimentId) {
            const experiment = experimentsMap.get(verification.experimentId.toString());
            if (experiment) {
              return {
                ...verification,
                experiment
              };
            }
          }
          return verification;
        });
      }
    }
    
    return verifications;
  }
});

/**
 * Get verification statistics
 */
export const getVerificationStats = query({
  args: {
    dateRange: v.optional(v.object({
      start: v.optional(v.number()),
      end: v.optional(v.number())
    }))
  },
  handler: async (ctx, args) => {
    // Initialize stats object
    const stats = {
      statusCounts: {
        pending: 0,
        verified: 0,
        rejected: 0,
        needsRevision: 0,
        total: 0
      },
      qualityStats: {
        reproducibility: {
          average: 0,
          min: 10,
          max: 0
        },
        quality: {
          average: 0,
          min: 10,
          max: 0
        },
        documentation: {
          average: 0,
          min: 10,
          max: 0
        },
        overall: {
          average: 0,
          min: 10,
          max: 0
        }
      },
      totalVerifications: 0,
      verifiedPercentage: 0,
      averageVerificationTime: null,
      verificationsByMonth: {}
    };
    
    // Query all verifications with date range filter if provided
    let query = ctx.db.query("labVerifications");
    
    if (args.dateRange) {
      if (args.dateRange.start !== undefined) {
        query = query.filter(q => 
          q.gte(q.field("requestDate"), args.dateRange!.start!)
        );
      }
      
      if (args.dateRange.end !== undefined) {
        query = query.filter(q => 
          q.lte(q.field("requestDate"), args.dateRange!.end!)
        );
      }
    }
    
    const verifications = await query.collect();
    
    // Count verifications by status
    for (const verification of verifications) {
      // Increment total count
      stats.totalVerifications++;
      
      // Count by status
      switch (verification.verificationStatus) {
        case "pending":
          stats.statusCounts.pending++;
          break;
        case "verified":
          stats.statusCounts.verified++;
          break;
        case "rejected":
          stats.statusCounts.rejected++;
          break;
        case "needs_revision":
          stats.statusCounts.needsRevision++;
          break;
      }
      
      // Collect ratings for verified experiments
      if (verification.verificationStatus === "verified") {
        // Reproducibility rating
        if (verification.reproducibilityRating !== undefined) {
          stats.qualityStats.reproducibility.average += verification.reproducibilityRating;
          stats.qualityStats.reproducibility.min = Math.min(
            stats.qualityStats.reproducibility.min, 
            verification.reproducibilityRating
          );
          stats.qualityStats.reproducibility.max = Math.max(
            stats.qualityStats.reproducibility.max, 
            verification.reproducibilityRating
          );
        }
        
        // Quality rating
        if (verification.qualityRating !== undefined) {
          stats.qualityStats.quality.average += verification.qualityRating;
          stats.qualityStats.quality.min = Math.min(
            stats.qualityStats.quality.min, 
            verification.qualityRating
          );
          stats.qualityStats.quality.max = Math.max(
            stats.qualityStats.quality.max, 
            verification.qualityRating
          );
        }
        
        // Documentation rating
        if (verification.documentationRating !== undefined) {
          stats.qualityStats.documentation.average += verification.documentationRating;
          stats.qualityStats.documentation.min = Math.min(
            stats.qualityStats.documentation.min, 
            verification.documentationRating
          );
          stats.qualityStats.documentation.max = Math.max(
            stats.qualityStats.documentation.max, 
            verification.documentationRating
          );
        }
        
        // Overall rating
        if (verification.overallRating !== undefined) {
          stats.qualityStats.overall.average += verification.overallRating;
          stats.qualityStats.overall.min = Math.min(
            stats.qualityStats.overall.min, 
            verification.overallRating
          );
          stats.qualityStats.overall.max = Math.max(
            stats.qualityStats.overall.max, 
            verification.overallRating
          );
        }
      }
      
      // Calculate verification time
      if (verification.verificationDate && verification.requestDate) {
        const timeToVerify = (verification.verificationDate - verification.requestDate) / (1000 * 60 * 60); // hours
        if (stats.averageVerificationTime === null) {
          stats.averageVerificationTime = timeToVerify;
        } else {
          stats.averageVerificationTime = (stats.averageVerificationTime + timeToVerify) / 2;
        }
      }
      
      // Group by month
      if (verification.requestDate) {
        const date = new Date(verification.requestDate);
        const monthKey = `${date.getFullYear()}-${String(date.getMonth() + 1).padStart(2, '0')}`;
        
        if (!stats.verificationsByMonth[monthKey]) {
          stats.verificationsByMonth[monthKey] = 0;
        }
        
        stats.verificationsByMonth[monthKey]++;
      }
    }
    
    // Update status count total
    stats.statusCounts.total = stats.totalVerifications;
    
    // Calculate verified percentage
    if (stats.totalVerifications > 0) {
      stats.verifiedPercentage = (stats.statusCounts.verified / stats.totalVerifications) * 100;
    }
    
    // Calculate averages
    if (stats.statusCounts.verified > 0) {
      stats.qualityStats.reproducibility.average /= stats.statusCounts.verified;
      stats.qualityStats.quality.average /= stats.statusCounts.verified;
      stats.qualityStats.documentation.average /= stats.statusCounts.verified;
      stats.qualityStats.overall.average /= stats.statusCounts.verified;
    } else {
      stats.qualityStats.reproducibility = null;
      stats.qualityStats.quality = null;
      stats.qualityStats.documentation = null;
      stats.qualityStats.overall = null;
    }
    
    return stats;
  }
});

/**
 * Helper function to validate a status string
 */
function isValidStatus(status: string): boolean {
  return ["pending", "verified", "rejected", "needs_revision"].includes(status);
}

// Export previous API functions for backward compatibility
export const getByExperimentId = getVerificationByExperiment;
export const getAll = listVerifications;
export const getStats = getVerificationStats;
export const create = requestVerification;
export const update = updateVerification;
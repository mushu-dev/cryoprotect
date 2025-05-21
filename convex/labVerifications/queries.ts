/**
 * Lab Verification Query Functions
 * 
 * This module provides specialized query functions for retrieving lab verification data
 * in formats needed for scientific reporting and dashboard displays.
 */

import { query } from "../_generated/server";
import { v } from "convex/values";
import { Id } from "../_generated/dataModel";

/**
 * Get verification timeline (verifications over time)
 */
export const getVerificationTimeline = query({
  args: {
    startDate: v.number(),
    endDate: v.number(),
    interval: v.string(), // "day", "week", "month"
  },
  handler: async (ctx, args) => {
    const { startDate, endDate, interval } = args;
    
    const verifications = await ctx.db
      .query("labVerifications")
      .withIndex("by_date", (q) => 
        q.gte("verificationDate", startDate).lte("verificationDate", endDate)
      )
      .collect();
    
    // Group by time interval
    const timeline: Record<string, { date: string, count: number, verified: number, rejected: number, pending: number }> = {};
    
    for (const verification of verifications) {
      const date = new Date(verification.verificationDate);
      let key: string;
      
      if (interval === "day") {
        key = date.toISOString().split('T')[0]; // YYYY-MM-DD
      } else if (interval === "week") {
        // Get the Monday of the week
        const day = date.getDay();
        const diff = date.getDate() - day + (day === 0 ? -6 : 1); // Adjust for Sunday
        const monday = new Date(date);
        monday.setDate(diff);
        key = monday.toISOString().split('T')[0]; // YYYY-MM-DD of Monday
      } else if (interval === "month") {
        key = `${date.getFullYear()}-${String(date.getMonth() + 1).padStart(2, '0')}`; // YYYY-MM
      } else {
        // Default to day if invalid interval
        key = date.toISOString().split('T')[0]; // YYYY-MM-DD
      }
      
      if (!timeline[key]) {
        timeline[key] = {
          date: key,
          count: 0,
          verified: 0,
          rejected: 0,
          pending: 0
        };
      }
      
      timeline[key].count++;
      
      if (verification.verificationStatus === "verified") {
        timeline[key].verified++;
      } else if (verification.verificationStatus === "rejected") {
        timeline[key].rejected++;
      } else if (verification.verificationStatus === "pending") {
        timeline[key].pending++;
      }
    }
    
    // Convert to array and sort by date
    return Object.values(timeline).sort((a, b) => a.date.localeCompare(b.date));
  },
});

/**
 * Get verification with experiment details
 */
export const getVerificationWithExperimentDetails = query({
  args: {
    verificationId: v.id("labVerifications"),
  },
  handler: async (ctx, args) => {
    const verification = await ctx.db.get(args.verificationId);
    
    if (!verification) {
      return null;
    }
    
    // Get experiment details
    const experiment = await ctx.db.get(verification.experimentId);
    
    if (!experiment) {
      return {
        verification,
        experiment: null,
      };
    }
    
    // Get verifier details
    const verifier = await ctx.db.get(verification.verifierId);
    
    return {
      verification,
      experiment,
      verifier,
    };
  },
});

/**
 * Get verifications by quality score range
 */
export const getVerificationsByQualityScore = query({
  args: {
    minScore: v.optional(v.number()),
    maxScore: v.optional(v.number()),
  },
  handler: async (ctx, args) => {
    const { minScore = 1, maxScore = 10 } = args;
    
    const verifications = await ctx.db
      .query("labVerifications")
      .filter((q) => 
        q.and(
          q.neq(q.field("qualityScore"), undefined),
          q.gte(q.field("qualityScore"), minScore),
          q.lte(q.field("qualityScore"), maxScore)
        )
      )
      .collect();
    
    return verifications;
  },
});

/**
 * Get verification summary by experiment type
 */
export const getVerificationSummaryByExperimentType = query({
  handler: async (ctx) => {
    const verifications = await ctx.db.query("labVerifications").collect();
    
    const experimentIds = verifications.map(v => v.experimentId);
    
    // Get experiment details
    const experiments = await Promise.all(
      experimentIds.map(id => ctx.db.get(id))
    );
    
    // Group by experiment type
    const experimentTypes: Record<string, {
      typeName: string;
      totalCount: number;
      verifiedCount: number;
      verifiedPercentage: number;
      averageQualityScore: number | null;
    }> = {};
    
    for (let i = 0; i < verifications.length; i++) {
      const experiment = experiments[i];
      const verification = verifications[i];
      
      if (!experiment) continue;
      
      const typeName = experiment.experimentTypeId;
      
      if (!experimentTypes[typeName]) {
        experimentTypes[typeName] = {
          typeName,
          totalCount: 0,
          verifiedCount: 0,
          verifiedPercentage: 0,
          averageQualityScore: null,
        };
      }
      
      experimentTypes[typeName].totalCount++;
      
      if (verification.verificationStatus === "verified") {
        experimentTypes[typeName].verifiedCount++;
      }
      
      // Update quality score average
      if (verification.qualityScore !== undefined) {
        const currentAvg = experimentTypes[typeName].averageQualityScore;
        const currentCount = experimentTypes[typeName].totalCount;
        
        if (currentAvg === null) {
          experimentTypes[typeName].averageQualityScore = verification.qualityScore;
        } else {
          // Recalculate average
          experimentTypes[typeName].averageQualityScore = 
            (currentAvg * (currentCount - 1) + verification.qualityScore) / currentCount;
        }
      }
    }
    
    // Calculate verified percentages
    for (const type in experimentTypes) {
      const { totalCount, verifiedCount } = experimentTypes[type];
      experimentTypes[type].verifiedPercentage = (verifiedCount / totalCount) * 100;
    }
    
    return Object.values(experimentTypes);
  },
});
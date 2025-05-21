/**
 * Types for Lab Verification System
 * 
 * These types define the data structures for the lab verification system,
 * which ensures experimental reproducibility and quality control.
 */

import { Id } from "../_generated/dataModel";

/**
 * Verification status options
 */
export type VerificationStatus = "pending" | "verified" | "rejected" | "needs_revision";

/**
 * Lab verification object
 */
export interface LabVerification {
  _id: Id<"labVerifications">;
  experimentId: Id<"enhancedExperiments">;
  verificationStatus: VerificationStatus;
  verifierId?: Id<"users">;
  requesterId?: Id<"users">;
  
  // Equipment and methodology
  equipmentUsed: string[];
  methodologyDescription?: string;
  controlProcedures?: string;
  
  // Dates
  requestDate: number;
  verificationDate?: number;
  
  // Documentation
  verifierNotes?: string;
  requestNotes?: string;
  evidenceUrls?: string[];
  
  // Protocol verification
  reviewedProtocolSteps?: string[];
  
  // Quality ratings (1-10 scale)
  reproducibilityRating?: number;
  qualityRating?: number;
  documentationRating?: number;
  overallRating?: number;
  
  // Metadata
  createdAt: number;
  updatedAt: number;
}

/**
 * Lab verification with expanded experiment details
 */
export interface LabVerificationWithExperiment extends LabVerification {
  experiment: {
    _id: Id<"enhancedExperiments">;
    name: string;
    description?: string;
    status: string;
    protocolId?: Id<"protocols">;
    conductedBy?: Id<"users">;
    date?: number;
  };
}

/**
 * Stats for lab verifications
 */
export interface VerificationStats {
  statusCounts: {
    pending: number;
    verified: number;
    rejected: number;
    needsRevision: number;
    total: number;
  };
  
  qualityStats: {
    reproducibility: {
      average: number;
      min: number;
      max: number;
    } | null;
    quality: {
      average: number;
      min: number;
      max: number;
    } | null;
    documentation: {
      average: number;
      min: number;
      max: number;
    } | null;
    overall: {
      average: number;
      min: number;
      max: number;
    } | null;
  };
  
  totalVerifications: number;
  verifiedPercentage: number;
  
  // Time metrics
  averageVerificationTime: number | null; // In hours
  verificationsByMonth: Record<string, number>;
}

/**
 * Parameters for creating a verification request
 */
export interface CreateVerificationParams {
  experimentId: Id<"enhancedExperiments">;
  equipmentUsed: string[];
  requestNotes?: string;
  methodologyDescription?: string;
  controlProcedures?: string;
  evidenceUrls?: string[];
  reviewedProtocolSteps?: string[];
}

/**
 * Parameters for verifying an experiment
 */
export interface VerifyExperimentParams {
  verificationId: Id<"labVerifications">;
  verificationStatus: VerificationStatus;
  verifierNotes?: string;
  reproducibilityRating?: number;
  qualityRating?: number;
  documentationRating?: number;
  evidenceUrls?: string[];
}

/**
 * Parameters for updating a verification
 */
export interface UpdateVerificationParams {
  verificationId: Id<"labVerifications">;
  verificationStatus?: VerificationStatus;
  equipmentUsed?: string[];
  requestNotes?: string;
  verifierNotes?: string;
  methodologyDescription?: string;
  controlProcedures?: string;
  evidenceUrls?: string[];
  reviewedProtocolSteps?: string[];
  reproducibilityRating?: number;
  qualityRating?: number;
  documentationRating?: number;
}

/**
 * Filter options for lab verifications
 */
export interface LabVerificationFilter {
  status?: VerificationStatus;
  verifierId?: Id<"users">;
  requesterId?: Id<"users">;
  experimentId?: Id<"enhancedExperiments">;
  dateRange?: {
    start?: number;
    end?: number;
  };
}

/**
 * Query options for lab verifications
 */
export interface LabVerificationQueryOptions {
  limit?: number;
  cursor?: string;
  sortBy?: "requestDate" | "verificationDate" | "status" | "overallRating";
  sortDirection?: "asc" | "desc";
  includeExperiment?: boolean;
}
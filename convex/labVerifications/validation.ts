/**
 * Validation utilities for Lab Verification
 * 
 * This module provides validation functions for lab verification data
 * to ensure data integrity and consistency.
 */

import { DatabaseReader } from "../_generated/server";
import { Id } from "../_generated/dataModel";
import { 
  VerificationStatus, 
  CreateVerificationParams, 
  UpdateVerificationParams,
  VerifyExperimentParams
} from "./types";

/**
 * Validates a verification status value
 */
export function isValidVerificationStatus(status: string): status is VerificationStatus {
  return ["pending", "verified", "rejected", "needs_revision"].includes(status);
}

/**
 * Validates rating is in the valid range (1-10, integer)
 */
export function isValidRating(rating: number | undefined): boolean {
  if (rating === undefined) return true;
  return rating >= 1 && rating <= 10 && Number.isInteger(rating);
}

/**
 * Validates evidence URLs format
 */
export function validateEvidenceUrls(urls: string[] | undefined): string | null {
  if (!urls) return null;
  
  for (const url of urls) {
    try {
      new URL(url);
    } catch (e) {
      return `Invalid URL format: ${url}`;
    }
  }
  
  return null;
}

/**
 * Validate experiment exists
 */
export async function validateExperimentExists(
  db: DatabaseReader,
  experimentId: Id<"enhancedExperiments">
): Promise<void> {
  const experiment = await db.get(experimentId);
  if (!experiment) {
    throw new Error(`Experiment with ID ${experimentId} not found`);
  }
}

/**
 * Validate lab verification exists
 */
export async function validateVerificationExists(
  db: DatabaseReader,
  verificationId: Id<"labVerifications">
): Promise<void> {
  const verification = await db.get(verificationId);
  if (!verification) {
    throw new Error(`Lab verification with ID ${verificationId} not found`);
  }
}

/**
 * Validate verification creation parameters
 */
export function validateCreateVerification(
  params: CreateVerificationParams
): void {
  // Validate experiment ID
  if (!params.experimentId) {
    throw new Error("Experiment ID is required");
  }
  
  // Validate equipment used
  if (!params.equipmentUsed || !Array.isArray(params.equipmentUsed) || params.equipmentUsed.length === 0) {
    throw new Error("Equipment used is required and must be an array with at least one item");
  }
  
  // Validate evidence URLs if provided
  const urlError = validateEvidenceUrls(params.evidenceUrls);
  if (urlError) throw new Error(urlError);
  
  // Validate protocol steps if provided
  if (params.reviewedProtocolSteps && !Array.isArray(params.reviewedProtocolSteps)) {
    throw new Error("Reviewed protocol steps must be an array");
  }
}

/**
 * Validate verification update parameters
 */
export function validateUpdateVerification(
  params: UpdateVerificationParams
): void {
  // Validate verification ID
  if (!params.verificationId) {
    throw new Error("Verification ID is required");
  }
  
  // Validate status if provided
  if (params.verificationStatus && !isValidVerificationStatus(params.verificationStatus)) {
    throw new Error("Invalid verification status");
  }
  
  // Validate equipment used if provided
  if (params.equipmentUsed && (!Array.isArray(params.equipmentUsed) || params.equipmentUsed.length === 0)) {
    throw new Error("Equipment used must be an array with at least one item");
  }
  
  // Validate evidence URLs if provided
  const urlError = validateEvidenceUrls(params.evidenceUrls);
  if (urlError) throw new Error(urlError);
  
  // Validate protocol steps if provided
  if (params.reviewedProtocolSteps && !Array.isArray(params.reviewedProtocolSteps)) {
    throw new Error("Reviewed protocol steps must be an array");
  }
  
  // Validate ratings if provided
  if (params.reproducibilityRating !== undefined && !isValidRating(params.reproducibilityRating)) {
    throw new Error("Reproducibility rating must be an integer between 1 and 10");
  }
  
  if (params.qualityRating !== undefined && !isValidRating(params.qualityRating)) {
    throw new Error("Quality rating must be an integer between 1 and 10");
  }
  
  if (params.documentationRating !== undefined && !isValidRating(params.documentationRating)) {
    throw new Error("Documentation rating must be an integer between 1 and 10");
  }
}

/**
 * Validate verification parameters
 */
export function validateVerifyExperiment(
  params: VerifyExperimentParams
): void {
  // Validate verification ID
  if (!params.verificationId) {
    throw new Error("Verification ID is required");
  }
  
  // Validate status
  if (!params.verificationStatus) {
    throw new Error("Verification status is required");
  }
  
  if (!["verified", "rejected", "needs_revision"].includes(params.verificationStatus)) {
    throw new Error("Invalid verification status for verification action");
  }
  
  // Validate ratings if provided
  if (params.reproducibilityRating !== undefined && !isValidRating(params.reproducibilityRating)) {
    throw new Error("Reproducibility rating must be an integer between 1 and 10");
  }
  
  if (params.qualityRating !== undefined && !isValidRating(params.qualityRating)) {
    throw new Error("Quality rating must be an integer between 1 and 10");
  }
  
  if (params.documentationRating !== undefined && !isValidRating(params.documentationRating)) {
    throw new Error("Documentation rating must be an integer between 1 and 10");
  }
  
  // If status is 'verified', require ratings
  if (params.verificationStatus === "verified") {
    if (
      params.reproducibilityRating === undefined ||
      params.qualityRating === undefined ||
      params.documentationRating === undefined
    ) {
      throw new Error("All ratings (reproducibility, quality, documentation) are required for verified status");
    }
  }
  
  // Validate evidence URLs if provided
  const urlError = validateEvidenceUrls(params.evidenceUrls);
  if (urlError) throw new Error(urlError);
}

/**
 * Validate user has permission to verify experiments
 * 
 * This is a placeholder for a more complex permission system.
 * For a complete implementation, you would check user roles and permissions.
 */
export async function validateCanVerify(
  db: DatabaseReader,
  userId: Id<"users"> | undefined
): Promise<void> {
  if (!userId) {
    throw new Error("User authentication required");
  }
  
  // In a real system, you would check if the user has the 'verifier' role
  // For this demonstration, we'll allow any authenticated user to verify
}

/**
 * Calculate overall rating from individual ratings
 */
export function calculateOverallRating(
  reproducibilityRating?: number,
  qualityRating?: number,
  documentationRating?: number
): number | undefined {
  // Only calculate if all ratings are provided
  if (
    reproducibilityRating === undefined ||
    qualityRating === undefined ||
    documentationRating === undefined
  ) {
    return undefined;
  }
  
  // Weighted average: reproducibility (40%), quality (40%), documentation (20%)
  return Math.round(
    (reproducibilityRating * 0.4) +
    (qualityRating * 0.4) +
    (documentationRating * 0.2)
  );
}

/**
 * Check if verification status can be updated
 * 
 * Implements business rules for status transitions
 */
export function canUpdateVerificationStatus(
  currentStatus: VerificationStatus,
  newStatus: VerificationStatus
): boolean {
  // Rules for status transitions:
  // - 'pending' can transition to any other status
  // - 'needs_revision' can transition to 'pending' (after revisions) or 'verified'/'rejected' (final decision)
  // - 'verified' and 'rejected' are terminal states and cannot be changed
  
  if (currentStatus === "pending") {
    // Pending can transition to any state
    return true;
  } else if (currentStatus === "needs_revision") {
    // Needs revision can go to pending (after revisions) or verified/rejected (final decision)
    return ["pending", "verified", "rejected"].includes(newStatus);
  } else if (["verified", "rejected"].includes(currentStatus)) {
    // Terminal states cannot be changed
    return false;
  }
  
  return false;
}

/**
 * Validates a complete lab verification object
 * Returns error message if invalid, null if valid
 */
export function validateLabVerification(verification: {
  verificationStatus: string;
  equipmentUsed: string[];
  reproducibilityRating?: number;
  qualityRating?: number;
  documentationRating?: number;
  evidenceUrls?: string[];
}): string | null {
  if (!isValidVerificationStatus(verification.verificationStatus)) {
    return `Invalid verification status: ${verification.verificationStatus}`;
  }
  
  if (!verification.equipmentUsed || !Array.isArray(verification.equipmentUsed) || verification.equipmentUsed.length === 0) {
    return "Equipment used must be an array with at least one item";
  }
  
  if (verification.reproducibilityRating !== undefined && !isValidRating(verification.reproducibilityRating)) {
    return "Reproducibility rating must be between 1 and 10";
  }
  
  if (verification.qualityRating !== undefined && !isValidRating(verification.qualityRating)) {
    return "Quality rating must be between 1 and 10";
  }
  
  if (verification.documentationRating !== undefined && !isValidRating(verification.documentationRating)) {
    return "Documentation rating must be between 1 and 10";
  }
  
  const urlError = validateEvidenceUrls(verification.evidenceUrls);
  if (urlError) return urlError;
  
  return null;
}
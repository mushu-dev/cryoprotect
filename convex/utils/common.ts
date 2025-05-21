/**
 * Common utilities for Convex functions
 */

import { v } from "convex/values";
import { MutationCtx, QueryCtx } from "../_generated/server";

/**
 * Get the current timestamp in milliseconds
 */
export const getCurrentTimestamp = (): number => {
  return Date.now();
};

/**
 * Check if a user is authenticated
 * @param ctx - The Convex context (query or mutation)
 * @returns True if user is authenticated
 */
export const isAuthenticated = (ctx: QueryCtx | MutationCtx): boolean => {
  return ctx.auth.userId !== null;
};

/**
 * Validate that a user is authenticated
 * @param ctx - The Convex context (query or mutation)
 * @throws Error if user is not authenticated
 */
export const validateAuthenticated = (ctx: QueryCtx | MutationCtx): void => {
  if (!isAuthenticated(ctx)) {
    throw new Error("Authentication required");
  }
};

/**
 * Get a user's role from the database
 * @param ctx - The Convex context (query or mutation)
 * @returns The user's role or null if not found
 */
export const getUserRole = async (ctx: QueryCtx | MutationCtx): Promise<string | null> => {
  if (!isAuthenticated(ctx)) {
    return null;
  }
  
  const user = await ctx.db
    .query("users")
    .filter(q => q.eq(q.field("email"), ctx.auth.subject))
    .first();
  
  return user?.role || null;
};

/**
 * Validate that a user has one of the specified roles
 * @param ctx - The Convex context (query or mutation)
 * @param allowedRoles - Array of roles that are allowed
 * @throws Error if user doesn't have one of the allowed roles
 */
export const validateRole = async (
  ctx: QueryCtx | MutationCtx,
  allowedRoles: string[]
): Promise<void> => {
  validateAuthenticated(ctx);
  
  const role = await getUserRole(ctx);
  
  if (!role || !allowedRoles.includes(role)) {
    throw new Error("Insufficient permissions");
  }
};

/**
 * Common validators for function arguments
 */
export const validators = {
  id: v.id,
  optional: v.optional,
  string: v.string,
  number: v.number,
  boolean: v.boolean,
  any: v.any,
  
  // Custom validators
  pubchemCid: v.optional(v.string()),
  inchiKey: v.optional(v.string()),
  canonicalSmiles: v.optional(v.string()),
  formula: v.optional(v.string()),
  status: v.union(v.literal("active"), v.literal("deprecated"), v.literal("consolidated")),
};

/**
 * Functions for pagination
 */
export const pagination = {
  /**
   * Create parameters for pagination
   */
  params: {
    first: v.optional(v.number()),
    after: v.optional(v.string()),
  },
  
  /**
   * Default page size
   */
  defaultPageSize: 50,
  
  /**
   * Maximum page size
   */
  maxPageSize: 100,
};

/**
 * Create an audit log entry
 * @param ctx - The Convex mutation context
 * @param table - The table being modified
 * @param documentId - The ID of the document being modified
 * @param operation - The operation (create, update, delete)
 * @param previousValue - The previous value (for updates and deletes)
 * @param newValue - The new value (for creates and updates)
 * @param reason - Optional reason for the change
 */
export const createAuditLog = async (
  ctx: MutationCtx,
  table: string,
  documentId: string,
  operation: "create" | "update" | "delete",
  previousValue?: any,
  newValue?: any,
  reason?: string
): Promise<void> => {
  await ctx.db.insert("scientificDataAudit", {
    table,
    documentId: documentId as any, // Type cast to avoid ID type issues
    operation,
    userId: isAuthenticated(ctx) ? ctx.auth.userId : undefined,
    previousValue,
    newValue,
    timestamp: getCurrentTimestamp(),
    reason,
  });
};
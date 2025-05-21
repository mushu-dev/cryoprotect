/**
 * User management functions for authentication
 */

import { mutation, query } from "../_generated/server";
import { v } from "convex/values";
import { Id } from "../_generated/dataModel";
import { createAuditLog, getCurrentTimestamp, validateAuthenticated, validateRole } from "../utils/common";

/**
 * Get the current user
 */
export const getCurrentUser = query({
  args: {},
  handler: async (ctx): Promise<{
    user: any | null;
    isAuthenticated: boolean;
    roles: string[];
  }> => {
    const isAuthenticated = !!ctx.auth.userId;
    
    if (!isAuthenticated) {
      return {
        user: null,
        isAuthenticated: false,
        roles: [],
      };
    }
    
    // Get user by identity token subject
    const user = await ctx.db
      .query("users")
      .filter(q => q.eq(q.field("email"), ctx.auth.subject))
      .first();
    
    if (!user) {
      return {
        user: null,
        isAuthenticated: true,
        roles: [],
      };
    }
    
    // Update last login time
    await ctx.db.patch(user._id, {
      lastLogin: getCurrentTimestamp(),
    });
    
    return {
      user,
      isAuthenticated: true,
      roles: [user.role],
    };
  },
});

/**
 * Create or update user on login
 */
export const upsertUser = mutation({
  args: {
    email: v.string(),
    name: v.optional(v.string()),
  },
  handler: async (ctx, args): Promise<Id<"users">> => {
    validateAuthenticated(ctx);
    
    if (args.email !== ctx.auth.subject) {
      throw new Error("Email doesn't match authentication subject");
    }
    
    // Check if user already exists
    const existingUser = await ctx.db
      .query("users")
      .filter(q => q.eq(q.field("email"), args.email))
      .first();
    
    const now = getCurrentTimestamp();
    
    if (existingUser) {
      // Update existing user
      await ctx.db.patch(existingUser._id, {
        name: args.name || existingUser.name,
        lastLogin: now,
        updatedAt: now,
      });
      
      return existingUser._id;
    } else {
      // Create new user
      const userId = await ctx.db.insert("users", {
        email: args.email,
        name: args.name,
        role: "viewer", // Default role
        lastLogin: now,
        createdAt: now,
        updatedAt: now,
      });
      
      return userId;
    }
  },
});

/**
 * Update user role (admin only)
 */
export const updateUserRole = mutation({
  args: {
    userId: v.id("users"),
    role: v.union(
      v.literal("admin"),
      v.literal("scientist"),
      v.literal("viewer")
    ),
  },
  handler: async (ctx, args): Promise<any> => {
    // Only admins can update roles
    await validateRole(ctx, ["admin"]);
    
    // Get current user data for audit
    const currentUser = await ctx.db.get(args.userId);
    if (!currentUser) {
      throw new Error("User not found");
    }
    
    // Update role
    const updatedUser = await ctx.db.patch(args.userId, {
      role: args.role,
      updatedAt: getCurrentTimestamp(),
    });
    
    // Create audit log
    await createAuditLog(
      ctx,
      "users",
      args.userId.toString(),
      "update",
      currentUser,
      updatedUser,
      `Role updated from ${currentUser.role} to ${args.role}`
    );
    
    return updatedUser;
  },
});

/**
 * Get all users (admin only)
 */
export const listUsers = query({
  args: {},
  handler: async (ctx) => {
    // Only admins can list all users
    await validateRole(ctx, ["admin"]);
    
    // Get all users
    const users = await ctx.db.query("users").collect();
    return users;
  },
});

/**
 * Get user by ID (admin only)
 */
export const getUserById = query({
  args: {
    userId: v.id("users"),
  },
  handler: async (ctx, args) => {
    // Only admins can get other users by ID
    await validateRole(ctx, ["admin"]);
    
    // Get user
    const user = await ctx.db.get(args.userId);
    return user;
  },
});
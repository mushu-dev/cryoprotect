/**
 * CRUD operations for users
 */

import { mutation, query } from "../_generated/server";
import { v } from "convex/values";
import { Id } from "../_generated/dataModel";
import { 
  User, 
  CreateUserInput, 
  UpdateUserInput,
  UserFilter,
  UserQueryOptions
} from "./types";
import { 
  validateCreateUserInput, 
  validateUpdateUserInput,
  validateUserExists
} from "./validation";
import { 
  formatUserRole,
  expandUserWithProjects,
  isUserAdmin
} from "./helpers";
import { getCurrentUser } from "../auth/users";

/**
 * Create a new user
 */
export const createUser = mutation({
  args: {
    user: v.object({
      email: v.string(),
      name: v.optional(v.string()),
      role: v.optional(v.string())
    })
  },
  handler: async (ctx, args) => {
    // Validate input
    validateCreateUserInput(args.user);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Check if current user is an admin
    if (userId) {
      const isAdmin = await isUserAdmin(ctx, userId);
      
      // Only admins can create users with specific roles
      if (args.user.role && args.user.role !== "viewer" && !isAdmin) {
        throw new Error("Only administrators can create users with admin or scientist roles");
      }
    }
    
    // Check if user with same email already exists
    const existingUser = await ctx.db
      .query("users")
      .withIndex("by_email", q => q.eq("email", args.user.email))
      .first();
    
    if (existingUser) {
      throw new Error(`User with email ${args.user.email} already exists`);
    }
    
    // Prepare data for insertion
    const now = Date.now();
    const userData = {
      ...args.user,
      role: args.user.role || "viewer", // Default role is viewer
      createdAt: now,
      updatedAt: now
    };
    
    // Insert user
    const newUserId = await ctx.db.insert("users", userData);
    
    // Create audit log entry
    await ctx.db.insert("scientificDataAudit", {
      table: "users",
      documentId: newUserId,
      operation: "create",
      userId,
      newValue: userData,
      timestamp: now,
    });
    
    return newUserId;
  }
});

/**
 * Get a user by ID
 */
export const getUser = query({
  args: {
    userId: v.id("users"),
    includeProjects: v.optional(v.boolean())
  },
  handler: async (ctx, args) => {
    // Validate user exists
    validateUserExists(args.userId);
    
    // Get the user
    const user = await ctx.db.get(args.userId);
    if (!user) {
      return null;
    }
    
    // Include projects if requested
    if (args.includeProjects) {
      return expandUserWithProjects(ctx, user);
    }
    
    return user;
  }
});

/**
 * Update an existing user
 */
export const updateUser = mutation({
  args: {
    userId: v.id("users"),
    update: v.object({
      name: v.optional(v.string()),
      role: v.optional(v.string()),
      lastLogin: v.optional(v.number())
    })
  },
  handler: async (ctx, args) => {
    // Validate input
    validateUpdateUserInput(args.update);
    validateUserExists(args.userId);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const currentUserId = user?._id;
    
    // Only allow self-updates or admin updates
    if (!currentUserId) {
      throw new Error("Authentication required to update a user");
    }
    
    const isAdmin = await isUserAdmin(ctx, currentUserId);
    const isSelfUpdate = currentUserId.equals(args.userId);
    
    if (!isAdmin && !isSelfUpdate) {
      throw new Error("You do not have permission to update this user");
    }
    
    // Get existing user
    const existingUser = await ctx.db.get(args.userId);
    if (!existingUser) {
      throw new Error(`User with ID ${args.userId} not found`);
    }
    
    // Role changes require admin privileges
    if (args.update.role && args.update.role !== existingUser.role && !isAdmin) {
      throw new Error("Only administrators can change user roles");
    }
    
    // Prepare update data
    const now = Date.now();
    const updateData = {
      ...args.update,
      updatedAt: now
    };
    
    // Update user
    await ctx.db.patch(args.userId, updateData);
    
    // Create audit log entry
    await ctx.db.insert("scientificDataAudit", {
      table: "users",
      documentId: args.userId,
      operation: "update",
      userId: currentUserId,
      previousValue: existingUser,
      newValue: { ...existingUser, ...updateData },
      timestamp: now
    });
    
    return args.userId;
  }
});

/**
 * Delete a user (admin only)
 */
export const deleteUser = mutation({
  args: {
    userId: v.id("users")
  },
  handler: async (ctx, args) => {
    // Validate user exists
    validateUserExists(args.userId);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const currentUserId = user?._id;
    
    // Only admins can delete users
    if (!currentUserId) {
      throw new Error("Authentication required to delete a user");
    }
    
    const isAdmin = await isUserAdmin(ctx, currentUserId);
    
    if (!isAdmin) {
      throw new Error("Only administrators can delete users");
    }
    
    // Get existing user
    const existingUser = await ctx.db.get(args.userId);
    if (!existingUser) {
      throw new Error(`User with ID ${args.userId} not found`);
    }
    
    // Check for team memberships
    const teamMemberships = await ctx.db
      .query("teamMembers")
      .withIndex("by_user", q => q.eq("userId", args.userId))
      .collect();
    
    // Remove user from all teams
    for (const membership of teamMemberships) {
      await ctx.db.delete(membership._id);
    }
    
    // Check for owned projects
    const ownedProjects = await ctx.db
      .query("projects")
      .withIndex("by_owner", q => q.eq("ownerId", args.userId))
      .collect();
    
    if (ownedProjects.length > 0) {
      throw new Error(
        `Cannot delete user with ${ownedProjects.length} owned projects. ` +
        `Transfer ownership of these projects first.`
      );
    }
    
    // Create audit log entry
    const now = Date.now();
    await ctx.db.insert("scientificDataAudit", {
      table: "users",
      documentId: args.userId,
      operation: "delete",
      userId: currentUserId,
      previousValue: existingUser,
      timestamp: now
    });
    
    // Delete the user
    await ctx.db.delete(args.userId);
    
    return true;
  }
});

/**
 * List users with optional filtering
 */
export const listUsers = query({
  args: {
    filter: v.optional(v.object({
      email: v.optional(v.string()),
      role: v.optional(v.string()),
      includeProjects: v.optional(v.boolean())
    })),
    options: v.optional(v.object({
      limit: v.optional(v.number()),
      cursor: v.optional(v.string()),
      sortBy: v.optional(v.string()),
      sortDirection: v.optional(v.string())
    }))
  },
  handler: async (ctx, args) => {
    // Set up the query
    let query = ctx.db.query("users");
    
    // Apply filters
    if (args.filter) {
      if (args.filter.email) {
        query = query.withIndex("by_email", q => 
          q.eq("email", args.filter!.email!)
        );
      }
      
      if (args.filter.role) {
        query = query.withIndex("by_role", q => 
          q.eq("role", args.filter!.role!)
        );
      }
    }
    
    // Apply sorting
    if (args.options?.sortBy) {
      const sortDirection = args.options.sortDirection === "desc" ? "desc" : "asc";
      
      switch (args.options.sortBy) {
        case "email":
          query = query.order("email", sortDirection);
          break;
        case "role":
          query = query.order("role", sortDirection);
          break;
        case "createdAt":
          query = query.order("createdAt", sortDirection);
          break;
        default:
          query = query.order("email", "asc"); // default sort
      }
    } else {
      // Default sort by email ascending
      query = query.order("email", "asc");
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
    const users = await query.collect();
    
    // Include projects if requested
    if (args.filter?.includeProjects) {
      const usersWithProjects = [];
      
      for (const user of users) {
        usersWithProjects.push(await expandUserWithProjects(ctx, user));
      }
      
      return usersWithProjects;
    }
    
    return users;
  }
});

/**
 * Search users by email or name
 */
export const searchUsers = query({
  args: {
    query: v.string(),
    limit: v.optional(v.number())
  },
  handler: async (ctx, args) => {
    // Set up the query
    let query = ctx.db.query("users");
    
    // Apply search filter (simple implementation)
    // In a real implementation, we would use full-text search
    query = query.filter(q => 
      q.or(
        q.includes(q.field("email"), args.query),
        q.includes(q.field("name"), args.query)
      )
    );
    
    // Apply limit
    const limit = args.limit || 20;
    query = query.take(limit);
    
    // Execute query
    return query.collect();
  }
});

/**
 * Get the current user's profile
 */
export const getCurrentUserProfile = query({
  args: {
    includeProjects: v.optional(v.boolean())
  },
  handler: async (ctx, args) => {
    // Get current user
    const user = await getCurrentUser(ctx);
    
    if (!user) {
      return null;
    }
    
    // Include projects if requested
    if (args.includeProjects) {
      return expandUserWithProjects(ctx, user);
    }
    
    return user;
  }
});

/**
 * Update last login time for a user
 */
export const updateLastLogin = mutation({
  args: {
    userId: v.id("users")
  },
  handler: async (ctx, args) => {
    // Validate user exists
    validateUserExists(args.userId);
    
    // Get existing user
    const existingUser = await ctx.db.get(args.userId);
    if (!existingUser) {
      throw new Error(`User with ID ${args.userId} not found`);
    }
    
    // Update lastLogin
    const now = Date.now();
    await ctx.db.patch(args.userId, {
      lastLogin: now,
      updatedAt: now
    });
    
    return args.userId;
  }
});
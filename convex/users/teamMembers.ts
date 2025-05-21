/**
 * CRUD operations for team members
 */

import { mutation, query } from "../_generated/server";
import { v } from "convex/values";
import { Id } from "../_generated/dataModel";
import { 
  TeamMember, 
  CreateTeamMemberInput, 
  UpdateTeamMemberInput,
  TeamMemberFilter,
  TeamMemberQueryOptions
} from "./types";
import { 
  validateCreateTeamMemberInput, 
  validateUpdateTeamMemberInput,
  validateTeamMemberExists,
  validateProjectExists,
  validateUserExists,
  validateProjectAdminAccess
} from "./validation";
import { 
  formatTeamRole,
  getProjectTeamRoles
} from "./helpers";
import { getCurrentUser } from "../auth/users";

/**
 * Add a team member to a project
 */
export const addTeamMember = mutation({
  args: {
    teamMember: v.object({
      projectId: v.id("projects"),
      userId: v.id("users"),
      role: v.optional(v.string())
    })
  },
  handler: async (ctx, args) => {
    // Validate input
    validateCreateTeamMemberInput(args.teamMember);
    validateProjectExists(args.teamMember.projectId);
    validateUserExists(args.teamMember.userId);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const currentUserId = user?._id;
    
    // User must be logged in to add a team member
    if (!currentUserId) {
      throw new Error("Authentication required to add a team member");
    }
    
    // Get project
    const project = await ctx.db.get(args.teamMember.projectId);
    if (!project) {
      throw new Error(`Project with ID ${args.teamMember.projectId} not found`);
    }
    
    // Get team roles
    const teamRoles = await getProjectTeamRoles(ctx, args.teamMember.projectId);
    
    // Check admin access
    validateProjectAdminAccess(
      args.teamMember.projectId, 
      currentUserId, 
      project.ownerId,
      teamRoles
    );
    
    // Verify user exists
    const userToAdd = await ctx.db.get(args.teamMember.userId);
    if (!userToAdd) {
      throw new Error(`User with ID ${args.teamMember.userId} not found`);
    }
    
    // Check if user is already a team member
    const existingMember = await ctx.db
      .query("teamMembers")
      .withIndex("by_project", q => q.eq("projectId", args.teamMember.projectId))
      .filter(q => q.eq(q.field("userId"), args.teamMember.userId))
      .first();
    
    if (existingMember) {
      // Update role if different
      if (args.teamMember.role && existingMember.role !== args.teamMember.role) {
        await ctx.db.patch(existingMember._id, { 
          role: args.teamMember.role,
          updatedAt: Date.now()
        });
      }
      
      return existingMember._id;
    }
    
    // Prepare data for insertion
    const now = Date.now();
    const teamMemberData = {
      ...args.teamMember,
      role: args.teamMember.role || "viewer", // Default role is viewer
      createdAt: now,
      updatedAt: now
    };
    
    // Insert team member
    const teamMemberId = await ctx.db.insert("teamMembers", teamMemberData);
    
    // Create audit log entry
    await ctx.db.insert("scientificDataAudit", {
      table: "teamMembers",
      documentId: teamMemberId,
      operation: "create",
      userId: currentUserId,
      newValue: teamMemberData,
      timestamp: now,
    });
    
    return teamMemberId;
  }
});

/**
 * Get a team member by ID
 */
export const getTeamMember = query({
  args: {
    teamMemberId: v.id("teamMembers")
  },
  handler: async (ctx, args) => {
    // Validate team member exists
    validateTeamMemberExists(args.teamMemberId);
    
    // Get the team member
    const teamMember = await ctx.db.get(args.teamMemberId);
    if (!teamMember) {
      return null;
    }
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const currentUserId = user?._id;
    
    // Get project
    const project = await ctx.db.get(teamMember.projectId);
    if (!project) {
      return null;
    }
    
    // Get team roles
    const teamRoles = await getProjectTeamRoles(ctx, teamMember.projectId);
    
    // Check project access
    try {
      validateProjectAdminAccess(
        teamMember.projectId, 
        currentUserId, 
        project.ownerId,
        teamRoles
      );
    } catch (error) {
      // If not admin, check if user is viewing their own membership
      if (!currentUserId || !teamMember.userId.equals(currentUserId)) {
        throw new Error("You do not have permission to view this team member");
      }
    }
    
    return teamMember;
  }
});

/**
 * Update a team member's role
 */
export const updateTeamMember = mutation({
  args: {
    teamMemberId: v.id("teamMembers"),
    update: v.object({
      role: v.string()
    })
  },
  handler: async (ctx, args) => {
    // Validate input
    validateUpdateTeamMemberInput(args.update);
    validateTeamMemberExists(args.teamMemberId);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const currentUserId = user?._id;
    
    // User must be logged in to update a team member
    if (!currentUserId) {
      throw new Error("Authentication required to update a team member");
    }
    
    // Get existing team member
    const existingMember = await ctx.db.get(args.teamMemberId);
    if (!existingMember) {
      throw new Error(`Team member with ID ${args.teamMemberId} not found`);
    }
    
    // Get project
    const project = await ctx.db.get(existingMember.projectId);
    if (!project) {
      throw new Error(`Project with ID ${existingMember.projectId} not found`);
    }
    
    // Get team roles
    const teamRoles = await getProjectTeamRoles(ctx, existingMember.projectId);
    
    // Check admin access
    validateProjectAdminAccess(
      existingMember.projectId, 
      currentUserId, 
      project.ownerId,
      teamRoles
    );
    
    // Prepare update data
    const now = Date.now();
    const updateData = {
      ...args.update,
      updatedAt: now
    };
    
    // Update team member
    await ctx.db.patch(args.teamMemberId, updateData);
    
    // Create audit log entry
    await ctx.db.insert("scientificDataAudit", {
      table: "teamMembers",
      documentId: args.teamMemberId,
      operation: "update",
      userId: currentUserId,
      previousValue: existingMember,
      newValue: { ...existingMember, ...updateData },
      timestamp: now
    });
    
    return args.teamMemberId;
  }
});

/**
 * Remove a team member from a project
 */
export const removeTeamMember = mutation({
  args: {
    teamMemberId: v.id("teamMembers")
  },
  handler: async (ctx, args) => {
    // Validate team member exists
    validateTeamMemberExists(args.teamMemberId);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const currentUserId = user?._id;
    
    // User must be logged in to remove a team member
    if (!currentUserId) {
      throw new Error("Authentication required to remove a team member");
    }
    
    // Get existing team member
    const existingMember = await ctx.db.get(args.teamMemberId);
    if (!existingMember) {
      throw new Error(`Team member with ID ${args.teamMemberId} not found`);
    }
    
    // Get project
    const project = await ctx.db.get(existingMember.projectId);
    if (!project) {
      throw new Error(`Project with ID ${existingMember.projectId} not found`);
    }
    
    // Get team roles
    const teamRoles = await getProjectTeamRoles(ctx, existingMember.projectId);
    
    // Check admin access (project owner or admin team member)
    try {
      validateProjectAdminAccess(
        existingMember.projectId, 
        currentUserId, 
        project.ownerId,
        teamRoles
      );
    } catch (error) {
      // If not admin, check if user is removing themselves
      if (!existingMember.userId.equals(currentUserId)) {
        throw new Error("You do not have permission to remove this team member");
      }
    }
    
    // Cannot remove the owner from their own project
    if (existingMember.userId.equals(project.ownerId)) {
      throw new Error("Cannot remove the project owner from the team");
    }
    
    // Create audit log entry
    const now = Date.now();
    await ctx.db.insert("scientificDataAudit", {
      table: "teamMembers",
      documentId: args.teamMemberId,
      operation: "delete",
      userId: currentUserId,
      previousValue: existingMember,
      timestamp: now
    });
    
    // Delete the team member
    await ctx.db.delete(args.teamMemberId);
    
    return true;
  }
});

/**
 * List team members for a project
 */
export const listTeamMembers = query({
  args: {
    projectId: v.id("projects"),
    includeOwner: v.optional(v.boolean())
  },
  handler: async (ctx, args) => {
    // Validate project exists
    validateProjectExists(args.projectId);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const currentUserId = user?._id;
    
    // Get project
    const project = await ctx.db.get(args.projectId);
    if (!project) {
      throw new Error(`Project with ID ${args.projectId} not found`);
    }
    
    // Get team roles
    const teamRoles = await getProjectTeamRoles(ctx, args.projectId);
    
    // Check project access
    validateProjectAdminAccess(
      args.projectId, 
      currentUserId, 
      project.ownerId,
      teamRoles
    );
    
    // Get team members
    const teamMembers = await ctx.db
      .query("teamMembers")
      .withIndex("by_project", q => q.eq("projectId", args.projectId))
      .collect();
    
    // Include project owner if requested and not already in team
    if (args.includeOwner) {
      const ownerInTeam = teamMembers.some(member => member.userId.equals(project.ownerId));
      
      if (!ownerInTeam) {
        const now = Date.now();
        teamMembers.push({
          _id: project._id as Id<"teamMembers">, // This is a hack for the frontend
          _creationTime: project._creationTime,
          projectId: project._id,
          userId: project.ownerId,
          role: "owner",
          createdAt: project.createdAt,
          updatedAt: now
        });
      }
    }
    
    return teamMembers;
  }
});

/**
 * Batch add multiple team members to a project
 */
export const batchAddTeamMembers = mutation({
  args: {
    projectId: v.id("projects"),
    members: v.array(
      v.object({
        userId: v.id("users"),
        role: v.optional(v.string())
      })
    )
  },
  handler: async (ctx, args) => {
    // Validate project exists
    validateProjectExists(args.projectId);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const currentUserId = user?._id;
    
    // User must be logged in to add team members
    if (!currentUserId) {
      throw new Error("Authentication required to add team members");
    }
    
    // Get project
    const project = await ctx.db.get(args.projectId);
    if (!project) {
      throw new Error(`Project with ID ${args.projectId} not found`);
    }
    
    // Get team roles
    const teamRoles = await getProjectTeamRoles(ctx, args.projectId);
    
    // Check admin access
    validateProjectAdminAccess(
      args.projectId, 
      currentUserId, 
      project.ownerId,
      teamRoles
    );
    
    // Record time for all operations
    const now = Date.now();
    
    // Get existing team members to avoid duplicates
    const existingMembers = await ctx.db
      .query("teamMembers")
      .withIndex("by_project", q => q.eq("projectId", args.projectId))
      .collect();
    
    // Create a map of userId -> team member
    const existingMembersMap = new Map<string, TeamMember>();
    
    for (const member of existingMembers) {
      existingMembersMap.set(member.userId.toString(), member);
    }
    
    // Validate and insert/update each team member
    const teamMemberIds: Id<"teamMembers">[] = [];
    
    for (const member of args.members) {
      // Validate team member input
      validateCreateTeamMemberInput({
        ...member,
        projectId: args.projectId
      });
      
      // Verify user exists
      const userToAdd = await ctx.db.get(member.userId);
      if (!userToAdd) {
        throw new Error(`User with ID ${member.userId} not found`);
      }
      
      // Skip adding the project owner as a regular team member
      if (member.userId.equals(project.ownerId)) {
        continue;
      }
      
      // Check if user is already a team member
      const existingMember = existingMembersMap.get(member.userId.toString());
      
      if (existingMember) {
        // Update role if different
        if (member.role && existingMember.role !== member.role) {
          await ctx.db.patch(existingMember._id, { 
            role: member.role,
            updatedAt: now
          });
        }
        
        teamMemberIds.push(existingMember._id);
      } else {
        // Prepare team member data
        const teamMemberData = {
          projectId: args.projectId,
          userId: member.userId,
          role: member.role || "viewer", // Default role is viewer
          createdAt: now,
          updatedAt: now
        };
        
        // Insert team member
        const teamMemberId = await ctx.db.insert("teamMembers", teamMemberData);
        
        // Create audit log entry
        await ctx.db.insert("scientificDataAudit", {
          table: "teamMembers",
          documentId: teamMemberId,
          operation: "create",
          userId: currentUserId,
          newValue: teamMemberData,
          timestamp: now
        });
        
        teamMemberIds.push(teamMemberId);
      }
    }
    
    return teamMemberIds;
  }
});
/**
 * Helper functions for users, projects, and team members
 */

import { Id } from "../_generated/dataModel";
import { QueryCtx, MutationCtx } from "../_generated/server";
import { 
  User, 
  Project,
  TeamMember,
  UserWithProjects,
  ProjectWithTeam
} from "./types";

/**
 * Format project status for display
 */
export function formatProjectStatus(status: string): string {
  switch (status) {
    case "active":
      return "Active";
    case "archived":
      return "Archived";
    case "completed":
      return "Completed";
    default:
      return status.charAt(0).toUpperCase() + status.slice(1);
  }
}

/**
 * Format user role for display
 */
export function formatUserRole(role: string): string {
  switch (role) {
    case "admin":
      return "Administrator";
    case "scientist":
      return "Scientist";
    case "viewer":
      return "Viewer";
    default:
      return role.charAt(0).toUpperCase() + role.slice(1);
  }
}

/**
 * Format team member role for display
 */
export function formatTeamRole(role: string): string {
  switch (role) {
    case "owner":
      return "Owner";
    case "editor":
      return "Editor";
    case "viewer":
      return "Viewer";
    default:
      return role.charAt(0).toUpperCase() + role.slice(1);
  }
}

/**
 * Get team roles for a project as a map of userId -> role
 */
export async function getProjectTeamRoles(
  ctx: QueryCtx,
  projectId: Id<"projects">
): Promise<Record<string, string>> {
  const teamMembers = await ctx.db
    .query("teamMembers")
    .withIndex("by_project", q => q.eq("projectId", projectId))
    .collect();
  
  // Create a map of userId -> role
  const teamRoles: Record<string, string> = {};
  
  for (const member of teamMembers) {
    teamRoles[member.userId.toString()] = member.role;
  }
  
  return teamRoles;
}

/**
 * Expand a user with their projects
 */
export async function expandUserWithProjects(
  ctx: QueryCtx,
  user: User
): Promise<UserWithProjects> {
  // Get projects owned by user
  const ownedProjects = await ctx.db
    .query("projects")
    .withIndex("by_owner", q => q.eq("ownerId", user._id))
    .collect();
  
  // Get team memberships
  const teamMemberships = await ctx.db
    .query("teamMembers")
    .withIndex("by_user", q => q.eq("userId", user._id))
    .collect();
  
  // Get projects user is a member of
  const memberProjects = [];
  
  for (const membership of teamMemberships) {
    const project = await ctx.db.get(membership.projectId);
    
    if (project) {
      memberProjects.push({
        project,
        role: membership.role
      });
    }
  }
  
  return {
    ...user,
    projects: {
      owned: ownedProjects,
      member: memberProjects
    }
  };
}

/**
 * Expand a project with its team
 */
export async function expandProjectWithTeam(
  ctx: QueryCtx,
  project: Project
): Promise<ProjectWithTeam> {
  // Get team members
  const teamMembers = await ctx.db
    .query("teamMembers")
    .withIndex("by_project", q => q.eq("projectId", project._id))
    .collect();
  
  // Get user details for each team member
  const team = [];
  
  for (const member of teamMembers) {
    const user = await ctx.db.get(member.userId);
    
    if (user) {
      team.push({
        userId: user._id,
        userEmail: user.email,
        userName: user.name,
        role: member.role
      });
    }
  }
  
  // Add owner if not already in team
  const ownerInTeam = team.some(member => member.userId.equals(project.ownerId));
  
  if (!ownerInTeam) {
    const owner = await ctx.db.get(project.ownerId);
    
    if (owner) {
      team.push({
        userId: owner._id,
        userEmail: owner.email,
        userName: owner.name,
        role: "owner"
      });
    }
  }
  
  return {
    ...project,
    team
  };
}

/**
 * Ensure a user is a team member (with specified role)
 */
export async function ensureTeamMember(
  ctx: MutationCtx,
  projectId: Id<"projects">,
  userId: Id<"users">,
  role: string
): Promise<Id<"teamMembers">> {
  // Check if user is already a team member
  const existing = await ctx.db
    .query("teamMembers")
    .withIndex("by_project", q => q.eq("projectId", projectId))
    .filter(q => q.eq(q.field("userId"), userId))
    .first();
  
  if (existing) {
    // Update role if different
    if (existing.role !== role) {
      await ctx.db.patch(existing._id, { 
        role,
        updatedAt: Date.now()
      });
    }
    
    return existing._id;
  }
  
  // Otherwise, create a new team member
  const now = Date.now();
  
  return await ctx.db.insert("teamMembers", {
    projectId,
    userId,
    role,
    createdAt: now,
    updatedAt: now
  });
}

/**
 * Check if user is a system administrator
 */
export async function isUserAdmin(
  ctx: QueryCtx,
  userId: Id<"users">
): Promise<boolean> {
  const user = await ctx.db.get(userId);
  return !!user && user.role === "admin";
}

/**
 * Get all projects a user has access to
 */
export async function getUserAccessibleProjects(
  ctx: QueryCtx,
  userId: Id<"users">
): Promise<Project[]> {
  // Get projects owned by user
  const ownedProjects = await ctx.db
    .query("projects")
    .withIndex("by_owner", q => q.eq("ownerId", userId))
    .collect();
  
  // Get team memberships
  const teamMemberships = await ctx.db
    .query("teamMembers")
    .withIndex("by_user", q => q.eq("userId", userId))
    .collect();
  
  // Get projects user is a member of
  const memberProjectIds = new Set<string>();
  
  for (const membership of teamMemberships) {
    memberProjectIds.add(membership.projectId.toString());
  }
  
  // Get public projects
  const publicProjects = await ctx.db
    .query("projects")
    .withIndex("by_public", q => q.eq("public", true))
    .collect();
  
  // Combine all projects, avoiding duplicates
  const allProjects = [...ownedProjects];
  const addedIds = new Set(ownedProjects.map(p => p._id.toString()));
  
  for (const project of publicProjects) {
    const projectId = project._id.toString();
    
    if (!addedIds.has(projectId)) {
      allProjects.push(project);
      addedIds.add(projectId);
    }
  }
  
  // Add projects where user is a team member
  for (const membership of teamMemberships) {
    const projectId = membership.projectId.toString();
    
    if (!addedIds.has(projectId)) {
      const project = await ctx.db.get(membership.projectId);
      
      if (project) {
        allProjects.push(project);
        addedIds.add(projectId);
      }
    }
  }
  
  return allProjects;
}
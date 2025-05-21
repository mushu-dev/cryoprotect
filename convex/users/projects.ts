/**
 * CRUD operations for projects
 */

import { mutation, query } from "../_generated/server";
import { v } from "convex/values";
import { Id } from "../_generated/dataModel";
import { 
  Project, 
  CreateProjectInput, 
  UpdateProjectInput,
  ProjectFilter,
  ProjectQueryOptions
} from "./types";
import { 
  validateCreateProjectInput, 
  validateUpdateProjectInput,
  validateProjectExists,
  validateProjectAccess,
  validateProjectEditAccess,
  validateProjectAdminAccess
} from "./validation";
import { 
  formatProjectStatus,
  expandProjectWithTeam,
  ensureTeamMember,
  getProjectTeamRoles,
  getUserAccessibleProjects
} from "./helpers";
import { getCurrentUser } from "../auth/users";

/**
 * Create a new project
 */
export const createProject = mutation({
  args: {
    project: v.object({
      name: v.string(),
      description: v.optional(v.string()),
      status: v.optional(v.string()),
      public: v.optional(v.boolean())
    })
  },
  handler: async (ctx, args) => {
    // Validate input
    validateCreateProjectInput(args.project);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // User must be logged in to create a project
    if (!userId) {
      throw new Error("Authentication required to create a project");
    }
    
    // Prepare data for insertion
    const now = Date.now();
    const projectData = {
      ...args.project,
      ownerId: userId,
      status: args.project.status || "active", // Default status is active
      createdAt: now,
      updatedAt: now,
      public: args.project.public ?? false
    };
    
    // Insert project
    const projectId = await ctx.db.insert("projects", projectData);
    
    // Add owner as a team member
    await ensureTeamMember(ctx, projectId, userId, "owner");
    
    // Create audit log entry
    await ctx.db.insert("scientificDataAudit", {
      table: "projects",
      documentId: projectId,
      operation: "create",
      userId,
      newValue: projectData,
      timestamp: now,
    });
    
    return projectId;
  }
});

/**
 * Get a project by ID
 */
export const getProject = query({
  args: {
    projectId: v.id("projects"),
    includeTeam: v.optional(v.boolean())
  },
  handler: async (ctx, args) => {
    // Validate project exists
    validateProjectExists(args.projectId);
    
    // Get the project
    const project = await ctx.db.get(args.projectId);
    if (!project) {
      return null;
    }
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Get team roles
    const teamRoles = await getProjectTeamRoles(ctx, args.projectId);
    
    // Check access
    validateProjectAccess(
      args.projectId, 
      userId, 
      project.public, 
      project.ownerId,
      teamRoles
    );
    
    // Include team if requested
    if (args.includeTeam) {
      return expandProjectWithTeam(ctx, project);
    }
    
    return project;
  }
});

/**
 * Update an existing project
 */
export const updateProject = mutation({
  args: {
    projectId: v.id("projects"),
    update: v.object({
      name: v.optional(v.string()),
      description: v.optional(v.string()),
      status: v.optional(v.string()),
      public: v.optional(v.boolean())
    })
  },
  handler: async (ctx, args) => {
    // Validate input
    validateUpdateProjectInput(args.update);
    validateProjectExists(args.projectId);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Get existing project
    const existingProject = await ctx.db.get(args.projectId);
    if (!existingProject) {
      throw new Error(`Project with ID ${args.projectId} not found`);
    }
    
    // Get team roles
    const teamRoles = await getProjectTeamRoles(ctx, args.projectId);
    
    // Check edit access
    validateProjectEditAccess(
      args.projectId, 
      userId, 
      existingProject.ownerId,
      teamRoles
    );
    
    // Prepare update data
    const now = Date.now();
    const updateData = {
      ...args.update,
      updatedAt: now
    };
    
    // Update project
    await ctx.db.patch(args.projectId, updateData);
    
    // Create audit log entry
    await ctx.db.insert("scientificDataAudit", {
      table: "projects",
      documentId: args.projectId,
      operation: "update",
      userId,
      previousValue: existingProject,
      newValue: { ...existingProject, ...updateData },
      timestamp: now
    });
    
    return args.projectId;
  }
});

/**
 * Transfer project ownership
 */
export const transferProjectOwnership = mutation({
  args: {
    projectId: v.id("projects"),
    newOwnerId: v.id("users")
  },
  handler: async (ctx, args) => {
    // Validate project exists
    validateProjectExists(args.projectId);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Get existing project
    const existingProject = await ctx.db.get(args.projectId);
    if (!existingProject) {
      throw new Error(`Project with ID ${args.projectId} not found`);
    }
    
    // Get team roles
    const teamRoles = await getProjectTeamRoles(ctx, args.projectId);
    
    // Check admin access
    validateProjectAdminAccess(
      args.projectId, 
      userId, 
      existingProject.ownerId,
      teamRoles
    );
    
    // Verify new owner exists
    const newOwner = await ctx.db.get(args.newOwnerId);
    if (!newOwner) {
      throw new Error(`User with ID ${args.newOwnerId} not found`);
    }
    
    // Prepare update data
    const now = Date.now();
    const updateData = {
      ownerId: args.newOwnerId,
      updatedAt: now
    };
    
    // Update project
    await ctx.db.patch(args.projectId, updateData);
    
    // Make sure new owner is in the team with owner role
    await ensureTeamMember(ctx, args.projectId, args.newOwnerId, "owner");
    
    // Create audit log entry
    await ctx.db.insert("scientificDataAudit", {
      table: "projects",
      documentId: args.projectId,
      operation: "update",
      userId,
      previousValue: existingProject,
      newValue: { ...existingProject, ...updateData },
      timestamp: now,
      reason: `Ownership transferred from ${existingProject.ownerId} to ${args.newOwnerId}`
    });
    
    return args.projectId;
  }
});

/**
 * Delete a project (owner only)
 */
export const deleteProject = mutation({
  args: {
    projectId: v.id("projects")
  },
  handler: async (ctx, args) => {
    // Validate project exists
    validateProjectExists(args.projectId);
    
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Get existing project
    const existingProject = await ctx.db.get(args.projectId);
    if (!existingProject) {
      throw new Error(`Project with ID ${args.projectId} not found`);
    }
    
    // Only the owner can delete a project
    if (!userId || !existingProject.ownerId.equals(userId)) {
      throw new Error("Only the project owner can delete a project");
    }
    
    // Get team members
    const teamMembers = await ctx.db
      .query("teamMembers")
      .withIndex("by_project", q => q.eq("projectId", args.projectId))
      .collect();
    
    // Delete all team members
    for (const member of teamMembers) {
      await ctx.db.delete(member._id);
    }
    
    // Check for project resources (experiments, mixtures, etc.)
    // In a real implementation, we would check for related resources
    // and either delete them or prevent deletion if resources exist
    
    // Create audit log entry
    const now = Date.now();
    await ctx.db.insert("scientificDataAudit", {
      table: "projects",
      documentId: args.projectId,
      operation: "delete",
      userId,
      previousValue: existingProject,
      timestamp: now
    });
    
    // Delete the project
    await ctx.db.delete(args.projectId);
    
    return true;
  }
});

/**
 * List projects with optional filtering
 */
export const listProjects = query({
  args: {
    filter: v.optional(v.object({
      name: v.optional(v.string()),
      ownerId: v.optional(v.id("users")),
      status: v.optional(v.string()),
      public: v.optional(v.boolean()),
      includeTeam: v.optional(v.boolean()),
      memberOf: v.optional(v.id("users"))
    })),
    options: v.optional(v.object({
      limit: v.optional(v.number()),
      cursor: v.optional(v.string()),
      sortBy: v.optional(v.string()),
      sortDirection: v.optional(v.string())
    }))
  },
  handler: async (ctx, args) => {
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Set up the query
    let query = ctx.db.query("projects");
    
    // Apply filters
    if (args.filter) {
      if (args.filter.name) {
        query = query.filter(q => 
          q.includes(q.field("name"), args.filter!.name!)
        );
      }
      
      if (args.filter.ownerId) {
        query = query.withIndex("by_owner", q => 
          q.eq("ownerId", args.filter!.ownerId!)
        );
      }
      
      if (args.filter.status) {
        query = query.withIndex("by_status", q => 
          q.eq("status", args.filter!.status!)
        );
      }
      
      if (args.filter.public !== undefined) {
        query = query.withIndex("by_public", q => 
          q.eq("public", args.filter!.public!)
        );
      }
      
      // If memberOf filter is applied, we handle it separately below
    }
    
    // Filter for projects accessible to the current user
    if (!userId) {
      // If no user is logged in, only show public projects
      query = query.withIndex("by_public", q => q.eq("public", true));
    } else if (!args.filter?.ownerId && !args.filter?.memberOf) {
      // If no specific owner or membership filter, include all accessible projects
      // This gets complex with document db, so we'll post-process below
    }
    
    // Apply sorting
    if (args.options?.sortBy) {
      const sortDirection = args.options.sortDirection === "desc" ? "desc" : "asc";
      
      switch (args.options.sortBy) {
        case "name":
          query = query.order("name", sortDirection);
          break;
        case "status":
          query = query.order("status", sortDirection);
          break;
        case "updatedAt":
          query = query.order("updatedAt", sortDirection);
          break;
        default:
          query = query.order("name", "asc"); // default sort
      }
    } else {
      // Default sort by name ascending
      query = query.order("name", "asc");
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
    
    // Fetch projects based on the current query
    let projects: Project[] = [];
    
    // If memberOf filter is applied, get projects directly
    if (args.filter?.memberOf) {
      projects = await getUserAccessibleProjects(ctx, args.filter.memberOf);
      
      // Apply any post-filtering for name or status
      if (args.filter.name) {
        projects = projects.filter(p => 
          p.name.toLowerCase().includes(args.filter!.name!.toLowerCase())
        );
      }
      
      if (args.filter.status) {
        projects = projects.filter(p => 
          p.status === args.filter!.status
        );
      }
      
      if (args.filter.public !== undefined) {
        projects = projects.filter(p => 
          p.public === args.filter!.public
        );
      }
    } else if (userId && !args.filter?.ownerId) {
      // If user is logged in and no specific owner filter, get all accessible projects
      projects = await getUserAccessibleProjects(ctx, userId);
      
      // Apply any post-filtering for name or status
      if (args.filter?.name) {
        projects = projects.filter(p => 
          p.name.toLowerCase().includes(args.filter.name!.toLowerCase())
        );
      }
      
      if (args.filter?.status) {
        projects = projects.filter(p => 
          p.status === args.filter.status
        );
      }
      
      if (args.filter?.public !== undefined) {
        projects = projects.filter(p => 
          p.public === args.filter.public
        );
      }
    } else {
      // Otherwise, execute the base query
      projects = await query.collect();
    }
    
    // Include team if requested
    if (args.filter?.includeTeam) {
      const projectsWithTeam = [];
      
      for (const project of projects) {
        projectsWithTeam.push(await expandProjectWithTeam(ctx, project));
      }
      
      return projectsWithTeam;
    }
    
    return projects;
  }
});

/**
 * Search projects by name or description
 */
export const searchProjects = query({
  args: {
    query: v.string(),
    limit: v.optional(v.number())
  },
  handler: async (ctx, args) => {
    // Get current user
    const user = await getCurrentUser(ctx);
    const userId = user?._id;
    
    // Set up the query
    let query = ctx.db.query("projects")
      .filter(q => 
        q.or(
          q.includes(q.field("name"), args.query),
          q.includes(q.field("description"), args.query)
        )
      );
    
    // Filter based on access control
    if (!userId) {
      // If no user is logged in, only show public projects
      query = query.filter(q => q.eq(q.field("public"), true));
    } else {
      // If user is logged in, show public projects and their own
      query = query.filter(q => 
        q.or(
          q.eq(q.field("public"), true),
          q.eq(q.field("ownerId"), userId)
        )
      );
      
      // For team membership, we need to post-process
    }
    
    // Apply limit
    const limit = args.limit || 20;
    query = query.take(limit);
    
    // Execute query
    let projects = await query.collect();
    
    // Post-process to include projects where user is a team member
    if (userId) {
      const teamMemberships = await ctx.db
        .query("teamMembers")
        .withIndex("by_user", q => q.eq("userId", userId))
        .collect();
      
      const teamProjectIds = new Set(teamMemberships.map(m => m.projectId.toString()));
      
      // Check if there are any team projects that matched the search but weren't included
      if (teamMemberships.length > 0) {
        const teamProjects = await ctx.db
          .query("projects")
          .filter(q => 
            q.and(
              q.or(
                q.includes(q.field("name"), args.query),
                q.includes(q.field("description"), args.query)
              ),
              q.eq(q.field("public"), false),
              q.neq(q.field("ownerId"), userId)
            )
          )
          .collect();
        
        for (const project of teamProjects) {
          if (teamProjectIds.has(project._id.toString())) {
            projects.push(project);
          }
        }
        
        // Limit the total results
        if (projects.length > limit) {
          projects = projects.slice(0, limit);
        }
      }
    }
    
    return projects;
  }
});

/**
 * Get user's recent projects
 */
export const getRecentProjects = query({
  args: {
    limit: v.optional(v.number())
  },
  handler: async (ctx, args) => {
    // Get current user
    const user = await getCurrentUser(ctx);
    
    if (!user) {
      return [];
    }
    
    const userId = user._id;
    const limit = args.limit || 5;
    
    // Get projects owned by user
    const ownedProjects = await ctx.db
      .query("projects")
      .withIndex("by_owner", q => q.eq("ownerId", userId))
      .order("updatedAt", "desc")
      .take(limit)
      .collect();
    
    // Get team memberships
    const teamMemberships = await ctx.db
      .query("teamMembers")
      .withIndex("by_user", q => q.eq("userId", userId))
      .collect();
    
    const memberProjectIds = new Set(teamMemberships.map(m => m.projectId.toString()));
    
    // Get projects user is a member of
    const memberProjects = [];
    
    for (const projectId of memberProjectIds) {
      const project = await ctx.db.get(projectId as Id<"projects">);
      
      if (project) {
        memberProjects.push(project);
      }
    }
    
    // Sort member projects by updatedAt
    memberProjects.sort((a, b) => b.updatedAt - a.updatedAt);
    
    // Combine and limit results
    const result = [...ownedProjects];
    const addedIds = new Set(ownedProjects.map(p => p._id.toString()));
    
    for (const project of memberProjects) {
      const projectId = project._id.toString();
      
      if (!addedIds.has(projectId)) {
        result.push(project);
        addedIds.add(projectId);
        
        // Stop if we've reached the limit
        if (result.length >= limit) {
          break;
        }
      }
    }
    
    return result;
  }
});
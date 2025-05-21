/**
 * Type definitions for users, teams, and projects
 */

import { Id } from "../_generated/dataModel";

/**
 * User document in the database
 */
export interface User {
  _id: Id<"users">;
  _creationTime: number;
  
  // Core fields
  email: string;
  name?: string;
  role: "admin" | "scientist" | "viewer";
  lastLogin?: number;
  
  // Metadata
  createdAt: number;
  updatedAt: number;
}

/**
 * Input for creating a new user
 */
export interface CreateUserInput {
  email: string;
  name?: string;
  role?: "admin" | "scientist" | "viewer";
}

/**
 * Input for updating a user
 */
export interface UpdateUserInput {
  name?: string;
  role?: "admin" | "scientist" | "viewer";
  lastLogin?: number;
}

/**
 * Project document in the database
 */
export interface Project {
  _id: Id<"projects">;
  _creationTime: number;
  
  // Core fields
  name: string;
  description?: string;
  ownerId: Id<"users">;
  status: "active" | "archived" | "completed";
  
  // Metadata
  createdAt: number;
  updatedAt: number;
  public: boolean;
}

/**
 * Input for creating a new project
 */
export interface CreateProjectInput {
  name: string;
  description?: string;
  status?: "active" | "archived" | "completed";
  public?: boolean;
}

/**
 * Input for updating a project
 */
export interface UpdateProjectInput {
  name?: string;
  description?: string;
  ownerId?: Id<"users">;
  status?: "active" | "archived" | "completed";
  public?: boolean;
}

/**
 * Team member document in the database
 */
export interface TeamMember {
  _id: Id<"teamMembers">;
  _creationTime: number;
  
  // Core fields
  projectId: Id<"projects">;
  userId: Id<"users">;
  role: "owner" | "editor" | "viewer";
  
  // Metadata
  createdAt: number;
  updatedAt: number;
}

/**
 * Input for creating a new team member
 */
export interface CreateTeamMemberInput {
  projectId: Id<"projects">;
  userId: Id<"users">;
  role?: "owner" | "editor" | "viewer";
}

/**
 * Input for updating a team member
 */
export interface UpdateTeamMemberInput {
  role?: "owner" | "editor" | "viewer";
}

/**
 * User with expanded project data
 */
export interface UserWithProjects extends User {
  projects: {
    owned: Project[];
    member: Array<{
      project: Project;
      role: string;
    }>;
  };
}

/**
 * Project with expanded team data
 */
export interface ProjectWithTeam extends Project {
  team: Array<{
    userId: Id<"users">;
    userEmail: string;
    userName?: string;
    role: string;
  }>;
}

/**
 * Filter conditions for querying users
 */
export interface UserFilter {
  email?: string;
  role?: "admin" | "scientist" | "viewer";
  includeProjects?: boolean;
}

/**
 * Filter conditions for querying projects
 */
export interface ProjectFilter {
  name?: string;
  ownerId?: Id<"users">;
  status?: "active" | "archived" | "completed";
  public?: boolean;
  includeTeam?: boolean;
  memberOf?: Id<"users">;
}

/**
 * Filter conditions for querying team members
 */
export interface TeamMemberFilter {
  projectId?: Id<"projects">;
  userId?: Id<"users">;
  role?: "owner" | "editor" | "viewer";
}

/**
 * Options for user queries
 */
export interface UserQueryOptions {
  limit?: number;
  cursor?: string;
  sortBy?: "email" | "role" | "createdAt";
  sortDirection?: "asc" | "desc";
}

/**
 * Options for project queries
 */
export interface ProjectQueryOptions {
  limit?: number;
  cursor?: string;
  sortBy?: "name" | "status" | "updatedAt";
  sortDirection?: "asc" | "desc";
}

/**
 * Options for team member queries
 */
export interface TeamMemberQueryOptions {
  limit?: number;
  cursor?: string;
}
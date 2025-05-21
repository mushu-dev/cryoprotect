/**
 * Validation functions for users, projects, and team members
 */

import { ConvexError } from "convex/values";
import { Id } from "../_generated/dataModel";
import { 
  CreateUserInput, 
  UpdateUserInput,
  CreateProjectInput,
  UpdateProjectInput,
  CreateTeamMemberInput,
  UpdateTeamMemberInput
} from "./types";

/**
 * Validates the input for creating a user
 */
export function validateCreateUserInput(input: CreateUserInput): void {
  // Email is required and should be a valid email
  if (!input.email || input.email.trim() === "" || !isValidEmail(input.email)) {
    throw new ConvexError("Valid email address is required");
  }
  
  // Name should be a string if provided
  if (input.name !== undefined && typeof input.name !== "string") {
    throw new ConvexError("Name must be a string");
  }
  
  // Role should be a valid value if provided
  if (input.role !== undefined) {
    const validRoles = ["admin", "scientist", "viewer"];
    if (!validRoles.includes(input.role)) {
      throw new ConvexError(`Role must be one of: ${validRoles.join(", ")}`);
    }
  }
}

/**
 * Validates the input for updating a user
 */
export function validateUpdateUserInput(input: UpdateUserInput): void {
  // Name should be a string if provided
  if (input.name !== undefined && typeof input.name !== "string") {
    throw new ConvexError("Name must be a string");
  }
  
  // Role should be a valid value if provided
  if (input.role !== undefined) {
    const validRoles = ["admin", "scientist", "viewer"];
    if (!validRoles.includes(input.role)) {
      throw new ConvexError(`Role must be one of: ${validRoles.join(", ")}`);
    }
  }
  
  // LastLogin should be a valid timestamp if provided
  if (input.lastLogin !== undefined && typeof input.lastLogin !== "number") {
    throw new ConvexError("Last login must be a valid timestamp");
  }
}

/**
 * Validates the input for creating a project
 */
export function validateCreateProjectInput(input: CreateProjectInput): void {
  // Name is required and should be a non-empty string
  if (!input.name || input.name.trim() === "") {
    throw new ConvexError("Project name is required");
  }
  
  // Description should be a string if provided
  if (input.description !== undefined && typeof input.description !== "string") {
    throw new ConvexError("Description must be a string");
  }
  
  // Status should be a valid value if provided
  if (input.status !== undefined) {
    const validStatuses = ["active", "archived", "completed"];
    if (!validStatuses.includes(input.status)) {
      throw new ConvexError(`Status must be one of: ${validStatuses.join(", ")}`);
    }
  }
  
  // Public should be a boolean if provided
  if (input.public !== undefined && typeof input.public !== "boolean") {
    throw new ConvexError("Public flag must be a boolean");
  }
}

/**
 * Validates the input for updating a project
 */
export function validateUpdateProjectInput(input: UpdateProjectInput): void {
  // Name should be a non-empty string if provided
  if (input.name !== undefined && (typeof input.name !== "string" || input.name.trim() === "")) {
    throw new ConvexError("Name must be a non-empty string");
  }
  
  // Description should be a string if provided
  if (input.description !== undefined && typeof input.description !== "string") {
    throw new ConvexError("Description must be a string");
  }
  
  // OwnerId should be a valid ID if provided
  if (input.ownerId !== undefined && typeof input.ownerId !== "object") {
    throw new ConvexError("Owner ID must be a valid user reference");
  }
  
  // Status should be a valid value if provided
  if (input.status !== undefined) {
    const validStatuses = ["active", "archived", "completed"];
    if (!validStatuses.includes(input.status)) {
      throw new ConvexError(`Status must be one of: ${validStatuses.join(", ")}`);
    }
  }
  
  // Public should be a boolean if provided
  if (input.public !== undefined && typeof input.public !== "boolean") {
    throw new ConvexError("Public flag must be a boolean");
  }
}

/**
 * Validates the input for creating a team member
 */
export function validateCreateTeamMemberInput(input: CreateTeamMemberInput): void {
  // ProjectId is required
  if (!input.projectId) {
    throw new ConvexError("Project ID is required");
  }
  
  // UserId is required
  if (!input.userId) {
    throw new ConvexError("User ID is required");
  }
  
  // Role should be a valid value if provided
  if (input.role !== undefined) {
    const validRoles = ["owner", "editor", "viewer"];
    if (!validRoles.includes(input.role)) {
      throw new ConvexError(`Role must be one of: ${validRoles.join(", ")}`);
    }
  }
}

/**
 * Validates the input for updating a team member
 */
export function validateUpdateTeamMemberInput(input: UpdateTeamMemberInput): void {
  // Role should be a valid value if provided
  if (input.role !== undefined) {
    const validRoles = ["owner", "editor", "viewer"];
    if (!validRoles.includes(input.role)) {
      throw new ConvexError(`Role must be one of: ${validRoles.join(", ")}`);
    }
  }
}

/**
 * Validates that a user ID exists
 */
export function validateUserExists(userId: Id<"users"> | undefined): void {
  if (!userId) {
    throw new ConvexError("User ID is required");
  }
}

/**
 * Validates that a project ID exists
 */
export function validateProjectExists(projectId: Id<"projects"> | undefined): void {
  if (!projectId) {
    throw new ConvexError("Project ID is required");
  }
}

/**
 * Validates that a team member ID exists
 */
export function validateTeamMemberExists(teamMemberId: Id<"teamMembers"> | undefined): void {
  if (!teamMemberId) {
    throw new ConvexError("Team member ID is required");
  }
}

/**
 * Validates that a user has access to a project
 */
export function validateProjectAccess(
  projectId: Id<"projects">,
  userId: Id<"users"> | null,
  isPublic: boolean,
  ownerId: Id<"users">,
  teamRoles?: Record<string, string>
): void {
  // If project is public, anyone can access
  if (isPublic) {
    return;
  }
  
  // If user is not logged in and project is not public, deny access
  if (!userId) {
    throw new ConvexError("Authentication required to access this project");
  }
  
  // If project belongs to user, allow access
  if (ownerId.equals(userId)) {
    return;
  }
  
  // If user is a team member, allow access
  if (teamRoles && userId.toString() in teamRoles) {
    return;
  }
  
  // Otherwise deny access
  throw new ConvexError("You do not have permission to access this project");
}

/**
 * Validates that a user has edit access to a project
 */
export function validateProjectEditAccess(
  projectId: Id<"projects">,
  userId: Id<"users"> | null,
  ownerId: Id<"users">,
  teamRoles?: Record<string, string>
): void {
  // If user is not logged in, deny access
  if (!userId) {
    throw new ConvexError("Authentication required to modify this project");
  }
  
  // If project belongs to user, allow access
  if (ownerId.equals(userId)) {
    return;
  }
  
  // If user is a team member with editor role, allow access
  if (teamRoles && 
      userId.toString() in teamRoles && 
      (teamRoles[userId.toString()] === "owner" || teamRoles[userId.toString()] === "editor")) {
    return;
  }
  
  // Otherwise deny access
  throw new ConvexError("You do not have permission to modify this project");
}

/**
 * Validates that a user has admin access to a project
 */
export function validateProjectAdminAccess(
  projectId: Id<"projects">,
  userId: Id<"users"> | null,
  ownerId: Id<"users">,
  teamRoles?: Record<string, string>
): void {
  // If user is not logged in, deny access
  if (!userId) {
    throw new ConvexError("Authentication required to administrate this project");
  }
  
  // If project belongs to user, allow access
  if (ownerId.equals(userId)) {
    return;
  }
  
  // If user is a team member with owner role, allow access
  if (teamRoles && 
      userId.toString() in teamRoles && 
      teamRoles[userId.toString()] === "owner") {
    return;
  }
  
  // Otherwise deny access
  throw new ConvexError("You do not have permission to administrate this project");
}

/**
 * Helper to validate email format
 */
function isValidEmail(email: string): boolean {
  const emailRegex = /^[^\s@]+@[^\s@]+\.[^\s@]+$/;
  return emailRegex.test(email);
}
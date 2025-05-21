/**
 * Export all user and project-related functionality
 */

// Export from users.ts
export { 
  createUser,
  getUser,
  updateUser,
  deleteUser,
  listUsers,
  searchUsers,
  getCurrentUserProfile,
  updateLastLogin
} from "./users";

// Export from projects.ts
export {
  createProject,
  getProject,
  updateProject,
  transferProjectOwnership,
  deleteProject,
  listProjects,
  searchProjects,
  getRecentProjects
} from "./projects";

// Export from teamMembers.ts
export {
  addTeamMember,
  getTeamMember,
  updateTeamMember,
  removeTeamMember,
  listTeamMembers,
  batchAddTeamMembers
} from "./teamMembers";

// Export types
export type * from "./types";

// Export helpers
export {
  formatProjectStatus,
  formatUserRole,
  formatTeamRole,
  getProjectTeamRoles,
  expandUserWithProjects,
  expandProjectWithTeam,
  ensureTeamMember,
  isUserAdmin,
  getUserAccessibleProjects
} from "./helpers";
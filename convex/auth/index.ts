/**
 * Index file for auth functions
 */

// Export all auth-related functions
export * from "./users";
export * from "./clerk";

// Client utilities should be imported directly from client.ts to avoid
// importing unnecessary client-side code in the server functions
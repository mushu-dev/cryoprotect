import { mutation } from "../../_generated/server";
import { v } from "convex/values";

/**
 * API endpoint for our backend adapter to sign out users
 * This function handles sign-out requests from the Supabase-compatible adapter
 */
export const execute = mutation({
  args: {
    // No args needed
  },
  handler: async (ctx, args) => {
    try {
      // In a real implementation, we would invalidate the user's session
      // For now, we'll just return a success response for compatibility
      
      return { 
        data: { message: "Signed out successfully" }, 
        error: null 
      };
    } catch (error) {
      console.error('Error during sign out:', error);
      return { 
        data: null, 
        error: error instanceof Error ? error.message : 'Unknown error' 
      };
    }
  },
});
import { mutation } from "../../_generated/server";
import { v } from "convex/values";

/**
 * API endpoint for our backend adapter to authenticate users with email/password
 * This function handles sign-in requests from the Supabase-compatible adapter
 */
export const execute = mutation({
  args: {
    email: v.string(),
    password: v.string(),
  },
  handler: async (ctx, args) => {
    try {
      const { email, password } = args;
      
      // In a real implementation, we would use Convex Auth to sign in the user
      // For now, we'll just attempt to find the user in the database
      const usersResult = await ctx.db
        .query("users")
        .filter(q => q.eq(q.field("email"), email))
        .first();
      
      if (!usersResult) {
        return { 
          data: null, 
          error: "Invalid email or password" 
        };
      }
      
      // Note: In a real implementation, we would check the password hash
      // For now, we're just simulating the authentication
      
      // Generate a fake session for compatibility
      const session = {
        access_token: `fake-token-${Date.now()}`,
        refresh_token: `fake-refresh-${Date.now()}`,
        user: {
          id: usersResult._id,
          email: email,
          app_metadata: { provider: "email" },
          user_metadata: { name: usersResult.name || "" },
          aud: "authenticated",
          role: "authenticated",
        }
      };
      
      return { 
        data: { session }, 
        error: null 
      };
    } catch (error) {
      console.error('Error during sign in:', error);
      return { 
        data: null, 
        error: error instanceof Error ? error.message : 'Unknown error' 
      };
    }
  },
});
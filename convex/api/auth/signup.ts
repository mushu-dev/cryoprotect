import { mutation } from "../../_generated/server";
import { v } from "convex/values";

/**
 * API endpoint for our backend adapter to register new users
 * This function handles sign-up requests from the Supabase-compatible adapter
 */
export const execute = mutation({
  args: {
    email: v.string(),
    password: v.string(),
    name: v.optional(v.string()),
  },
  handler: async (ctx, args) => {
    try {
      const { email, password, name } = args;
      
      // Check if the user already exists
      const existingUser = await ctx.db
        .query("users")
        .filter(q => q.eq(q.field("email"), email))
        .first();
      
      if (existingUser) {
        return { 
          data: null, 
          error: "User with this email already exists" 
        };
      }
      
      // Create a new user
      // Note: In a real implementation, we would hash the password
      const timestamp = Date.now();
      const userId = await ctx.db.insert("users", {
        email,
        // We would normally store a hashed password
        // This is just a placeholder for the adapter compatibility
        passwordHash: `placeholder_${password.length}`,
        name: name || "",
        createdAt: timestamp,
        updatedAt: timestamp,
        role: "user",
      });
      
      // Generate a fake session for compatibility
      const session = {
        access_token: `fake-token-${Date.now()}`,
        refresh_token: `fake-refresh-${Date.now()}`,
        user: {
          id: userId,
          email: email,
          app_metadata: { provider: "email" },
          user_metadata: { name: name || "" },
          aud: "authenticated",
          role: "authenticated",
        }
      };
      
      return { 
        data: { 
          user: {
            id: userId,
            email,
            name: name || "",
          },
          session
        }, 
        error: null 
      };
    } catch (error) {
      console.error('Error during sign up:', error);
      return { 
        data: null, 
        error: error instanceof Error ? error.message : 'Unknown error' 
      };
    }
  },
});
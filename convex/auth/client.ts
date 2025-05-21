/**
 * Client-side auth utilities for Convex with Clerk
 * 
 * NOTE: This file is meant to be used as a reference for frontend implementation.
 * It should NOT be imported directly by Convex functions as it contains
 * client-side dependencies that are not available in the Convex runtime.
 * 
 * To use these utilities, copy this code into your frontend codebase.
 */

/*
// --- FRONTEND CODE REFERENCE ---
// Import these in your frontend application:
// 
// import { ConvexReactClient } from "convex/react";
// import { ClerkProvider, useAuth } from "@clerk/clerk-react";
// import { ConvexProviderWithClerk } from "convex/react-clerk";

/**
 * Initialize Convex with Clerk authentication
 * 
 * Usage in your app entry point:
 * ```
 * const convex = createConvexClientWithAuth();
 * 
 * function App() {
 *   return (
 *     <ConvexProviderWithClerk client={convex} useAuth={useAuth}>
 *       <ClerkProvider publishableKey={process.env.NEXT_PUBLIC_CLERK_PUBLISHABLE_KEY}>
 *         <YourApp />
 *       </ClerkProvider>
 *     </ConvexProviderWithClerk>
 *   );
 * }
 * ```
 */
/*
function createConvexClientWithAuth(convexUrl) {
  const url = convexUrl || process.env.NEXT_PUBLIC_CONVEX_URL || "";
  return new ConvexReactClient(url);
}
*/

/**
 * Get the current auth status and user information
 * 
 * Usage in a component:
 * ```
 * function Profile() {
 *   const { isAuthenticated, user, isLoading } = useConvexAuth();
 *   
 *   if (isLoading) return <div>Loading...</div>;
 *   
 *   if (!isAuthenticated) {
 *     return <div>Please log in</div>;
 *   }
 *   
 *   return <div>Welcome, {user.name}</div>;
 * }
 * ```
 */
/*
function useConvexAuth() {
  const { isSignedIn, isLoaded, user } = useAuth();
  
  return {
    isAuthenticated: isSignedIn,
    isLoading: !isLoaded,
    user: user ? {
      id: user.id,
      email: user.primaryEmailAddress?.emailAddress,
      name: user.fullName,
      imageUrl: user.imageUrl,
    } : null,
  };
}
*/
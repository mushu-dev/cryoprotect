/**
 * Convex client configuration for the browser
 * 
 * This file sets up the Convex client for use in the browser,
 * enabling real-time data synchronization and reactive queries.
 */
import { ConvexReactClient } from "convex/react";

// Create a client using the Convex deployment URL from the environment variable
// If not defined (during development), use a default value
export const convex = new ConvexReactClient(
  process.env.NEXT_PUBLIC_CONVEX_URL || "https://dynamic-mink-63.convex.cloud"
);
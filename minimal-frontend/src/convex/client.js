// Convex client setup
import { ConvexReactClient } from 'convex/react';

/**
 * ConvexReactClient instance for interacting with Convex from React components.
 * 
 * This client handles all communication with the Convex backend, including:
 * - Fetching data with queries
 * - Modifying data with mutations
 * - Real-time subscriptions to data changes
 * - Authentication
 */
export const convex = new ConvexReactClient(
  process.env.NEXT_PUBLIC_CONVEX_URL || "https://primary-meerkat-478.convex.cloud"
);

/**
 * Health check method to verify if the client is properly configured and connected.
 * 
 * @returns {Promise<boolean>} - True if the client is connected, false otherwise
 */
convex.health = async function() {
  try {
    // Try to make a basic call to Convex to verify connection
    // This will fail if the client isn't properly initialized
    await this.getAuth();
    return true;
  } catch (error) {
    console.error("Convex health check failed:", error);
    return false;
  }
};

/**
 * Get all environment variables related to Convex
 * @returns {Object} - Object containing all Convex-related environment variables
 */
convex.getEnvironment = function() {
  return {
    url: process.env.NEXT_PUBLIC_CONVEX_URL || "https://primary-meerkat-478.convex.cloud",
    useConvex: process.env.NEXT_PUBLIC_USE_CONVEX === 'true' || false,
    environment: process.env.NODE_ENV || 'development'
  };
};
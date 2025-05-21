/**
 * Unified Convex client for CryoProtect.
 * 
 * This client can be used by both the main and minimal frontends,
 * providing a consistent interface to Convex functionality.
 */
import { ConvexReactClient } from "convex/react";

// Determine if Convex is enabled based on environment variable
const isConvexEnabled = process.env.NEXT_PUBLIC_USE_CONVEX === 'true';

// Get Convex URL from environment variable with fallback to production
const convexUrl = process.env.NEXT_PUBLIC_CONVEX_URL || "https://upbeat-parrot-866.convex.cloud";

// Create and export the Convex client instance with production configuration
export const convex = new ConvexReactClient(convexUrl, {
  // Production configuration options
  unsavedChangesWarning: false, // Disable unsaved changes warning
  networkTimeout: 10000, // 10 seconds timeout for network requests
  initialPresence: {}, // Start with empty presence state
});

// Add health check method for diagnostics
convex.health = async function() {
  try {
    await this.getAuth();
    return true;
  } catch (error) {
    console.error("Convex health check failed:", error);
    return false;
  }
};

// Export convenience function to check if Convex is enabled
export const isEnabled = () => isConvexEnabled;

// Export connection URL for diagnostics
export const getConnectionUrl = () => convexUrl;
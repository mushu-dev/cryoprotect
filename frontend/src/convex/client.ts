import { ConvexReactClient } from "convex/react";

// Determine if Convex is enabled based on environment variable
const isConvexEnabled = process.env.NEXT_PUBLIC_USE_CONVEX === 'true';

// Get Convex URL from environment variable with fallback
const convexUrl = process.env.NEXT_PUBLIC_CONVEX_URL || "https://upbeat-parrot-866.convex.cloud";

// Create and export the Convex client instance with production configuration
export const convex = new ConvexReactClient(convexUrl, {
  // Production configuration options
  unsavedChangesWarning: false, // Disable unsaved changes warning
  networkTimeout: 10000, // 10 seconds timeout for network requests
  initialPresence: {}, // Start with empty presence state
});

// Log Convex initialization in development
if (process.env.NODE_ENV === 'development') {
  console.log(`Convex client initialized with URL: ${convexUrl}`);
  console.log(`Convex is ${isConvexEnabled ? 'enabled' : 'disabled'}`);
}
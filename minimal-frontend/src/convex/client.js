/**
 * Convex client setup
 * Enhanced with better health check, resilience, and configuration
 */
import { ConvexReactClient } from 'convex/react';

// Get Convex URL from environment variables
const CONVEX_URL = process.env.NEXT_PUBLIC_CONVEX_URL || "https://upbeat-parrot-866.convex.cloud";

/**
 * ConvexReactClient instance for interacting with Convex from React components.
 * 
 * This client handles all communication with the Convex backend, including:
 * - Fetching data with queries
 * - Modifying data with mutations
 * - Real-time subscriptions to data changes
 * - Authentication
 */
export const convex = new ConvexReactClient(CONVEX_URL);

// Error tracking for connection issues
let connectionErrorCount = 0;
const MAX_CONNECTION_ERRORS = 5;
let lastConnectionAttempt = 0;
const CONNECTION_RETRY_DELAY = 15000; // 15 seconds

/**
 * Health check method to verify if the client is properly configured and connected.
 * Includes retry mechanism and error tracking for improved reliability.
 * 
 * @returns {Promise<Object>} - Health status object
 */
convex.health = async function() {
  const now = Date.now();
  
  // If we've had too many errors, wait before trying again
  if (connectionErrorCount >= MAX_CONNECTION_ERRORS && 
      now - lastConnectionAttempt < CONNECTION_RETRY_DELAY) {
    return { 
      connected: false, 
      status: 'rate_limited',
      message: `Too many connection failures. Retry after ${Math.ceil((CONNECTION_RETRY_DELAY - (now - lastConnectionAttempt)) / 1000)}s`,
      errors: connectionErrorCount
    };
  }
  
  lastConnectionAttempt = now;
  
  try {
    // Try to make a basic call to Convex to verify connection
    await this.getAuth();
    
    // Reset error count on success
    connectionErrorCount = 0;
    
    return { 
      connected: true, 
      status: 'connected',
      message: 'Convex connection established successfully' 
    };
  } catch (error) {
    // Increment error count on failure
    connectionErrorCount++;
    
    console.error("Convex health check failed:", error);
    
    return { 
      connected: false, 
      status: 'error',
      message: error.message || 'Failed to connect to Convex',
      errors: connectionErrorCount
    };
  }
};

/**
 * Get all environment variables related to Convex
 * @returns {Object} - Object containing all Convex-related environment variables
 */
convex.getEnvironment = function() {
  const isClient = typeof window !== 'undefined';
  return {
    url: CONVEX_URL,
    useConvex: process.env.NEXT_PUBLIC_USE_CONVEX === 'true' || false,
    environment: process.env.NODE_ENV || 'development',
    isClient
  };
};

/**
 * Get connection status for Convex
 * @returns {Object} - Connection status object
 */
convex.getConnectionStatus = function() {
  return {
    url: CONVEX_URL,
    enabled: process.env.NEXT_PUBLIC_USE_CONVEX === 'true',
    errorCount: connectionErrorCount,
    lastAttempt: lastConnectionAttempt
  };
};

// Export the convex client
export default convex;
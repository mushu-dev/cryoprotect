/**
 * Convex Client Provider
 * 
 * This component provides the Convex client context to the application.
 * It should be used at the root level of the application to provide
 * Convex client access to all components that need it.
 */
import React from 'react';
import { ConvexProvider, ConvexReactClient } from 'convex/react';

// Create a Convex client using the Convex URL from environment variables
const convexUrl = process.env.NEXT_PUBLIC_CONVEX_URL || 'https://reserved-catfish-123.convex.cloud';
const convex = new ConvexReactClient(convexUrl);

/**
 * Provides Convex client context to the application
 * 
 * @param {Object} props - Component props
 * @param {React.ReactNode} props.children - Child components to render
 * @returns {React.ReactElement} Rendered component
 */
export default function ConvexClientProvider({ children }) {
  return (
    <ConvexProvider client={convex}>
      {children}
    </ConvexProvider>
  );
}
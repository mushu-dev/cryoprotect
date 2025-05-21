/**
 * Conditional Convex Provider
 * 
 * This component conditionally wraps the application with the Convex provider
 * based on whether Convex is enabled in the environment.
 */
import React, { useEffect, useState } from 'react';
import dynamic from 'next/dynamic';

// Dynamically import the ConvexClientProvider to avoid server-side rendering issues
const ConvexClientProvider = dynamic(
  () => import('./ConvexClientProvider'),
  { ssr: false }
);

/**
 * Conditionally wraps children with the Convex provider if Convex is enabled
 * 
 * @param {Object} props - Component props
 * @param {React.ReactNode} props.children - Child components to render
 * @returns {React.ReactElement} Rendered component
 */
export default function ConditionalConvexProvider({ children }) {
  const [isConvexEnabled, setIsConvexEnabled] = useState(false);
  const [isClient, setIsClient] = useState(false);

  useEffect(() => {
    setIsClient(true);
    setIsConvexEnabled(process.env.NEXT_PUBLIC_USE_CONVEX === 'true');
  }, []);

  // Don't render anything on the server
  if (!isClient) {
    return null;
  }

  // Only wrap with the Convex provider if Convex is enabled
  if (isConvexEnabled) {
    return <ConvexClientProvider>{children}</ConvexClientProvider>;
  }

  // Otherwise, just render the children
  return <>{children}</>;
}
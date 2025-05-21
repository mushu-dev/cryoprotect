import { ConvexProvider } from 'convex/react';
import { convex } from './client';

/**
 * ConvexClientProvider component
 * 
 * This component wraps the application with the Convex client provider.
 * It's used when NEXT_PUBLIC_USE_CONVEX=true to enable direct Convex integration.
 * 
 * @param {Object} props - Component props
 * @param {React.ReactNode} props.children - Child components
 */
export function ConvexClientProvider({ children }) {
  // Render the children with the Convex provider
  return (
    <ConvexProvider client={convex}>
      {children}
    </ConvexProvider>
  );
}
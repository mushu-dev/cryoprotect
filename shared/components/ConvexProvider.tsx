/**
 * Unified ConvexProvider component for CryoProtect.
 * 
 * This component can be used by both the main and minimal frontends to
 * conditionally wrap the application with ConvexProvider when enabled.
 */
import React from 'react';
import { ConvexProvider } from "convex/react";
import { convex, isEnabled } from '../convex/client';

type ConvexWrapperProps = {
  children: React.ReactNode;
};

/**
 * Conditionally wraps children with ConvexProvider when Convex is enabled.
 * When Convex is disabled, renders children directly without wrapping.
 */
export function ConvexWrapper({ children }: ConvexWrapperProps) {
  // Only wrap with ConvexProvider if Convex is enabled
  if (!isEnabled()) {
    return <>{children}</>;
  }

  return (
    <ConvexProvider client={convex}>
      {children}
    </ConvexProvider>
  );
}

/**
 * React hook to check if Convex is enabled.
 * Can be used in components to conditionally render Convex-specific features.
 */
export function useConvexStatus() {
  return {
    enabled: isEnabled()
  };
}
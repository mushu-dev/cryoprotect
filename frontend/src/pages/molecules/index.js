import React from 'react';
import dynamic from 'next/dynamic';

// Import the standard API version
import DefaultMoleculesPage from './index.api';

// Dynamically import the Convex version to avoid build issues
const ConvexMoleculesPage = dynamic(
  () => import('./index.convex'),
  { 
    ssr: false,
    loading: () => <div className="container mx-auto px-4 py-8 text-center">Loading Convex components...</div>
  }
);

/**
 * Router component that selects between Convex and API implementations
 * based on the NEXT_PUBLIC_USE_CONVEX environment variable
 */
export default function MoleculesPage() {
  // Check if Convex is enabled
  const isConvexEnabled = process.env.NEXT_PUBLIC_USE_CONVEX === 'true';
  
  // Return the appropriate implementation based on the environment setting
  return isConvexEnabled ? <ConvexMoleculesPage /> : <DefaultMoleculesPage />;
}
import React from 'react';
import Head from 'next/head';
import dynamic from 'next/dynamic';

// Dynamically import Convex provider to avoid issues during build
const ConvexClientProvider = dynamic(
  () => import('../convex/ConvexClientProvider').then(mod => mod.ConvexClientProvider),
  { ssr: false }
);

function MyApp({ Component, pageProps }) {
  // Check if Convex is enabled from environment
  // This will be 'true' when running with NEXT_PUBLIC_USE_CONVEX=true
  const isConvexEnabled = process.env.NEXT_PUBLIC_USE_CONVEX === 'true';
  
  // If Convex is enabled, wrap in Convex provider
  return (
    <>
      <Head>
        <meta name="viewport" content="width=device-width, initial-scale=1" />
      </Head>
      {isConvexEnabled ? (
        <ConvexClientProvider>
          <Component {...pageProps} />
        </ConvexClientProvider>
      ) : (
        <Component {...pageProps} />
      )}
    </>
  );
}

export default MyApp;

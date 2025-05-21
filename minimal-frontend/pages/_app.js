import React from 'react';
import '../styles/globals.css';
import { ConvexClientProvider } from '../src/convex/ConvexClientProvider';

// Check if Convex should be enabled
const useConvex = process.env.NEXT_PUBLIC_USE_CONVEX === 'true';

function MyApp({ Component, pageProps }) {
  // If Convex is enabled, wrap the app with ConvexClientProvider
  if (useConvex) {
    return (
      <ConvexClientProvider>
        <Component {...pageProps} />
      </ConvexClientProvider>
    );
  }
  
  // Otherwise, just render the component
  return <Component {...pageProps} />;
}

export default MyApp;
import React from 'react';
import '../styles/globals.css';

// Disable Convex for minimal deployment
function MyApp({ Component, pageProps }) {
  return <Component {...pageProps} />;
}

export default MyApp;
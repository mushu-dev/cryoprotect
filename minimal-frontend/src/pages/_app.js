/**
 * Custom App Component
 * Provides global configuration and layout for all pages
 */
import React, { useEffect } from 'react';
import Head from 'next/head';
import { QueryClient, QueryClientProvider } from '@tanstack/react-query';
import { CircuitProvider } from '../components/circuit-breaker';
import '../styles/globals.css'; // Make sure this exists or create it

// Determine if Convex is enabled
const isConvexEnabled = () => {
  if (typeof process !== 'undefined') {
    return process.env.NEXT_PUBLIC_USE_CONVEX === 'true';
  }
  
  if (typeof window !== 'undefined') {
    return window.localStorage.getItem('use_convex') === 'true';
  }
  
  return false;
};

// Create a client
const queryClient = new QueryClient({
  defaultOptions: {
    queries: {
      refetchOnWindowFocus: false,
      retry: 1,
    },
  },
});

export default function App({ Component, pageProps }) {
  // Initialize Convex if enabled
  useEffect(() => {
    const initServices = async () => {
      if (isConvexEnabled()) {
        try {
          // Initialize molecule service
          const moleculeService = await import('../services/molecules/molecule-service');
          
          // Initialize mixture service
          const mixtureService = await import('../services/mixtures/mixture-service');
          
          // Assuming window.ConvexClient would be defined elsewhere if Convex is enabled
          if (typeof window !== 'undefined' && window.ConvexClient) {
            const convexClient = window.ConvexClient;
            moleculeService.initConvex(convexClient.molecules);
            mixtureService.initConvex(convexClient.mixtures);
            console.log('Convex services initialized');
          }
        } catch (error) {
          console.error('Error initializing Convex services:', error);
        }
      }
    };
    
    initServices();
  }, []);
  
  return (
    <>
      <Head>
        <title>CryoProtect</title>
        <meta name="description" content="A platform for exploring and designing cryoprotectants" />
        <meta name="viewport" content="width=device-width, initial-scale=1" />
        <link rel="icon" href="/favicon.ico" />
      </Head>
      
      <QueryClientProvider client={queryClient}>
        <CircuitProvider>
          {/* Add a simple navigation bar */}
          <nav className="bg-primary/10 border-b border-primary/20 py-4">
            <div className="container mx-auto flex justify-between items-center">
              <a href="/" className="text-2xl font-bold">CryoProtect</a>
              <div className="space-x-6">
                <a href="/molecules" className="hover:text-primary">Molecules</a>
                <a href="/mixtures" className="hover:text-primary">Mixtures</a>
              </div>
            </div>
          </nav>
          
          <Component {...pageProps} />
        </CircuitProvider>
      </QueryClientProvider>
    </>
  );
}
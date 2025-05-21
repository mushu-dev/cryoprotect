/**
 * Mixtures Page
 * Main page for displaying mixture list with search and filter capabilities
 */
import React from 'react';
import { MixturesList } from '../../components/mixtures';
import { QueryClient, QueryClientProvider } from '@tanstack/react-query';

// Create a client
const queryClient = new QueryClient({
  defaultOptions: {
    queries: {
      refetchOnWindowFocus: false,
      retry: 1,
    },
  },
});

export default function MixturesPage() {
  return (
    <QueryClientProvider client={queryClient}>
      <div className="container mx-auto py-8">
        <h1 className="text-3xl font-bold mb-6">Mixtures</h1>
        <p className="mb-8 text-muted-foreground">
          Browse our database of mixtures, including cryoprotectant formulations. 
          Use the search and filters to find specific mixtures.
        </p>
        
        <MixturesList />
      </div>
    </QueryClientProvider>
  );
}
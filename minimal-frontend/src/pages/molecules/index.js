/**
 * Molecules Page
 * Main page for displaying molecule list with search and filter capabilities
 */
import React from 'react';
import { MoleculesList } from '../../components/molecules';
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

export default function MoleculesPage() {
  return (
    <QueryClientProvider client={queryClient}>
      <div className="container mx-auto py-8">
        <h1 className="text-3xl font-bold mb-6">Molecules</h1>
        <p className="mb-8 text-muted-foreground">
          Browse our database of molecules, including cryoprotectants and related compounds. 
          Use the search and filters to find specific molecules.
        </p>
        
        <MoleculesList />
      </div>
    </QueryClientProvider>
  );
}
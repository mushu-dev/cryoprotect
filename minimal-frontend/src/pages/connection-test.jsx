import React, { useState, useEffect } from 'react';
import Head from 'next/head';
import ConnectionVerificationPanel from '../components/ConnectionVerificationPanel';

// Mock implementation of Convex client if needed
const createMockConvexClient = () => {
  return {
    query: (name) => {
      // Mock query function
      return async () => {
        if (name === 'system.healthCheck') {
          return { status: 'ok', timestamp: new Date().toISOString() };
        }
        if (name === 'molecules.getFirstMolecule') {
          return { id: 'mock-molecule-1', name: 'Mock Molecule' };
        }
        throw new Error(`Unknown query: ${name}`);
      };
    },
    table: (tableName) => {
      // Mock table API
      return {
        select: () => ({
          limit: () => ({
            execute: async () => ({ data: [{ id: 'mock-data' }] })
          })
        })
      };
    }
  };
};

export default function ConnectionTestPage() {
  const [convexClient, setConvexClient] = useState(null);
  const [useConvex, setUseConvex] = useState(false);
  const [useMock, setUseMock] = useState(false);
  const [authToken, setAuthToken] = useState(null);

  // Effect to initialize clients based on settings
  useEffect(() => {
    // Clean up function
    let cleanupFn = () => {};
    
    const initClients = async () => {
      if (useConvex) {
        if (useMock) {
          // Use mock client
          setConvexClient(createMockConvexClient());
        } else {
          try {
            // Try to dynamically import Convex client - this assumes it exists in the project
            // You may need to adjust this based on actual client implementation
            const importedModule = await import('../utils/convex-client');
            if (importedModule && importedModule.default) {
              const client = importedModule.default.createClient();
              setConvexClient(client);
              cleanupFn = () => {
                if (client && client.close) {
                  client.close();
                }
              };
            } else {
              console.error('Convex client module exists but has no default export');
              // Fall back to mock client
              setConvexClient(createMockConvexClient());
            }
          } catch (error) {
            console.error('Error importing Convex client:', error);
            // Fall back to mock client
            setConvexClient(createMockConvexClient());
          }
        }
      } else {
        setConvexClient(null);
      }
    };
    
    initClients();
    
    // Cleanup function
    return cleanupFn;
  }, [useConvex, useMock]);

  return (
    <div className="min-h-screen bg-gray-50 py-8">
      <Head>
        <title>Connection Test Dashboard</title>
        <meta name="description" content="Test connections between frontend and backend systems" />
      </Head>
      
      <div className="max-w-4xl mx-auto">
        <div className="bg-white shadow-md rounded-lg p-6 mb-8">
          <h1 className="text-2xl font-bold text-gray-900 mb-4">Connection Test Dashboard</h1>
          <p className="text-gray-600 mb-6">
            This dashboard allows you to verify connectivity between frontend and backend systems,
            including API, database, and authentication services.
          </p>
          
          <div className="mb-8 p-4 bg-gray-100 rounded-md flex flex-wrap gap-4">
            <div>
              <label className="flex items-center">
                <input
                  type="checkbox"
                  checked={useConvex}
                  onChange={(e) => setUseConvex(e.target.checked)}
                  className="mr-2 h-4 w-4"
                />
                <span>Enable Convex Client</span>
              </label>
            </div>
            
            {useConvex && (
              <div>
                <label className="flex items-center">
                  <input
                    type="checkbox"
                    checked={useMock}
                    onChange={(e) => setUseMock(e.target.checked)}
                    className="mr-2 h-4 w-4"
                  />
                  <span>Use Mock Convex Client</span>
                </label>
              </div>
            )}
            
            <div className="ml-auto">
              <button
                onClick={() => window.location.reload()}
                className="px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700"
              >
                Reload Page
              </button>
            </div>
          </div>
          
          <ConnectionVerificationPanel 
            convexClient={convexClient}
            authToken={authToken}
            autoCheck={true}
          />
        </div>
        
        <div className="bg-white shadow-md rounded-lg p-6">
          <h2 className="text-xl font-semibold mb-4">Connection Status Legend</h2>
          <div className="space-y-4">
            <div className="flex items-center">
              <div className="w-24 px-3 py-1 bg-green-100 text-green-800 border-green-200 rounded-full text-sm font-medium mr-4">
                <span className="mr-1">✓</span> Connected
              </div>
              <p>The connection was successfully established and verified</p>
            </div>
            
            <div className="flex items-center">
              <div className="w-24 px-3 py-1 bg-red-100 text-red-800 border-red-200 rounded-full text-sm font-medium mr-4">
                <span className="mr-1">✗</span> Failed
              </div>
              <p>The connection attempt failed</p>
            </div>
            
            <div className="flex items-center">
              <div className="w-24 px-3 py-1 bg-yellow-100 text-yellow-800 border-yellow-200 rounded-full text-sm font-medium mr-4">
                <span className="mr-1">⟳</span> Checking...
              </div>
              <p>Connection verification is in progress</p>
            </div>
            
            <div className="flex items-center">
              <div className="w-24 px-3 py-1 bg-gray-100 text-gray-800 border-gray-200 rounded-full text-sm font-medium mr-4">
                <span className="mr-1">?</span> Unknown
              </div>
              <p>Connection status has not been determined</p>
            </div>
          </div>
          
          <div className="mt-6 pt-6 border-t">
            <h3 className="font-semibold mb-2">Troubleshooting Tips</h3>
            <ul className="list-disc pl-5 space-y-2">
              <li>If API connections fail, check that the backend server is running and accessible</li>
              <li>For Convex connection failures, verify that environment variables are properly configured</li>
              <li>Authentication issues might be related to expired or invalid tokens</li>
              <li>Database connection problems often indicate configuration or access permission issues</li>
              <li>Service role connection failures usually mean the service role key is missing or invalid</li>
            </ul>
          </div>
        </div>
      </div>
    </div>
  );
}
import React, { useEffect, useState } from 'react';
import apiClient from '../services/api-client';

interface ApiStatus {
  connected: boolean;
  loading: boolean;
  error: string | null;
  details: any | null;
}

const ApiConnectivityTest: React.FC = () => {
  const [status, setStatus] = useState<ApiStatus>({
    connected: false,
    loading: true,
    error: null,
    details: null,
  });

  useEffect(() => {
    const checkConnectivity = async () => {
      try {
        setStatus(prev => ({ ...prev, loading: true, error: null }));
        // Try the API connectivity endpoint with proper path
        const response = await apiClient.get('/api/v1/health/connectivity');
        
        setStatus({
          connected: response.data.status === 'connected',
          loading: false,
          error: null,
          details: response.data,
        });
      } catch (error: any) {
        console.error('API connectivity check failed:', error);
        
        // Try the fallback endpoint for netlify deployment
        try {
          const netlifyApiUrl = process.env.NEXT_PUBLIC_API_URL || 
            'https://cryoprotect-8030e4025428.herokuapp.com/v1';
          
          // Make direct fetch request as fallback
          const directResponse = await fetch(`${netlifyApiUrl}/health/connectivity`);
          
          if (directResponse.ok) {
            const data = await directResponse.json();
            setStatus({
              connected: data.status === 'connected',
              loading: false,
              error: null,
              details: data,
            });
          } else {
            throw new Error(`HTTP error! status: ${directResponse.status}`);
          }
        } catch (fallbackError: any) {
          console.error('API fallback connectivity check failed:', fallbackError);
          
          setStatus({
            connected: false,
            loading: false,
            error: error.message || 'Failed to connect to API',
            details: null,
          });
        }
      }
    };

    checkConnectivity();
  }, []);

  return (
    <div className="p-4 max-w-md mx-auto bg-white rounded-lg border border-gray-200 shadow-md">
      <h2 className="text-xl font-semibold mb-4">API Connectivity Check</h2>
      
      {status.loading ? (
        <div className="flex items-center space-x-2">
          <div className="animate-spin rounded-full h-4 w-4 border-b-2 border-blue-500"></div>
          <span>Checking connection...</span>
        </div>
      ) : status.connected ? (
        <div className="bg-green-100 border-l-4 border-green-500 p-4">
          <div className="flex items-start">
            <div className="flex-shrink-0">
              <svg className="h-5 w-5 text-green-500" fill="currentColor" viewBox="0 0 20 20">
                <path fillRule="evenodd" d="M10 18a8 8 0 100-16 8 8 0 000 16zm3.707-9.293a1 1 0 00-1.414-1.414L9 10.586 7.707 9.293a1 1 0 00-1.414 1.414l2 2a1 1 0 001.414 0l4-4z" clipRule="evenodd" />
              </svg>
            </div>
            <div className="ml-3">
              <p className="text-sm font-medium text-green-800">
                Connected to API successfully
              </p>
            </div>
          </div>
        </div>
      ) : (
        <div className="bg-red-100 border-l-4 border-red-500 p-4">
          <div className="flex items-start">
            <div className="flex-shrink-0">
              <svg className="h-5 w-5 text-red-500" fill="currentColor" viewBox="0 0 20 20">
                <path fillRule="evenodd" d="M10 18a8 8 0 100-16 8 8 0 000 16zM8.707 7.293a1 1 0 00-1.414 1.414L8.586 10l-1.293 1.293a1 1 0 101.414 1.414L10 11.414l1.293 1.293a1 1 0 001.414-1.414L11.414 10l1.293-1.293a1 1 0 00-1.414-1.414L10 8.586 8.707 7.293z" clipRule="evenodd" />
              </svg>
            </div>
            <div className="ml-3">
              <p className="text-sm font-medium text-red-800">
                Failed to connect to API: {status.error}
              </p>
            </div>
          </div>
        </div>
      )}

      {status.details && (
        <div className="mt-4">
          <h3 className="text-md font-semibold mb-2">API Details:</h3>
          <pre className="bg-gray-100 p-2 text-xs rounded overflow-auto max-h-56">
            {JSON.stringify(status.details, null, 2)}
          </pre>
        </div>
      )}
      
      <div className="mt-4 pt-4 border-t border-gray-200">
        <h3 className="text-sm font-semibold mb-2">Environment:</h3>
        <div className="text-xs text-gray-700">
          <p><strong>API URL:</strong> {process.env.NEXT_PUBLIC_API_URL || 'Not configured'}</p>
          <p><strong>Environment:</strong> {process.env.NEXT_PUBLIC_ENVIRONMENT || 'development'}</p>
          <p><strong>Deployment:</strong> {process.env.NEXT_PUBLIC_VERCEL_ENV || process.env.CONTEXT || 'local'}</p>
        </div>
      </div>
    </div>
  );
};

export default ApiConnectivityTest;
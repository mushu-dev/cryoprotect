import React from 'react';
import ApiConnectivityTest from '../components/ApiConnectivityTest';

const ApiTestPage: React.FC = () => {
  return (
    <div className="container mx-auto py-8 px-4">
      <h1 className="text-2xl font-bold mb-6 text-center">API Connection Test</h1>
      
      <div className="max-w-lg mx-auto">
        <ApiConnectivityTest />
      </div>

      <div className="mt-8 max-w-lg mx-auto bg-white rounded-lg border border-gray-200 p-4 shadow-md">
        <h2 className="text-xl font-semibold mb-4">Troubleshooting</h2>
        
        <div className="space-y-4">
          <div>
            <h3 className="font-medium">If connection fails:</h3>
            <ul className="list-disc ml-5 text-sm mt-2">
              <li>Check that both the backend and frontend are deployed and running</li>
              <li>Verify that the <code className="bg-gray-100 px-1">NEXT_PUBLIC_API_URL</code> environment variable is set correctly</li>
              <li>Check CORS settings in the backend to allow requests from this frontend</li>
              <li>Check browser console for more detailed error messages</li>
            </ul>
          </div>
          
          <div>
            <h3 className="font-medium">Environment variables:</h3>
            <div className="bg-gray-100 p-2 text-xs rounded mt-2">
              <p><strong>NEXT_PUBLIC_API_URL:</strong> {process.env.NEXT_PUBLIC_API_URL || 'Not set'}</p>
              <p><strong>NEXTAUTH_URL:</strong> {process.env.NEXTAUTH_URL || 'Not set'}</p>
              <p><strong>NEXT_PUBLIC_ENVIRONMENT:</strong> {process.env.NEXT_PUBLIC_ENVIRONMENT || 'development'}</p>
            </div>
          </div>
        </div>
      </div>
    </div>
  );
};

export default ApiTestPage;
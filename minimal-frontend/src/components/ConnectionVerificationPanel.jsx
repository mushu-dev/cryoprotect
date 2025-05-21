/**
 * Connection Verification Panel
 * 
 * A component for checking and displaying the status of all connections
 * between the frontend and backend systems.
 */

import React, { useState, useEffect } from 'react';
import { verifyAllConnections } from '../utils/connection-verifier';

const statusColor = {
  success: 'bg-green-100 text-green-800 border-green-200',
  error: 'bg-red-100 text-red-800 border-red-200',
  pending: 'bg-yellow-100 text-yellow-800 border-yellow-200',
  unknown: 'bg-gray-100 text-gray-800 border-gray-200'
};

/**
 * Status badge component
 */
const StatusBadge = ({ status }) => {
  const getStatusProps = () => {
    if (status === 'success') {
      return {
        color: statusColor.success,
        text: 'Connected',
        icon: '✓'
      };
    } else if (status === 'error') {
      return {
        color: statusColor.error,
        text: 'Failed',
        icon: '✗'
      };
    } else if (status === 'pending') {
      return {
        color: statusColor.pending,
        text: 'Checking...',
        icon: '⟳'
      };
    } else {
      return {
        color: statusColor.unknown,
        text: 'Unknown',
        icon: '?'
      };
    }
  };

  const { color, text, icon } = getStatusProps();

  return (
    <div className={`px-3 py-1 inline-flex items-center rounded-full text-sm font-medium ${color}`}>
      <span className="mr-1">{icon}</span> {text}
    </div>
  );
};

/**
 * Connection Verification Panel Component
 */
export default function ConnectionVerificationPanel({ 
  convexClient = null, 
  authToken = null, 
  autoCheck = true 
}) {
  // State for verification results
  const [results, setResults] = useState({
    api: { status: 'unknown' },
    convex: { status: 'unknown' },
    auth: { status: 'unknown' },
    serviceRole: { status: 'unknown' },
    database: { status: 'unknown' }
  });
  
  // State for showing details
  const [showDetails, setShowDetails] = useState({
    api: false,
    convex: false,
    auth: false,
    serviceRole: false,
    database: false
  });
  
  // State for verification in progress
  const [isVerifying, setIsVerifying] = useState(false);
  
  // Overall status calculation
  const overallStatus = Object.values(results).every(r => r.status === 'success') 
    ? 'success' 
    : Object.values(results).some(r => r.status === 'error') 
      ? 'error' 
      : Object.values(results).some(r => r.status === 'pending') 
        ? 'pending' 
        : 'unknown';
  
  // Run verification
  const runVerification = async () => {
    // Reset results to pending
    setResults({
      api: { status: 'pending' },
      convex: { status: 'pending' },
      auth: { status: 'pending' },
      serviceRole: { status: 'pending' },
      database: { status: 'pending' }
    });
    
    setIsVerifying(true);
    
    try {
      const verificationResults = await verifyAllConnections(convexClient, authToken);
      
      // Update state with verification results
      setResults({
        api: { 
          status: verificationResults.api.success ? 'success' : 'error',
          message: verificationResults.api.message,
          details: verificationResults.api.details
        },
        convex: { 
          status: verificationResults.convex.success ? 'success' : 'error',
          message: verificationResults.convex.message,
          details: verificationResults.convex.details
        },
        auth: { 
          status: verificationResults.auth.success ? 'success' : 'error',
          message: verificationResults.auth.message,
          details: verificationResults.auth.details
        },
        serviceRole: { 
          status: verificationResults.serviceRole.success ? 'success' : 'error',
          message: verificationResults.serviceRole.message,
          details: verificationResults.serviceRole.details
        },
        database: { 
          status: verificationResults.database.success ? 'success' : 'error',
          message: verificationResults.database.message,
          details: verificationResults.database.details
        }
      });
    } catch (error) {
      console.error('Verification error:', error);
      
      // Update state with error
      setResults({
        api: { status: 'error', message: 'Verification failed', details: { error: error.message } },
        convex: { status: 'error', message: 'Verification failed', details: { error: error.message } },
        auth: { status: 'error', message: 'Verification failed', details: { error: error.message } },
        serviceRole: { status: 'error', message: 'Verification failed', details: { error: error.message } },
        database: { status: 'error', message: 'Verification failed', details: { error: error.message } }
      });
    } finally {
      setIsVerifying(false);
    }
  };
  
  // Toggle details visibility
  const toggleDetails = (service) => {
    setShowDetails({
      ...showDetails,
      [service]: !showDetails[service]
    });
  };
  
  // Run verification on mount if autoCheck is true
  useEffect(() => {
    if (autoCheck) {
      runVerification();
    }
  }, [autoCheck, convexClient, authToken]);
  
  return (
    <div className="border rounded-lg shadow-sm overflow-hidden">
      <div className="p-4 bg-gray-50 border-b flex justify-between items-center">
        <div>
          <h2 className="text-lg font-medium text-gray-900">Connection Status</h2>
          <p className="text-sm text-gray-500">Verify connectivity between frontend and backend systems</p>
        </div>
        <StatusBadge status={overallStatus} />
      </div>
      
      <div className="divide-y">
        {/* API Connection */}
        <div className="p-4">
          <div className="flex justify-between items-center mb-2">
            <div>
              <h3 className="font-medium">API Connection</h3>
              {results.api.message && (
                <p className="text-sm text-gray-500">{results.api.message}</p>
              )}
            </div>
            <div className="flex items-center space-x-3">
              <StatusBadge status={results.api.status} />
              <button 
                className="text-blue-600 hover:text-blue-800 text-sm"
                onClick={() => toggleDetails('api')}
              >
                {showDetails.api ? 'Hide' : 'Show'} Details
              </button>
            </div>
          </div>
          
          {showDetails.api && results.api.details && (
            <div className="mt-2 p-3 bg-gray-50 rounded-md text-sm font-mono overflow-x-auto">
              <pre>{JSON.stringify(results.api.details, null, 2)}</pre>
            </div>
          )}
        </div>
        
        {/* Convex Connection */}
        <div className="p-4">
          <div className="flex justify-between items-center mb-2">
            <div>
              <h3 className="font-medium">Convex Connection</h3>
              {results.convex.message && (
                <p className="text-sm text-gray-500">{results.convex.message}</p>
              )}
            </div>
            <div className="flex items-center space-x-3">
              <StatusBadge status={results.convex.status} />
              <button 
                className="text-blue-600 hover:text-blue-800 text-sm"
                onClick={() => toggleDetails('convex')}
              >
                {showDetails.convex ? 'Hide' : 'Show'} Details
              </button>
            </div>
          </div>
          
          {showDetails.convex && results.convex.details && (
            <div className="mt-2 p-3 bg-gray-50 rounded-md text-sm font-mono overflow-x-auto">
              <pre>{JSON.stringify(results.convex.details, null, 2)}</pre>
            </div>
          )}
        </div>
        
        {/* Authentication Connection */}
        <div className="p-4">
          <div className="flex justify-between items-center mb-2">
            <div>
              <h3 className="font-medium">Authentication Connection</h3>
              {results.auth.message && (
                <p className="text-sm text-gray-500">{results.auth.message}</p>
              )}
            </div>
            <div className="flex items-center space-x-3">
              <StatusBadge status={results.auth.status} />
              <button 
                className="text-blue-600 hover:text-blue-800 text-sm"
                onClick={() => toggleDetails('auth')}
              >
                {showDetails.auth ? 'Hide' : 'Show'} Details
              </button>
            </div>
          </div>
          
          {showDetails.auth && results.auth.details && (
            <div className="mt-2 p-3 bg-gray-50 rounded-md text-sm font-mono overflow-x-auto">
              <pre>{JSON.stringify(results.auth.details, null, 2)}</pre>
            </div>
          )}
        </div>
        
        {/* Service Role Connection */}
        <div className="p-4">
          <div className="flex justify-between items-center mb-2">
            <div>
              <h3 className="font-medium">Service Role Connection</h3>
              {results.serviceRole.message && (
                <p className="text-sm text-gray-500">{results.serviceRole.message}</p>
              )}
            </div>
            <div className="flex items-center space-x-3">
              <StatusBadge status={results.serviceRole.status} />
              <button 
                className="text-blue-600 hover:text-blue-800 text-sm"
                onClick={() => toggleDetails('serviceRole')}
              >
                {showDetails.serviceRole ? 'Hide' : 'Show'} Details
              </button>
            </div>
          </div>
          
          {showDetails.serviceRole && results.serviceRole.details && (
            <div className="mt-2 p-3 bg-gray-50 rounded-md text-sm font-mono overflow-x-auto">
              <pre>{JSON.stringify(results.serviceRole.details, null, 2)}</pre>
            </div>
          )}
        </div>
        
        {/* Database Connection */}
        <div className="p-4">
          <div className="flex justify-between items-center mb-2">
            <div>
              <h3 className="font-medium">Database Connection</h3>
              {results.database.message && (
                <p className="text-sm text-gray-500">{results.database.message}</p>
              )}
            </div>
            <div className="flex items-center space-x-3">
              <StatusBadge status={results.database.status} />
              <button 
                className="text-blue-600 hover:text-blue-800 text-sm"
                onClick={() => toggleDetails('database')}
              >
                {showDetails.database ? 'Hide' : 'Show'} Details
              </button>
            </div>
          </div>
          
          {showDetails.database && results.database.details && (
            <div className="mt-2 p-3 bg-gray-50 rounded-md text-sm font-mono overflow-x-auto">
              <pre>{JSON.stringify(results.database.details, null, 2)}</pre>
            </div>
          )}
        </div>
      </div>
      
      <div className="p-4 bg-gray-50 border-t flex justify-end">
        <button
          onClick={runVerification}
          disabled={isVerifying}
          className={`px-4 py-2 rounded-md text-white ${
            isVerifying ? 'bg-gray-400 cursor-not-allowed' : 'bg-blue-600 hover:bg-blue-700'
          }`}
        >
          {isVerifying ? 'Verifying...' : 'Verify Connections'}
        </button>
      </div>
    </div>
  );
}
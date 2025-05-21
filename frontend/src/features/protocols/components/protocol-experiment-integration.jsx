import React from 'react';

export function ProtocolExperimentIntegration({ protocolId, convexClient }) {
  return (
    <div className="bg-white dark:bg-gray-800 rounded-lg shadow-sm p-6">
      <div className="mb-6">
        <h2 className="text-xl font-semibold mb-2">Protocol-Experiment Integration</h2>
        <p className="text-gray-500 dark:text-gray-400">
          Create and manage experiments based on this protocol template.
        </p>
      </div>
      
      <div className="border rounded-md p-4 mb-6">
        <h3 className="text-lg font-medium mb-2">Protocol ID: {protocolId}</h3>
        <p className="text-sm text-gray-500 dark:text-gray-400 mb-4">
          This integration allows you to create standardized experiments based on this protocol.
        </p>
        
        <div className="flex flex-col space-y-2">
          <button 
            className="inline-flex justify-center items-center px-4 py-2 border border-transparent text-sm font-medium rounded-md text-white bg-blue-600 hover:bg-blue-700"
          >
            Create New Experiment
          </button>
          
          <button
            className="inline-flex justify-center items-center px-4 py-2 border border-gray-300 dark:border-gray-600 text-sm font-medium rounded-md text-gray-700 dark:text-gray-300 bg-white dark:bg-gray-700 hover:bg-gray-50 dark:hover:bg-gray-600"
          >
            View Associated Experiments
          </button>
        </div>
      </div>
      
      <div className="border rounded-md p-4">
        <h3 className="text-lg font-medium mb-2">Recent Experiments</h3>
        <p className="text-gray-500 dark:text-gray-400 italic">
          No experiments have been created from this protocol yet.
        </p>
      </div>
    </div>
  );
}
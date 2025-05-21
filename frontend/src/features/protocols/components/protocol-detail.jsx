import React, { useState } from 'react';
import { useRouter } from 'next/router';

/**
 * Protocol Detail Component
 * 
 * Component for displaying detailed protocol information
 */
export default function ProtocolDetail({ protocol }) {
  const router = useRouter();
  const [activeTab, setActiveTab] = useState('steps');
  const [showVersionHistory, setShowVersionHistory] = useState(false);
  
  if (!protocol) {
    return (
      <div className="flex justify-center items-center h-64">
        <div className="text-center">
          <div className="animate-spin rounded-full h-12 w-12 border-t-2 border-b-2 border-primary mx-auto mb-4"></div>
          <p className="text-muted-foreground">Loading protocol details...</p>
        </div>
      </div>
    );
  }
  
  return (
    <div className="bg-white dark:bg-gray-800 rounded-md shadow p-6">
      {/* Protocol Header */}
      <div className="mb-6">
        <h1 className="text-2xl font-bold mb-2">{protocol.name}</h1>
        <p className="text-gray-500 dark:text-gray-400 mb-4">{protocol.description}</p>
        
        <div className="flex flex-wrap gap-2 text-sm">
          <span className="px-2 py-1 bg-blue-100 dark:bg-blue-900/30 text-blue-800 dark:text-blue-300 rounded-md">
            {protocol.category || 'General'}
          </span>
          <span className="px-2 py-1 bg-gray-100 dark:bg-gray-700 text-gray-800 dark:text-gray-300 rounded-md">
            Version {protocol.version || '1.0'}
          </span>
          {protocol.author && (
            <span className="px-2 py-1 bg-gray-100 dark:bg-gray-700 text-gray-800 dark:text-gray-300 rounded-md">
              Created by {protocol.author}
            </span>
          )}
        </div>
      </div>
      
      {/* Tabs Navigation */}
      <div className="border-b border-gray-200 dark:border-gray-700 mb-6">
        <nav className="-mb-px flex space-x-8">
          <button
            onClick={() => setActiveTab('steps')}
            className={`pb-4 px-1 border-b-2 font-medium text-sm ${
              activeTab === 'steps'
                ? 'border-primary text-primary'
                : 'border-transparent text-gray-500 hover:text-gray-700 dark:text-gray-400 dark:hover:text-gray-300'
            }`}
          >
            Protocol Steps
          </button>
          <button
            onClick={() => setActiveTab('materials')}
            className={`pb-4 px-1 border-b-2 font-medium text-sm ${
              activeTab === 'materials'
                ? 'border-primary text-primary'
                : 'border-transparent text-gray-500 hover:text-gray-700 dark:text-gray-400 dark:hover:text-gray-300'
            }`}
          >
            Materials
          </button>
          <button
            onClick={() => setActiveTab('experiments')}
            className={`pb-4 px-1 border-b-2 font-medium text-sm ${
              activeTab === 'experiments'
                ? 'border-primary text-primary'
                : 'border-transparent text-gray-500 hover:text-gray-700 dark:text-gray-400 dark:hover:text-gray-300'
            }`}
          >
            Related Experiments
          </button>
          <button
            onClick={() => setActiveTab('versions')}
            className={`pb-4 px-1 border-b-2 font-medium text-sm ${
              activeTab === 'versions'
                ? 'border-primary text-primary'
                : 'border-transparent text-gray-500 hover:text-gray-700 dark:text-gray-400 dark:hover:text-gray-300'
            }`}
          >
            Version History
          </button>
        </nav>
      </div>
      
      {/* Tab Content */}
      <div>
        {/* Steps Tab */}
        {activeTab === 'steps' && (
          <div>
            <h2 className="text-lg font-semibold mb-4">Protocol Steps</h2>
            
            {!protocol.steps || protocol.steps.length === 0 ? (
              <p className="text-gray-500 dark:text-gray-400 italic">No steps defined for this protocol</p>
            ) : (
              <ol className="space-y-4">
                {protocol.steps.map((step, index) => (
                  <li key={index} className="p-4 bg-gray-50 dark:bg-gray-700 rounded-md">
                    <div className="flex items-center mb-2">
                      <span className="bg-blue-600 text-white w-6 h-6 rounded-full flex items-center justify-center mr-2">
                        {index + 1}
                      </span>
                      <h3 className="font-medium">{step.title || 'Step ' + (index + 1)}</h3>
                    </div>
                    
                    <p className="mb-2 text-sm">{step.description}</p>
                    
                    <div className="flex flex-wrap gap-2 text-xs text-gray-500 dark:text-gray-400">
                      {step.duration && <span>Duration: {step.duration}</span>}
                      {step.temperature && <span>Temperature: {step.temperature}</span>}
                    </div>
                    
                    {step.notes && (
                      <div className="mt-2 p-2 bg-yellow-50 dark:bg-yellow-900/20 text-sm rounded">
                        <strong>Notes:</strong> {step.notes}
                      </div>
                    )}
                  </li>
                ))}
              </ol>
            )}
          </div>
        )}
        
        {/* Materials Tab */}
        {activeTab === 'materials' && (
          <div>
            <h2 className="text-lg font-semibold mb-4">Required Materials</h2>
            
            {!protocol.materials || protocol.materials.length === 0 ? (
              <p className="text-gray-500 dark:text-gray-400 italic">No materials listed for this protocol</p>
            ) : (
              <ul className="space-y-2">
                {protocol.materials.map((material, index) => (
                  <li key={index} className="flex items-center p-2 bg-gray-50 dark:bg-gray-700 rounded-md">
                    <span className="h-2 w-2 bg-blue-600 rounded-full mr-2"></span>
                    <span>{typeof material === 'string' ? material : material.name || 'Unnamed material'}</span>
                  </li>
                ))}
              </ul>
            )}
          </div>
        )}
        
        {/* Experiments Tab */}
        {activeTab === 'experiments' && (
          <div>
            <div className="flex justify-between items-center mb-4">
              <h2 className="text-lg font-semibold">Related Experiments</h2>
              <button 
                onClick={() => router.push(`/protocols/${protocol.id}/experiment-integration`)}
                className="px-3 py-1.5 text-sm bg-blue-600 text-white rounded-md hover:bg-blue-700"
              >
                Create Experiment
              </button>
            </div>
            
            <p className="text-gray-500 dark:text-gray-400 italic">No experiments have been created from this protocol yet</p>
          </div>
        )}
        
        {/* Versions Tab */}
        {activeTab === 'versions' && (
          <div>
            <h2 className="text-lg font-semibold mb-4">Version History</h2>
            
            <div className="space-y-3">
              <div className="p-3 bg-blue-50 dark:bg-blue-900/20 rounded-md border-l-4 border-blue-600">
                <div className="flex justify-between">
                  <div>
                    <p className="font-medium">Version {protocol.version || '1.0'} (Current)</p>
                    <p className="text-sm text-gray-500 dark:text-gray-400">
                      {protocol.updated_at ? new Date(protocol.updated_at).toLocaleDateString() : 'No date available'}
                    </p>
                  </div>
                  <div>
                    <button
                      className="text-sm text-blue-600 dark:text-blue-400 underline"
                      onClick={() => {/* View this version */}}
                    >
                      Current Version
                    </button>
                  </div>
                </div>
              </div>
              
              {/* Mock older versions */}
              <div className="p-3 bg-gray-50 dark:bg-gray-700 rounded-md">
                <div className="flex justify-between">
                  <div>
                    <p className="font-medium">Version 1.0</p>
                    <p className="text-sm text-gray-500 dark:text-gray-400">
                      {new Date(Date.now() - 30 * 24 * 60 * 60 * 1000).toLocaleDateString()}
                    </p>
                  </div>
                  <div>
                    <button
                      className="text-sm text-blue-600 dark:text-blue-400 underline"
                      onClick={() => {/* View this version */}}
                    >
                      View
                    </button>
                  </div>
                </div>
              </div>
            </div>
          </div>
        )}
      </div>
    </div>
  );
}
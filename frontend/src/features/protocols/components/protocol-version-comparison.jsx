import React from 'react';

export default function ProtocolVersionComparison({ 
  protocolA, // Older version
  protocolB, // Newer version (current)
  onViewVersion
}) {
  // If protocols aren't passed, show loading state
  if (!protocolA || !protocolB) {
    return (
      <div className="flex justify-center items-center h-64">
        <div className="text-center">
          <div className="animate-spin rounded-full h-12 w-12 border-t-2 border-b-2 border-primary mx-auto mb-4"></div>
          <p className="text-muted-foreground">Loading version comparison...</p>
        </div>
      </div>
    );
  }
  
  // Find differences between protocols
  const differences = {
    basicInfo: {
      name: protocolA.name !== protocolB.name,
      description: protocolA.description !== protocolB.description,
      objective: protocolA.objective !== protocolB.objective
    },
    materials: !arraysEqual(protocolA.materials || [], protocolB.materials || []),
    steps: compareSteps(protocolA.steps || [], protocolB.steps || [])
  };
  
  // Check if there are any differences at all
  const hasDifferences = 
    differences.basicInfo.name || 
    differences.basicInfo.description || 
    differences.basicInfo.objective ||
    differences.materials ||
    differences.steps.added.length > 0 ||
    differences.steps.removed.length > 0 ||
    differences.steps.modified.length > 0;
  
  return (
    <div className="bg-white dark:bg-gray-800 rounded-md shadow p-6">
      <div className="mb-6">
        <h2 className="text-xl font-semibold mb-1">Version Comparison</h2>
        <p className="text-gray-500 dark:text-gray-400">
          Comparing version {protocolA.version || '1.0'} with version {protocolB.version || '1.1'}
        </p>
      </div>
      
      {!hasDifferences ? (
        <div className="p-4 bg-yellow-50 dark:bg-yellow-900/20 rounded-md mb-6">
          <p className="font-medium">No differences detected</p>
          <p className="text-sm text-gray-600 dark:text-gray-300">
            These protocol versions are identical.
          </p>
        </div>
      ) : (
        <div className="space-y-6">
          {/* Basic information changes */}
          <div>
            <h3 className="text-lg font-medium mb-3">Basic Information</h3>
            <div className="space-y-2">
              {differences.basicInfo.name && (
                <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                  <div className="p-3 bg-red-50 dark:bg-red-900/20 rounded-md">
                    <p className="font-medium text-sm text-red-800 dark:text-red-300">Previous</p>
                    <p>{protocolA.name}</p>
                  </div>
                  <div className="p-3 bg-green-50 dark:bg-green-900/20 rounded-md">
                    <p className="font-medium text-sm text-green-800 dark:text-green-300">Current</p>
                    <p>{protocolB.name}</p>
                  </div>
                </div>
              )}
              
              {differences.basicInfo.description && (
                <div>
                  <h4 className="text-sm font-medium mb-2">Description</h4>
                  <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                    <div className="p-3 bg-red-50 dark:bg-red-900/20 rounded-md">
                      <p className="font-medium text-sm text-red-800 dark:text-red-300">Previous</p>
                      <p className="text-sm">{protocolA.description}</p>
                    </div>
                    <div className="p-3 bg-green-50 dark:bg-green-900/20 rounded-md">
                      <p className="font-medium text-sm text-green-800 dark:text-green-300">Current</p>
                      <p className="text-sm">{protocolB.description}</p>
                    </div>
                  </div>
                </div>
              )}
              
              {differences.basicInfo.objective && (
                <div>
                  <h4 className="text-sm font-medium mb-2">Objective</h4>
                  <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                    <div className="p-3 bg-red-50 dark:bg-red-900/20 rounded-md">
                      <p className="font-medium text-sm text-red-800 dark:text-red-300">Previous</p>
                      <p className="text-sm">{protocolA.objective}</p>
                    </div>
                    <div className="p-3 bg-green-50 dark:bg-green-900/20 rounded-md">
                      <p className="font-medium text-sm text-green-800 dark:text-green-300">Current</p>
                      <p className="text-sm">{protocolB.objective}</p>
                    </div>
                  </div>
                </div>
              )}
            </div>
          </div>
          
          {/* Materials changes */}
          {differences.materials && (
            <div>
              <h3 className="text-lg font-medium mb-3">Materials</h3>
              <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                <div className="p-3 bg-red-50 dark:bg-red-900/20 rounded-md">
                  <p className="font-medium text-sm text-red-800 dark:text-red-300">Previous Materials</p>
                  {protocolA.materials && protocolA.materials.length > 0 ? (
                    <ul className="mt-2 list-disc pl-5">
                      {protocolA.materials.map((material, index) => (
                        <li key={index} className="text-sm">
                          {typeof material === 'string' ? material : material.name || 'Unnamed material'}
                        </li>
                      ))}
                    </ul>
                  ) : (
                    <p className="text-sm italic">No materials listed</p>
                  )}
                </div>
                <div className="p-3 bg-green-50 dark:bg-green-900/20 rounded-md">
                  <p className="font-medium text-sm text-green-800 dark:text-green-300">Current Materials</p>
                  {protocolB.materials && protocolB.materials.length > 0 ? (
                    <ul className="mt-2 list-disc pl-5">
                      {protocolB.materials.map((material, index) => (
                        <li key={index} className="text-sm">
                          {typeof material === 'string' ? material : material.name || 'Unnamed material'}
                        </li>
                      ))}
                    </ul>
                  ) : (
                    <p className="text-sm italic">No materials listed</p>
                  )}
                </div>
              </div>
            </div>
          )}
          
          {/* Steps changes */}
          <div>
            <h3 className="text-lg font-medium mb-3">Protocol Steps</h3>
            
            {/* Added steps */}
            {differences.steps.added.length > 0 && (
              <div className="mb-4">
                <h4 className="text-sm font-medium text-green-800 dark:text-green-300 mb-2">
                  Added Steps ({differences.steps.added.length})
                </h4>
                <div className="space-y-2">
                  {differences.steps.added.map((step, index) => (
                    <div key={index} className="p-3 bg-green-50 dark:bg-green-900/20 rounded-md">
                      <h5 className="font-medium">{step.title || 'Unnamed Step'}</h5>
                      <p className="text-sm">{step.description}</p>
                      <div className="flex flex-wrap gap-2 mt-1 text-xs text-gray-500 dark:text-gray-400">
                        {step.duration && <span>Duration: {step.duration}</span>}
                        {step.temperature && <span>Temperature: {step.temperature}</span>}
                      </div>
                    </div>
                  ))}
                </div>
              </div>
            )}
            
            {/* Removed steps */}
            {differences.steps.removed.length > 0 && (
              <div className="mb-4">
                <h4 className="text-sm font-medium text-red-800 dark:text-red-300 mb-2">
                  Removed Steps ({differences.steps.removed.length})
                </h4>
                <div className="space-y-2">
                  {differences.steps.removed.map((step, index) => (
                    <div key={index} className="p-3 bg-red-50 dark:bg-red-900/20 rounded-md">
                      <h5 className="font-medium">{step.title || 'Unnamed Step'}</h5>
                      <p className="text-sm">{step.description}</p>
                      <div className="flex flex-wrap gap-2 mt-1 text-xs text-gray-500 dark:text-gray-400">
                        {step.duration && <span>Duration: {step.duration}</span>}
                        {step.temperature && <span>Temperature: {step.temperature}</span>}
                      </div>
                    </div>
                  ))}
                </div>
              </div>
            )}
            
            {/* Modified steps */}
            {differences.steps.modified.length > 0 && (
              <div>
                <h4 className="text-sm font-medium text-blue-800 dark:text-blue-300 mb-2">
                  Modified Steps ({differences.steps.modified.length})
                </h4>
                <div className="space-y-4">
                  {differences.steps.modified.map((mod, index) => (
                    <div key={index} className="border border-blue-300 dark:border-blue-700 rounded-md">
                      <div className="p-2 bg-blue-50 dark:bg-blue-900/20 border-b border-blue-300 dark:border-blue-700">
                        <h5 className="font-medium">Step {mod.index + 1}: {mod.current.title || 'Unnamed Step'}</h5>
                      </div>
                      <div className="grid grid-cols-1 md:grid-cols-2 gap-4 p-3">
                        <div className="p-3 bg-red-50 dark:bg-red-900/20 rounded-md">
                          <p className="font-medium text-sm text-red-800 dark:text-red-300">Previous</p>
                          <p className="text-sm">{mod.previous.description}</p>
                          <div className="flex flex-wrap gap-2 mt-1 text-xs text-gray-500 dark:text-gray-400">
                            {mod.previous.duration && <span>Duration: {mod.previous.duration}</span>}
                            {mod.previous.temperature && <span>Temperature: {mod.previous.temperature}</span>}
                          </div>
                        </div>
                        <div className="p-3 bg-green-50 dark:bg-green-900/20 rounded-md">
                          <p className="font-medium text-sm text-green-800 dark:text-green-300">Current</p>
                          <p className="text-sm">{mod.current.description}</p>
                          <div className="flex flex-wrap gap-2 mt-1 text-xs text-gray-500 dark:text-gray-400">
                            {mod.current.duration && <span>Duration: {mod.current.duration}</span>}
                            {mod.current.temperature && <span>Temperature: {mod.current.temperature}</span>}
                          </div>
                        </div>
                      </div>
                    </div>
                  ))}
                </div>
              </div>
            )}
            
            {/* No changes to steps */}
            {differences.steps.added.length === 0 && 
             differences.steps.removed.length === 0 && 
             differences.steps.modified.length === 0 && (
              <p className="text-sm italic">No changes to protocol steps</p>
            )}
          </div>
        </div>
      )}
      
      <div className="flex justify-between items-center mt-8 pt-4 border-t">
        <button
          onClick={() => onViewVersion && onViewVersion(protocolA)}
          className="px-4 py-2 border border-gray-300 dark:border-gray-600 rounded-md text-sm font-medium"
        >
          View Previous Version
        </button>
        <button
          onClick={() => onViewVersion && onViewVersion(protocolB)}
          className="px-4 py-2 bg-blue-600 text-white rounded-md text-sm font-medium hover:bg-blue-700"
        >
          View Current Version
        </button>
      </div>
    </div>
  );
}

// Helper function to compare arrays
function arraysEqual(a, b) {
  if (a === b) return true;
  if (a == null || b == null) return false;
  if (a.length !== b.length) return false;
  
  // If they are array of objects, convert to strings for comparison
  const aStrings = a.map(item => typeof item === 'string' ? item : JSON.stringify(item));
  const bStrings = b.map(item => typeof item === 'string' ? item : JSON.stringify(item));
  
  // Sort and join for comparison
  return aStrings.sort().join(',') === bStrings.sort().join(',');
}

// Helper function to compare steps
function compareSteps(stepsA, stepsB) {
  const result = {
    added: [],
    removed: [],
    modified: []
  };
  
  // Simple algorithm for matching steps (can be improved for better matching)
  const usedBIndices = new Set();
  
  // Find modified and removed steps
  stepsA.forEach((stepA, indexA) => {
    // Try to find matching step in B
    const indexB = stepsB.findIndex((stepB, indexB) => 
      !usedBIndices.has(indexB) && stepB.title && stepA.title && stepB.title === stepA.title
    );
    
    if (indexB !== -1) {
      // Step exists in both versions, check if modified
      const stepB = stepsB[indexB];
      usedBIndices.add(indexB);
      
      if (
        stepA.description !== stepB.description ||
        stepA.duration !== stepB.duration ||
        stepA.temperature !== stepB.temperature ||
        stepA.notes !== stepB.notes
      ) {
        result.modified.push({
          index: indexA,
          previous: stepA,
          current: stepB
        });
      }
    } else {
      // Step was removed
      result.removed.push(stepA);
    }
  });
  
  // Find added steps
  stepsB.forEach((stepB, indexB) => {
    if (!usedBIndices.has(indexB)) {
      result.added.push(stepB);
    }
  });
  
  return result;
}
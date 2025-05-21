import React from 'react';
import Link from 'next/link';

/**
 * Card component for displaying experiment summary with selection capability
 */
export default function ExperimentCard({ 
  experiment, 
  isSelected = false, 
  onSelect, 
  showSelectOption = false 
}) {
  if (!experiment) return null;

  // Extract key metrics
  const viability = experiment.results?.viability || 'N/A';
  const recovery = experiment.results?.recovery || 'N/A';

  // Get status styling
  const getStatusStyles = () => {
    if (experiment.status === 'Completed') {
      return 'bg-green-100 text-green-800';
    } else if (experiment.status === 'In Progress') {
      return 'bg-blue-100 text-blue-800';
    } else {
      return 'bg-gray-100 text-gray-800';
    }
  };

  // Handle card selection for comparison
  const handleSelect = (e) => {
    e.preventDefault();
    e.stopPropagation();
    if (onSelect) {
      onSelect(experiment.id);
    }
  };

  return (
    <div 
      className="experiment-card rounded-lg border bg-card p-6 shadow-sm transition-all duration-200"
      data-testid="experiment-card"
    >
      <div className="flex items-center justify-between mb-4">
        <h3 className="text-lg font-semibold">{experiment.title}</h3>
        <span className={`inline-flex items-center rounded-full px-2.5 py-0.5 text-xs font-medium ${getStatusStyles()}`}>
          {experiment.status}
        </span>
      </div>
      
      <p className="text-muted-foreground mb-4 line-clamp-2">
        {experiment.description}
      </p>
      
      {/* Key metrics */}
      <div className="grid grid-cols-2 gap-2 mb-4" data-testid="experiment-chart">
        <div className="border rounded-md p-2">
          <div className="text-xs text-gray-500">Viability</div>
          <div className="font-medium">{viability}</div>
          <div className="w-full bg-gray-200 h-2 mt-1 rounded-full overflow-hidden">
            <div 
              className="bg-green-500 h-full" 
              style={{ width: viability !== 'N/A' ? viability : '0%' }}>
            </div>
          </div>
        </div>
        <div className="border rounded-md p-2">
          <div className="text-xs text-gray-500">Recovery</div>
          <div className="font-medium">{recovery}</div>
          <div className="w-full bg-gray-200 h-2 mt-1 rounded-full overflow-hidden">
            <div 
              className="bg-blue-500 h-full" 
              style={{ width: recovery !== 'N/A' ? recovery : '0%' }}>
            </div>
          </div>
        </div>
      </div>
      
      {/* Quick info */}
      <div className="flex flex-col gap-1 mb-4 text-sm">
        <div className="flex justify-between">
          <span className="text-gray-500">Date:</span>
          <time>{experiment.date}</time>
        </div>
        <div className="flex justify-between">
          <span className="text-gray-500">Cell Type:</span>
          <span className="truncate ml-2 max-w-[150px]">{experiment.cellType}</span>
        </div>
        <div className="flex justify-between">
          <span className="text-gray-500">Protocol:</span>
          <span>{experiment.protocol.name}</span>
        </div>
      </div>
      
      <div className="flex justify-between items-center mt-4">
        {showSelectOption ? (
          <button 
            onClick={handleSelect}
            className={`text-sm font-medium ${
              isSelected ? 'text-blue-700' : 'text-blue-600 hover:text-blue-800'
            }`}
            data-testid="experiment-select"
          >
            {isSelected ? '✓ Selected for comparison' : 'Select for comparison'}
          </button>
        ) : (
          <span className="text-sm text-muted-foreground">
            {experiment.cryoprotectants?.map(cp => cp.name).join(', ')}
          </span>
        )}
        
        <Link href={`/experiments/${experiment.id}`}>
          <a className="text-primary hover:underline text-sm">View Details</a>
        </Link>
      </div>
    </div>
  );
}
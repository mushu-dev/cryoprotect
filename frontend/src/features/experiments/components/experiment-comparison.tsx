import React from 'react';
import Link from 'next/link';

interface ExperimentComparisonProps {
  experiments: any[];
  onRemove: (id: string | string[]) => void;
}

/**
 * Component for comparing multiple experiments
 */
export default function ExperimentComparison({ experiments, onRemove }: ExperimentComparisonProps) {
  // If no experiments selected
  if (!experiments || experiments.length === 0) {
    return (
      <div className="bg-muted p-6 rounded-lg border text-center">
        <p className="text-muted-foreground mb-2">No experiments selected for comparison</p>
        <p className="text-sm text-muted-foreground">
          Select experiments from the list below to compare their properties and results
        </p>
      </div>
    );
  }
  
  // Helper function to calculate average for a metric
  const calculateAverage = (metric: string): string => {
    const values = experiments
      .map(exp => parseInt(exp.results[metric], 10))
      .filter(val => !isNaN(val));
    
    if (values.length === 0) return 'N/A';
    
    const avg = values.reduce((sum, val) => sum + val, 0) / values.length;
    return `${Math.round(avg)}%`;
  };
  
  return (
    <div className="bg-card rounded-lg border shadow-sm overflow-hidden">
      <div className="p-4 bg-muted border-b flex justify-between items-center">
        <h3 className="font-medium">Comparing {experiments.length} Experiments</h3>
        <div>
          <button 
            onClick={() => onRemove(experiments.map(exp => exp.id))}
            className="text-xs text-muted-foreground hover:text-foreground"
          >
            Clear All
          </button>
        </div>
      </div>
      
      <div className="p-4">
        <div className="overflow-x-auto">
          <table className="w-full divide-y divide-border">
            <thead>
              <tr>
                <th className="px-4 py-2 text-left text-xs font-medium text-muted-foreground uppercase tracking-wider">
                  Experiment
                </th>
                <th className="px-4 py-2 text-left text-xs font-medium text-muted-foreground uppercase tracking-wider">
                  Cell Type
                </th>
                <th className="px-4 py-2 text-left text-xs font-medium text-muted-foreground uppercase tracking-wider">
                  Cryoprotectants
                </th>
                <th className="px-4 py-2 text-left text-xs font-medium text-muted-foreground uppercase tracking-wider">
                  Viability
                </th>
                <th className="px-4 py-2 text-left text-xs font-medium text-muted-foreground uppercase tracking-wider">
                  Recovery
                </th>
                <th className="px-4 py-2 text-left text-xs font-medium text-muted-foreground uppercase tracking-wider">
                  Functionality
                </th>
                <th className="px-4 py-2 text-right text-xs font-medium text-muted-foreground uppercase tracking-wider">
                  Actions
                </th>
              </tr>
            </thead>
            <tbody className="divide-y divide-border">
              {experiments.map(exp => (
                <tr key={exp.id} className="hover:bg-muted/50">
                  <td className="px-4 py-2 whitespace-nowrap">
                    <Link href={`/experiments/${exp.id}`}>
                      <span className="text-primary hover:underline font-medium">{exp.title}</span>
                    </Link>
                  </td>
                  <td className="px-4 py-2 whitespace-nowrap text-sm">{exp.cellType}</td>
                  <td className="px-4 py-2 whitespace-nowrap text-sm">
                    {exp.cryoprotectants.map(cp => cp.name).join(', ')}
                  </td>
                  <td className="px-4 py-2 whitespace-nowrap text-sm font-medium">{exp.results.viability}</td>
                  <td className="px-4 py-2 whitespace-nowrap text-sm font-medium">{exp.results.recovery}</td>
                  <td className="px-4 py-2 whitespace-nowrap text-sm font-medium">{exp.results.functionality}</td>
                  <td className="px-4 py-2 whitespace-nowrap text-right text-sm">
                    <button 
                      onClick={() => onRemove(exp.id)}
                      className="text-red-500 hover:text-red-700"
                    >
                      Remove
                    </button>
                  </td>
                </tr>
              ))}
              
              {/* Summary row */}
              <tr className="bg-muted">
                <td colSpan={3} className="px-4 py-2 whitespace-nowrap text-sm font-medium">
                  Average/Summary
                </td>
                <td className="px-4 py-2 whitespace-nowrap text-sm font-medium">{calculateAverage('viability')}</td>
                <td className="px-4 py-2 whitespace-nowrap text-sm font-medium">{calculateAverage('recovery')}</td>
                <td className="px-4 py-2 whitespace-nowrap text-sm font-medium">{calculateAverage('functionality')}</td>
                <td className="px-4 py-2 whitespace-nowrap"></td>
              </tr>
            </tbody>
          </table>
        </div>
        
        {/* Visualization section */}
        <div className="mt-6">
          <h4 className="text-sm font-medium mb-3">Comparative Analysis</h4>
          <div className="bg-muted/30 rounded-lg p-4">
            {/* This would be replaced with a real chart in a production app */}
            <div className="space-y-4">
              {/* Viability comparison */}
              <div>
                <div className="flex justify-between items-center text-xs mb-1">
                  <span className="font-medium">Viability</span>
                </div>
                <div className="relative h-8 bg-muted rounded-lg">
                  {experiments.map((exp, index) => {
                    const viability = parseInt(exp.results.viability, 10) || 0;
                    const width = `${viability}%`;
                    
                    // Generate a color based on the index
                    const colors = [
                      'bg-primary', 'bg-blue-500', 'bg-green-500', 
                      'bg-yellow-500', 'bg-purple-500', 'bg-pink-500'
                    ];
                    const colorClass = colors[index % colors.length];
                    
                    return (
                      <div
                        key={exp.id}
                        className={`absolute top-0 h-8 ${colorClass} rounded-lg opacity-75`}
                        style={{ 
                          width, 
                          left: 0, 
                          zIndex: experiments.length - index // To ensure shorter bars are visible
                        }}
                      >
                        <span className="absolute inset-0 flex items-center justify-center text-white text-xs font-medium">
                          {exp.results.viability}
                        </span>
                      </div>
                    );
                  })}
                  
                  {/* Scale indicators */}
                  <div className="absolute inset-x-0 top-full mt-1 flex justify-between text-xs text-muted-foreground">
                    <span>0%</span>
                    <span>25%</span>
                    <span>50%</span>
                    <span>75%</span>
                    <span>100%</span>
                  </div>
                </div>
              </div>
              
              {/* Recovery comparison */}
              <div className="mt-8">
                <div className="flex justify-between items-center text-xs mb-1">
                  <span className="font-medium">Recovery</span>
                </div>
                <div className="relative h-8 bg-muted rounded-lg">
                  {experiments.map((exp, index) => {
                    const recovery = parseInt(exp.results.recovery, 10) || 0;
                    const width = `${recovery}%`;
                    
                    // Generate a color based on the index
                    const colors = [
                      'bg-primary', 'bg-blue-500', 'bg-green-500', 
                      'bg-yellow-500', 'bg-purple-500', 'bg-pink-500'
                    ];
                    const colorClass = colors[index % colors.length];
                    
                    return (
                      <div
                        key={exp.id}
                        className={`absolute top-0 h-8 ${colorClass} rounded-lg opacity-75`}
                        style={{ 
                          width, 
                          left: 0, 
                          zIndex: experiments.length - index 
                        }}
                      >
                        <span className="absolute inset-0 flex items-center justify-center text-white text-xs font-medium">
                          {exp.results.recovery}
                        </span>
                      </div>
                    );
                  })}
                  
                  {/* Scale indicators */}
                  <div className="absolute inset-x-0 top-full mt-1 flex justify-between text-xs text-muted-foreground">
                    <span>0%</span>
                    <span>25%</span>
                    <span>50%</span>
                    <span>75%</span>
                    <span>100%</span>
                  </div>
                </div>
              </div>
              
              {/* Functionality comparison */}
              <div className="mt-8">
                <div className="flex justify-between items-center text-xs mb-1">
                  <span className="font-medium">Functionality</span>
                </div>
                <div className="relative h-8 bg-muted rounded-lg">
                  {experiments.map((exp, index) => {
                    const functionality = parseInt(exp.results.functionality, 10) || 0;
                    const width = `${functionality}%`;
                    
                    // Generate a color based on the index
                    const colors = [
                      'bg-primary', 'bg-blue-500', 'bg-green-500', 
                      'bg-yellow-500', 'bg-purple-500', 'bg-pink-500'
                    ];
                    const colorClass = colors[index % colors.length];
                    
                    return (
                      <div
                        key={exp.id}
                        className={`absolute top-0 h-8 ${colorClass} rounded-lg opacity-75`}
                        style={{ 
                          width, 
                          left: 0, 
                          zIndex: experiments.length - index 
                        }}
                      >
                        <span className="absolute inset-0 flex items-center justify-center text-white text-xs font-medium">
                          {exp.results.functionality}
                        </span>
                      </div>
                    );
                  })}
                  
                  {/* Scale indicators */}
                  <div className="absolute inset-x-0 top-full mt-1 flex justify-between text-xs text-muted-foreground">
                    <span>0%</span>
                    <span>25%</span>
                    <span>50%</span>
                    <span>75%</span>
                    <span>100%</span>
                  </div>
                </div>
              </div>
            </div>
            
            {/* Legend */}
            <div className="mt-10 flex flex-wrap gap-4">
              {experiments.map((exp, index) => {
                const colors = [
                  'bg-primary', 'bg-blue-500', 'bg-green-500', 
                  'bg-yellow-500', 'bg-purple-500', 'bg-pink-500'
                ];
                const colorClass = colors[index % colors.length];
                
                return (
                  <div key={exp.id} className="flex items-center">
                    <div className={`w-3 h-3 rounded-full ${colorClass} mr-2`}></div>
                    <span className="text-xs">{exp.title}</span>
                  </div>
                );
              })}
            </div>
          </div>
        </div>
      </div>
    </div>
  );
}
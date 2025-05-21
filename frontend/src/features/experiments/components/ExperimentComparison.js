import React from 'react';
import { prepareExperimentComparisonData } from '../utils/data-transformation';

/**
 * Component for comparing multiple experiments
 */
export default function ExperimentComparison({ experiments, onRemove }) {
  if (!experiments || !Array.isArray(experiments) || experiments.length === 0) {
    return (
      <div className="text-center p-6 bg-gray-50 rounded-lg border" data-testid="empty-comparison">
        <p className="text-gray-500">No experiments selected for comparison</p>
        <p className="text-sm text-gray-400 mt-2">Select up to 4 experiments to compare</p>
      </div>
    );
  }

  const { comparisonData, labels } = prepareExperimentComparisonData(experiments);
  const colorMap = ['#3B82F6', '#10B981', '#F59E0B', '#8B5CF6']; // blue, green, amber, purple

  return (
    <div className="bg-white rounded-lg border shadow-sm overflow-hidden" data-testid="experiment-comparison">
      <div className="p-4 bg-gray-50 border-b font-medium">
        <div className="flex justify-between items-center">
          <h3>Experiment Comparison</h3>
          <button
            onClick={() => onRemove && onRemove([])} // Clear all selections
            className="text-xs text-gray-500 hover:text-red-500 hover:underline"
            data-testid="clear-comparison"
          >
            Clear All
          </button>
        </div>
      </div>

      <div className="p-4">
        {/* Legend */}
        <div className="mb-6 flex flex-wrap gap-4">
          {experiments.map((exp, index) => (
            <div key={exp.id} className="flex items-center" data-testid="comparison-item">
              <div 
                className="w-3 h-3 rounded-full mr-2" 
                style={{ backgroundColor: colorMap[index % colorMap.length] }}
              />
              <span className="text-sm font-medium">{exp.title}</span>
              <button
                onClick={() => onRemove && onRemove(exp.id)}
                className="ml-2 text-gray-400 hover:text-red-500"
                data-testid="remove-comparison-item"
              >
                <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
                  <line x1="18" y1="6" x2="6" y2="18"></line>
                  <line x1="6" y1="6" x2="18" y2="18"></line>
                </svg>
              </button>
            </div>
          ))}
        </div>

        {/* Chart */}
        <div className="h-64 flex items-end justify-around p-4 border-b border-l" data-testid="experiment-chart">
          {labels.map((label, labelIndex) => (
            <div key={label} className="flex items-end gap-2">
              {comparisonData.map((dataset, datasetIndex) => {
                // Calculate bar height as percentage of container height
                const value = dataset.data[labelIndex];
                const heightPercentage = value;
                
                return (
                  <div key={datasetIndex} className="w-10 flex flex-col items-center">
                    <div 
                      className="w-8 rounded-t-sm transition-all duration-500 ease-in-out"
                      style={{ 
                        height: `${heightPercentage}%`,
                        backgroundColor: colorMap[datasetIndex % colorMap.length]
                      }}
                    />
                    <span className="text-xs mt-2 rotate-45 origin-left">{value}%</span>
                  </div>
                );
              })}
              <div className="text-xs font-medium mt-8">{label}</div>
            </div>
          ))}
        </div>

        {/* Data Table */}
        <div className="mt-6 overflow-x-auto">
          <table className="w-full text-sm divide-y divide-gray-200" data-testid="comparison-table">
            <thead>
              <tr>
                <th className="px-4 py-2 text-left font-medium text-gray-500">Experiment</th>
                {labels.map(label => (
                  <th key={label} className="px-4 py-2 text-left font-medium text-gray-500">
                    {label}
                  </th>
                ))}
                <th className="px-4 py-2 text-left font-medium text-gray-500">
                  Protocol
                </th>
                <th className="px-4 py-2 text-left font-medium text-gray-500">
                  Cell Type
                </th>
              </tr>
            </thead>
            <tbody className="divide-y divide-gray-200">
              {experiments.map((exp, index) => (
                <tr key={exp.id} className="hover:bg-gray-50" data-testid="comparison-row">
                  <td className="px-4 py-2 font-medium">
                    <div className="flex items-center">
                      <div 
                        className="w-3 h-3 rounded-full mr-2" 
                        style={{ backgroundColor: colorMap[index % colorMap.length] }}
                      />
                      {exp.title}
                    </div>
                  </td>
                  <td className="px-4 py-2">{exp.results.viability}</td>
                  <td className="px-4 py-2">{exp.results.recovery}</td>
                  <td className="px-4 py-2">{exp.results.functionality}</td>
                  <td className="px-4 py-2">{exp.protocol.name}</td>
                  <td className="px-4 py-2">{exp.cellType}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </div>
    </div>
  );
}
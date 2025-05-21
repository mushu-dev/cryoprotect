import React, { useState } from 'react';
import { formatExperimentDataForCharts } from '../utils/data-transformation';

/**
 * Chart component for visualizing experiment results
 * This is a simplified version that renders bar charts with HTML/CSS
 * In a production app, you would likely use a library like Chart.js or Recharts
 */
export default function ExperimentChart({ experiment, chartType = 'bar' }) {
  const [selectedMetric, setSelectedMetric] = useState('combined');
  
  if (!experiment || !experiment.results) {
    return (
      <div className="p-6 text-center">
        <p>No data available for charting</p>
      </div>
    );
  }

  // Format data for charts
  const { viabilityData, recoveryData, functionalityData, combinedData } = 
    formatExperimentDataForCharts(experiment);

  // Get the appropriate data based on selected metric
  const getChartData = () => {
    switch (selectedMetric) {
      case 'viability':
        return viabilityData;
      case 'recovery':
        return recoveryData;
      case 'functionality':
        return functionalityData;
      case 'combined':
      default:
        return combinedData;
    }
  };

  const chartData = getChartData();
  const maxValue = 100; // Max percentage value

  return (
    <div className="bg-white rounded-lg border p-4 shadow-sm">
      <div className="flex justify-between items-center mb-4">
        <h3 className="text-lg font-semibold">Results Visualization</h3>
        <div className="flex gap-2">
          <select 
            value={selectedMetric}
            onChange={(e) => setSelectedMetric(e.target.value)}
            className="text-sm border rounded-md px-2 py-1"
          >
            <option value="combined">All Metrics</option>
            <option value="viability">Viability Only</option>
            <option value="recovery">Recovery Only</option>
            <option value="functionality">Functionality Only</option>
          </select>
          <select 
            value={chartType}
            disabled // In a real app, this would allow switching chart types
            className="text-sm border rounded-md px-2 py-1"
          >
            <option value="bar">Bar Chart</option>
            <option value="line">Line Chart</option>
            <option value="radar">Radar Chart</option>
          </select>
        </div>
      </div>

      <div className="h-64 flex items-end justify-around p-4 border-b border-l">
        {chartData.map((item, index) => {
          // Calculate bar height as percentage of container height
          const heightPercentage = (item.value / maxValue) * 100;
          
          return (
            <div key={index} className="flex flex-col items-center gap-2 w-16">
              <div 
                className="w-12 bg-blue-500 rounded-t-sm"
                style={{ 
                  height: `${heightPercentage}%`,
                  backgroundColor: 
                    item.name === 'Viability' ? '#3B82F6' : // blue
                    item.name === 'Recovery' ? '#10B981' : // green
                    '#8B5CF6' // purple for functionality
                }}
              />
              <span className="text-xs font-medium">{item.name}</span>
              <span className="text-sm font-bold">{item.value}%</span>
            </div>
          );
        })}
      </div>

      <div className="mt-4 text-center text-sm text-gray-500">
        {selectedMetric === 'combined' ? (
          <p>Comparing all performance metrics for this experiment</p>
        ) : (
          <p>Showing {selectedMetric} performance</p>
        )}
      </div>
    </div>
  );
}
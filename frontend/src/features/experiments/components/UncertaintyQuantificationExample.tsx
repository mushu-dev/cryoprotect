/**
 * Example usage of the UncertaintyQuantificationVisualizer component
 * 
 * This component demonstrates how to use the advanced uncertainty quantification
 * features with experimental data from the Convex backend.
 */

import React, { useState, useEffect } from 'react';
import { useQuery, useMutation } from 'convex/react';
import { api } from '../../../convex/_generated/api';
import { Id } from '../../../convex/_generated/dataModel';
import { UncertaintyQuantificationVisualizer, UncertaintyQuantificationMethod, UncertaintyDataPoint } from './UncertaintyQuantificationVisualizer';
import { useExperimentTimeSeries } from '../hooks/use-convex-experiments';

interface UncertaintyQuantificationExampleProps {
  experimentId: Id<"enhancedExperiments">;
}

export function UncertaintyQuantificationExample({ experimentId }: UncertaintyQuantificationExampleProps) {
  // Get experiment details
  const experiment = useQuery(
    api.experiments.enhanced_experiments.getEnhancedExperiment,
    { experimentId }
  );
  
  // Get time series for this experiment
  const { timeSeries, loading: timeSeriesLoading } = useExperimentTimeSeries(experimentId);
  
  // State for selected time series
  const [selectedTimeSeriesId, setSelectedTimeSeriesId] = useState<Id<"timeSeries"> | null>(null);
  
  // State for visualization method
  const [selectedMethod, setSelectedMethod] = useState<UncertaintyQuantificationMethod>('confidence_interval');
  
  // State for dark/light theme
  const [isDarkTheme, setIsDarkTheme] = useState(false);
  
  // State for confidence level
  const [confidenceLevel, setConfidenceLevel] = useState(0.95);
  
  // State for Monte Carlo samples
  const [monteCarloSamples, setMonteCarloSamples] = useState(100);
  
  // State for sensitivity parameters
  const [sensitivityParameters, setSensitivityParameters] = useState<string[]>(['temperature', 'concentration']);
  
  // Auto-select first time series if none selected
  useEffect(() => {
    if (!selectedTimeSeriesId && timeSeries && timeSeries.length > 0) {
      setSelectedTimeSeriesId(timeSeries[0]._id);
    }
  }, [timeSeries, selectedTimeSeriesId]);
  
  // Handle point click
  const handlePointClick = (data: UncertaintyDataPoint) => {
    console.log('Data point clicked:', data);
  };
  
  // Generate sample annotations
  const createSampleAnnotations = () => {
    if (!experiment?.date) return [];
    
    const baseTimestamp = new Date(experiment.date).getTime();
    
    return [
      {
        label: 'Critical Point',
        time: baseTimestamp + 3600000, // 1 hour after experiment start
        description: 'Significant change in experimental conditions',
        color: '#f43f5e'
      },
      {
        label: 'Anomaly Detected',
        time: baseTimestamp + 7200000, // 2 hours after experiment start
        description: 'Unexpected fluctuation in measurements',
        color: '#eab308'
      },
      {
        label: 'Steady State',
        time: baseTimestamp + 10800000, // 3 hours after experiment start
        endTime: baseTimestamp + 14400000, // 4 hours after experiment start
        description: 'Period of stable measurements',
        type: 'range' as const,
        color: '#3b82f6'
      }
    ];
  };
  
  // Handle method change
  const handleMethodChange = (method: UncertaintyQuantificationMethod) => {
    setSelectedMethod(method);
  };
  
  // Loading state
  if (timeSeriesLoading || !experiment) {
    return <div className="p-6">Loading experiment data...</div>;
  }
  
  // No time series data
  if (!timeSeries || timeSeries.length === 0) {
    return (
      <div className="p-6">
        <h2 className="text-2xl font-bold mb-4">{experiment.name}</h2>
        <p className="text-gray-700">No time series data available for uncertainty analysis.</p>
      </div>
    );
  }
  
  return (
    <div className="p-6">
      <div className="mb-6">
        <h2 className="text-2xl font-bold mb-2">{experiment.name}</h2>
        <p className="text-gray-700 mb-4">{experiment.description}</p>
        
        <div className="grid grid-cols-1 md:grid-cols-3 gap-4 mb-6">
          <div className="bg-gray-100 dark:bg-gray-800 p-4 rounded-lg">
            <h3 className="font-semibold mb-1">Status</h3>
            <p>{experiment.status}</p>
          </div>
          
          {experiment.date && (
            <div className="bg-gray-100 dark:bg-gray-800 p-4 rounded-lg">
              <h3 className="font-semibold mb-1">Date</h3>
              <p>{new Date(experiment.date).toLocaleDateString()}</p>
            </div>
          )}
          
          {experiment.temperature && (
            <div className="bg-gray-100 dark:bg-gray-800 p-4 rounded-lg">
              <h3 className="font-semibold mb-1">Temperature</h3>
              <p>{experiment.temperature} {experiment.temperatureUnit}</p>
            </div>
          )}
        </div>
      </div>
      
      <div className="mb-6 grid grid-cols-1 md:grid-cols-2 gap-6">
        {/* Time series selector */}
        <div>
          <label className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-2">
            Select time series for uncertainty analysis:
          </label>
          <select
            value={selectedTimeSeriesId?.toString() || ''}
            onChange={(e) => setSelectedTimeSeriesId(e.target.value as Id<"timeSeries">)}
            className="w-full p-2 border border-gray-300 rounded-md bg-white dark:bg-gray-800 dark:border-gray-700"
          >
            {timeSeries.map((series) => (
              <option key={series._id.toString()} value={series._id.toString()}>
                {series.name} {series.units ? `(${series.units})` : ''}
              </option>
            ))}
          </select>
        </div>
        
        {/* Method selector */}
        <div>
          <label className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-2">
            Select uncertainty quantification method:
          </label>
          <select
            value={selectedMethod}
            onChange={(e) => handleMethodChange(e.target.value as UncertaintyQuantificationMethod)}
            className="w-full p-2 border border-gray-300 rounded-md bg-white dark:bg-gray-800 dark:border-gray-700"
          >
            <option value="confidence_interval">Confidence Interval Analysis</option>
            <option value="monte_carlo">Monte Carlo Simulation</option>
            <option value="sensitivity">Sensitivity Analysis</option>
            <option value="error_propagation">Error Propagation</option>
          </select>
        </div>
      </div>
      
      {/* Visualization controls */}
      <div className="mb-6 grid grid-cols-1 md:grid-cols-3 gap-4">
        {/* Theme toggle */}
        <div className="flex items-center">
          <label className="flex items-center cursor-pointer">
            <input
              type="checkbox"
              checked={isDarkTheme}
              onChange={() => setIsDarkTheme(!isDarkTheme)}
              className="sr-only"
            />
            <div className={`relative w-10 h-5 transition-colors duration-200 ease-linear rounded-full ${isDarkTheme ? 'bg-blue-600' : 'bg-gray-300'}`}>
              <div className={`absolute left-0 w-5 h-5 transition-transform duration-200 ease-linear transform bg-white rounded-full ${isDarkTheme ? 'translate-x-5' : 'translate-x-0'}`}></div>
            </div>
            <span className="ml-2 text-sm font-medium text-gray-700 dark:text-gray-300">Dark Theme</span>
          </label>
        </div>
        
        {/* Confidence level selector */}
        {selectedMethod === 'confidence_interval' && (
          <div>
            <label className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-2">
              Confidence Level:
            </label>
            <select
              value={confidenceLevel}
              onChange={(e) => setConfidenceLevel(parseFloat(e.target.value))}
              className="w-full p-2 border border-gray-300 rounded-md bg-white dark:bg-gray-800 dark:border-gray-700"
            >
              <option value="0.99">99%</option>
              <option value="0.95">95%</option>
              <option value="0.90">90%</option>
              <option value="0.80">80%</option>
              <option value="0.50">50%</option>
            </select>
          </div>
        )}
        
        {/* Monte Carlo samples selector */}
        {selectedMethod === 'monte_carlo' && (
          <div>
            <label className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-2">
              Monte Carlo Samples:
            </label>
            <select
              value={monteCarloSamples}
              onChange={(e) => setMonteCarloSamples(parseInt(e.target.value))}
              className="w-full p-2 border border-gray-300 rounded-md bg-white dark:bg-gray-800 dark:border-gray-700"
            >
              <option value="20">20 samples</option>
              <option value="50">50 samples</option>
              <option value="100">100 samples</option>
              <option value="200">200 samples</option>
              <option value="500">500 samples</option>
            </select>
          </div>
        )}
        
        {/* Sensitivity parameters selector */}
        {selectedMethod === 'sensitivity' && (
          <div>
            <label className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-2">
              Sensitivity Parameters:
            </label>
            <div className="flex flex-wrap gap-2">
              {['temperature', 'concentration', 'pressure', 'ph', 'mixing_rate'].map((param) => (
                <label key={param} className="flex items-center">
                  <input
                    type="checkbox"
                    checked={sensitivityParameters.includes(param)}
                    onChange={() => {
                      if (sensitivityParameters.includes(param)) {
                        setSensitivityParameters(sensitivityParameters.filter(p => p !== param));
                      } else {
                        setSensitivityParameters([...sensitivityParameters, param]);
                      }
                    }}
                    className="mr-1.5"
                  />
                  <span className="text-sm text-gray-700 dark:text-gray-300">
                    {param.replace('_', ' ')}
                  </span>
                </label>
              ))}
            </div>
          </div>
        )}
      </div>
      
      {/* Uncertainty Visualization */}
      {selectedTimeSeriesId && (
        <div className="bg-white dark:bg-gray-900 rounded-lg shadow-md overflow-hidden">
          <UncertaintyQuantificationVisualizer
            timeSeriesId={selectedTimeSeriesId}
            experimentId={experimentId}
            method={selectedMethod}
            confidenceLevel={confidenceLevel}
            monteCarloSamples={monteCarloSamples}
            height={600}
            theme={isDarkTheme ? 'dark' : 'light'}
            sensitivityParameters={sensitivityParameters}
            annotations={createSampleAnnotations()}
            onPointClick={handlePointClick}
            actions={
              <button
                onClick={() => alert('Export uncertainty analysis data functionality would go here')}
                className={`px-3 py-1.5 text-sm rounded-md transition-colors ${
                  isDarkTheme
                    ? 'bg-gray-700 text-white hover:bg-gray-600'
                    : 'bg-gray-200 text-gray-800 hover:bg-gray-300'
                }`}
              >
                Export Analysis
              </button>
            }
          />
        </div>
      )}
      
      {/* Method details */}
      <div className="mt-6 bg-gray-100 dark:bg-gray-800 p-4 rounded-lg">
        <h3 className="text-lg font-semibold mb-2">About Uncertainty Quantification</h3>
        <p className="text-gray-700 dark:text-gray-300 mb-4">
          Uncertainty quantification helps researchers understand the reliability and precision of their experimental results. 
          It provides statistical tools to assess how variations in measurements, parameters, and conditions affect the conclusions drawn from the data.
        </p>
        
        <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
          <div>
            <h4 className="font-medium mb-1">Benefits:</h4>
            <ul className="list-disc list-inside text-sm text-gray-700 dark:text-gray-300 space-y-1">
              <li>Improves result reliability assessment</li>
              <li>Highlights potential sources of error</li>
              <li>Supports data-driven decision making</li>
              <li>Enhances reproducibility of experiments</li>
              <li>Guides experimental design improvements</li>
            </ul>
          </div>
          
          <div>
            <h4 className="font-medium mb-1">Scientific Applications:</h4>
            <ul className="list-disc list-inside text-sm text-gray-700 dark:text-gray-300 space-y-1">
              <li>Cryoprotectant performance evaluation</li>
              <li>Cell viability assessments</li>
              <li>Freezing protocol optimization</li>
              <li>Material property characterization</li>
              <li>Multi-factor experimental analysis</li>
            </ul>
          </div>
        </div>
      </div>
    </div>
  );
}
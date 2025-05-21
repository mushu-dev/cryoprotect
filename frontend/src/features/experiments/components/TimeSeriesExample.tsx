/**
 * Example usage of the TimeSeriesVisualizer component
 * 
 * This component demonstrates how to use the TimeSeriesVisualizer
 * with real data from the Convex backend.
 */

import React, { useState } from 'react';
import { useQuery, useMutation } from 'convex/react';
import { api } from '../../../convex/_generated/api';
import { Id } from '../../../convex/_generated/dataModel';
import { TimeSeriesVisualizer, TimeSeriesData } from './TimeSeriesVisualizer';
import { useExperimentTimeSeries } from '../hooks/use-convex-experiments';

interface TimeSeriesExampleProps {
  experimentId: Id<"enhancedExperiments">;
}

export function TimeSeriesExample({ experimentId }: TimeSeriesExampleProps) {
  // Get experiment details
  const { experiment, loading: experimentLoading } = useQuery(
    api.experiments.enhanced_experiments.getEnhancedExperiment,
    { experimentId }
  );
  
  // Get time series for this experiment
  const { timeSeries, loading: timeSeriesLoading } = useExperimentTimeSeries(experimentId);
  
  // State for selected time series
  const [selectedTimeSeriesId, setSelectedTimeSeriesId] = useState<Id<"timeSeries"> | null>(null);
  
  // State for dark/light theme
  const [isDarkTheme, setIsDarkTheme] = useState(false);
  
  // State for showing raw data
  const [showRawData, setShowRawData] = useState(false);
  
  // State for showing uncertainty
  const [showUncertainty, setShowUncertainty] = useState(true);
  
  // Handle point click
  const handlePointClick = (data: TimeSeriesData) => {
    console.log('Time series data point clicked:', data);
    alert(`Clicked data point: ${new Date(data.timestamp).toLocaleString()}, Value: ${data.value}`);
  };
  
  // Create a sample annotation
  const createSampleAnnotations = (seriesId: Id<"timeSeries">) => {
    const timeSeriesData = useQuery(
      api.experiments.enhanced_index.getTimeSeriesData,
      { timeSeriesId: seriesId }
    );
    
    if (!timeSeriesData || timeSeriesData.length === 0) {
      return [];
    }
    
    // Get a point in the middle of the data
    const midpointIndex = Math.floor(timeSeriesData.length / 2);
    const midTimestamp = timeSeriesData[midpointIndex].timestamp;
    
    // Get a range near the end of the data
    const rangeStartIndex = Math.floor(timeSeriesData.length * 0.7);
    const rangeEndIndex = Math.floor(timeSeriesData.length * 0.8);
    const rangeStartTimestamp = timeSeriesData[rangeStartIndex].timestamp;
    const rangeEndTimestamp = timeSeriesData[rangeEndIndex].timestamp;
    
    return [
      {
        label: 'Important Event',
        time: midTimestamp,
        description: 'This is a significant moment in the experiment',
        color: '#f43f5e'
      },
      {
        label: 'Cooling Phase',
        time: rangeStartTimestamp,
        endTime: rangeEndTimestamp,
        description: 'Period of controlled cooling during the experiment',
        type: 'range' as const,
        color: '#3b82f6'
      }
    ];
  };
  
  // Loading state
  if (experimentLoading || timeSeriesLoading) {
    return <div>Loading experiment data...</div>;
  }
  
  // No experiment data
  if (!experiment) {
    return <div>Experiment not found</div>;
  }
  
  // No time series data
  if (!timeSeries || timeSeries.length === 0) {
    return (
      <div>
        <h2>{experiment.name}</h2>
        <p>No time series data available for this experiment.</p>
      </div>
    );
  }
  
  // Auto-select first time series if none selected
  if (!selectedTimeSeriesId && timeSeries.length > 0) {
    setSelectedTimeSeriesId(timeSeries[0]._id);
  }
  
  return (
    <div style={{ padding: '20px' }}>
      <h2>{experiment.name}</h2>
      
      {/* Experiment metadata */}
      <div style={{ marginBottom: '20px' }}>
        <p>{experiment.description}</p>
        <p><strong>Status:</strong> {experiment.status}</p>
        {experiment.date && (
          <p><strong>Date:</strong> {new Date(experiment.date).toLocaleDateString()}</p>
        )}
      </div>
      
      {/* Time series selector */}
      <div style={{ marginBottom: '20px' }}>
        <label style={{ marginRight: '10px' }}>Select time series:</label>
        <select
          value={selectedTimeSeriesId?.toString() || ''}
          onChange={(e) => setSelectedTimeSeriesId(e.target.value as Id<"timeSeries">)}
          style={{ padding: '8px', borderRadius: '4px', minWidth: '200px' }}
        >
          {timeSeries.map((series) => (
            <option key={series._id.toString()} value={series._id.toString()}>
              {series.name} {series.units ? `(${series.units})` : ''}
            </option>
          ))}
        </select>
      </div>
      
      {/* Visualization controls */}
      <div style={{ 
        marginBottom: '20px', 
        display: 'flex', 
        gap: '15px',
        flexWrap: 'wrap'
      }}>
        <label style={{ display: 'flex', alignItems: 'center', cursor: 'pointer' }}>
          <input
            type="checkbox"
            checked={isDarkTheme}
            onChange={() => setIsDarkTheme(!isDarkTheme)}
            style={{ marginRight: '5px' }}
          />
          Dark Theme
        </label>
        
        <label style={{ display: 'flex', alignItems: 'center', cursor: 'pointer' }}>
          <input
            type="checkbox"
            checked={showUncertainty}
            onChange={() => setShowUncertainty(!showUncertainty)}
            style={{ marginRight: '5px' }}
          />
          Show Uncertainty
        </label>
        
        <label style={{ display: 'flex', alignItems: 'center', cursor: 'pointer' }}>
          <input
            type="checkbox"
            checked={showRawData}
            onChange={() => setShowRawData(!showRawData)}
            style={{ marginRight: '5px' }}
          />
          Show Raw Data
        </label>
      </div>
      
      {/* Time series visualizer */}
      {selectedTimeSeriesId && (
        <TimeSeriesVisualizer
          timeSeriesId={selectedTimeSeriesId}
          showUncertainty={showUncertainty}
          showRaw={showRawData}
          height={500}
          theme={isDarkTheme ? 'dark' : 'light'}
          annotations={createSampleAnnotations(selectedTimeSeriesId)}
          onPointClick={handlePointClick}
          compareWith={timeSeries
            .filter(series => series._id !== selectedTimeSeriesId)
            .slice(0, 2) // Limit to 2 comparison series for clarity
            .map((series, index) => ({
              id: series._id,
              label: series.name
            }))}
          actions={
            <button
              onClick={() => alert('Export functionality would go here')}
              style={{
                backgroundColor: isDarkTheme ? '#4b5563' : '#e5e7eb',
                color: isDarkTheme ? '#f9fafb' : '#111827',
                border: 'none',
                padding: '5px 10px',
                borderRadius: '4px',
                cursor: 'pointer'
              }}
            >
              Export Data
            </button>
          }
        />
      )}
      
      {/* Additional experiment information */}
      <div style={{ marginTop: '30px' }}>
        <h3>Experiment Details</h3>
        <div style={{ 
          display: 'grid', 
          gridTemplateColumns: 'repeat(auto-fill, minmax(300px, 1fr))',
          gap: '15px',
          marginTop: '10px'
        }}>
          <div>
            <strong>Temperature:</strong> {experiment.temperature} {experiment.temperatureUnit}
          </div>
          {experiment.coolingRate && (
            <div>
              <strong>Cooling Rate:</strong> {experiment.coolingRate} {experiment.coolingRateUnit}
            </div>
          )}
          {experiment.pressure && (
            <div>
              <strong>Pressure:</strong> {experiment.pressure} {experiment.pressureUnit}
            </div>
          )}
          {experiment.tags && experiment.tags.length > 0 && (
            <div>
              <strong>Tags:</strong> {experiment.tags.join(', ')}
            </div>
          )}
        </div>
      </div>
    </div>
  );
}
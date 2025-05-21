/**
 * Uncertainty Quantification Visualizer Component
 * 
 * This component provides advanced uncertainty visualization and analysis tools
 * for experimental data, including confidence intervals, error propagation,
 * sensitivity analysis, and Monte Carlo simulations.
 */

import React, { useState, useCallback, useMemo } from 'react';
import { useQuery, useMutation } from 'convex/react';
import { api } from '../../../convex/_generated/api';
import { Id } from '../../../convex/_generated/dataModel';
import {
  LineChart,
  Line,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  Legend,
  ResponsiveContainer,
  ReferenceArea,
  ReferenceLine,
  Area,
  Scatter,
  ScatterChart,
  Brush,
  ComposedChart
} from 'recharts';

// Types
export interface UncertaintyDataPoint {
  timestamp: number;
  value: number;
  standardDeviation?: number;
  confidenceInterval?: [number, number]; // [lower, upper]
  monteCarlo?: Array<number>;
  sensitivityScores?: Record<string, number>;
  metadata?: Record<string, any>;
}

export type UncertaintyQuantificationMethod = 'confidence_interval' | 'monte_carlo' | 'sensitivity' | 'error_propagation';

export interface UncertaintyQuantificationVisualizerProps {
  timeSeriesId: Id<"timeSeries">;
  experimentId?: Id<"enhancedExperiments">;
  method?: UncertaintyQuantificationMethod;
  confidenceLevel?: number; // between 0 and 1, e.g., 0.95 for 95% confidence interval
  monteCarloSamples?: number; // number of Monte Carlo samples to display
  height?: number;
  width?: number;
  theme?: 'light' | 'dark';
  sensitivityParameters?: string[]; // parameters to show sensitivity for
  annotations?: Array<{
    label: string;
    time: number;
    description?: string;
    type?: 'point' | 'range';
    endTime?: number;
    color?: string;
  }>;
  onPointClick?: (data: UncertaintyDataPoint) => void;
  actions?: React.ReactNode;
}

/**
 * Uncertainty Quantification Visualizer Component
 */
export function UncertaintyQuantificationVisualizer({
  timeSeriesId,
  experimentId,
  method = 'confidence_interval',
  confidenceLevel = 0.95,
  monteCarloSamples = 100,
  height = 500,
  width = '100%',
  theme = 'light',
  sensitivityParameters = [],
  annotations = [],
  onPointClick,
  actions
}: UncertaintyQuantificationVisualizerProps) {
  // State for zoom
  const [zoomDomain, setZoomDomain] = useState<[number, number] | null>(null);
  const [isZooming, setIsZooming] = useState(false);
  const [zoomStart, setZoomStart] = useState<number | null>(null);
  
  // State for visualization method
  const [selectedMethod, setSelectedMethod] = useState<UncertaintyQuantificationMethod>(method);
  
  // State for confidence level
  const [selectedConfidenceLevel, setSelectedConfidenceLevel] = useState(confidenceLevel);
  
  // State for Monte Carlo samples to display
  const [selectedMonteCarloSamples, setSelectedMonteCarloSamples] = useState(monteCarloSamples);
  
  // State for selected sensitivity parameters
  const [selectedSensitivityParameters, setSelectedSensitivityParameters] = useState<string[]>(sensitivityParameters);
  
  // Fetch the time series metadata
  const timeSeries = useQuery(api.experiments.enhanced_index.getTimeSeries, { timeSeriesId });
  
  // Fetch the time series data
  const rawTimeSeriesData = useQuery(api.experiments.enhanced_index.getTimeSeriesData, { timeSeriesId });
  
  // Fetch the experiment data if experimentId is provided
  const experiment = useQuery(api.experiments.enhanced_experiments.getEnhancedExperiment, 
    experimentId ? { experimentId } : "skip");
  
  // Fetch uncertainty analysis data
  const uncertaintyData = useQuery(api.experiments.uncertainty.getUncertaintyAnalysis, 
    { timeSeriesId, method: selectedMethod });

  // Calculate confidence intervals if needed and not already provided
  const processedData = useMemo(() => {
    if (!rawTimeSeriesData) return [];
    
    // If uncertainty data is available, use it
    if (uncertaintyData) {
      return uncertaintyData;
    }
    
    // Otherwise, perform basic uncertainty analysis based on the selected method
    return rawTimeSeriesData.map(point => {
      const dataPoint: UncertaintyDataPoint = {
        timestamp: point.timestamp,
        value: point.value,
        standardDeviation: point.uncertainty,
        metadata: point.metadata
      };
      
      // If standard deviation/uncertainty is available, calculate confidence interval
      if (point.uncertainty) {
        const zScore = getZScoreForConfidenceLevel(selectedConfidenceLevel);
        dataPoint.confidenceInterval = [
          point.value - zScore * point.uncertainty,
          point.value + zScore * point.uncertainty
        ];
        
        // Generate simple Monte Carlo samples if that method is selected
        if (selectedMethod === 'monte_carlo') {
          dataPoint.monteCarlo = Array.from({ length: selectedMonteCarloSamples }, () => 
            generateNormalRandom(point.value, point.uncertainty));
        }
      }
      
      return dataPoint;
    });
  }, [rawTimeSeriesData, uncertaintyData, selectedMethod, selectedConfidenceLevel, selectedMonteCarloSamples]);
  
  // Prepare the chart data with formatted dates and times
  const chartData = useMemo(() => {
    if (!processedData) return [];
    
    return processedData.map(dataPoint => ({
      ...dataPoint,
      time: new Date(dataPoint.timestamp).toLocaleTimeString(),
      date: new Date(dataPoint.timestamp).toLocaleDateString(),
    }));
  }, [processedData]);
  
  // Get available sensitivity parameters if any
  const availableSensitivityParameters = useMemo(() => {
    if (!processedData || processedData.length === 0) return [];
    
    const firstPoint = processedData[0];
    if (!firstPoint.sensitivityScores) return [];
    
    return Object.keys(firstPoint.sensitivityScores);
  }, [processedData]);
  
  // Calculate sensitivity chart data
  const sensitivityChartData = useMemo(() => {
    if (
      selectedMethod !== 'sensitivity' || 
      !processedData || 
      processedData.length === 0 || 
      selectedSensitivityParameters.length === 0
    ) {
      return [];
    }
    
    // Find the average sensitivity score for each parameter
    const paramScores: Record<string, number> = {};
    let count = 0;
    
    processedData.forEach(point => {
      if (point.sensitivityScores) {
        count++;
        selectedSensitivityParameters.forEach(param => {
          if (param in point.sensitivityScores!) {
            paramScores[param] = (paramScores[param] || 0) + point.sensitivityScores![param];
          }
        });
      }
    });
    
    // Calculate averages and convert to chart format
    return Object.entries(paramScores).map(([param, score]) => ({
      parameter: param,
      score: score / count,
      normalized: (score / count) / Math.max(...Object.values(paramScores).map(s => s / count))
    }));
  }, [processedData, selectedMethod, selectedSensitivityParameters]);
  
  // Get the domain for the time axis
  const domain = useMemo(() => {
    if (chartData.length === 0) return [0, 1]; // Default domain
    
    const timestamps = chartData.map(data => data.timestamp);
    const minTime = Math.min(...timestamps);
    const maxTime = Math.max(...timestamps);
    
    // If zoomed, use zoom domain
    if (zoomDomain) {
      return zoomDomain;
    }
    
    // Add some padding to the domain
    const padding = (maxTime - minTime) * 0.05;
    return [minTime - padding, maxTime + padding];
  }, [chartData, zoomDomain]);
  
  // Handle mouse down for zooming
  const handleMouseDown = useCallback((e: any) => {
    if (!e || !e.activeLabel) return;
    
    setIsZooming(true);
    setZoomStart(Number(e.activeLabel));
  }, []);
  
  // Handle mouse up for zooming
  const handleMouseUp = useCallback((e: any) => {
    if (!isZooming || !zoomStart || !e || !e.activeLabel) {
      setIsZooming(false);
      setZoomStart(null);
      return;
    }
    
    const zoomEnd = Number(e.activeLabel);
    
    // Ensure start < end
    const start = Math.min(zoomStart, zoomEnd);
    const end = Math.max(zoomStart, zoomEnd);
    
    // Only zoom if range is meaningful
    if (end - start > 1000) { // 1 second minimum zoom
      setZoomDomain([start, end]);
    }
    
    setIsZooming(false);
    setZoomStart(null);
  }, [isZooming, zoomStart]);
  
  // Handle reset zoom
  const handleResetZoom = useCallback(() => {
    setZoomDomain(null);
  }, []);
  
  // Handle click on a data point
  const handleClick = useCallback((data: any) => {
    if (onPointClick && data && data.activePayload) {
      onPointClick(data.activePayload[0].payload);
    }
  }, [onPointClick]);
  
  // Get colors based on theme
  const colors = useMemo(() => {
    return theme === 'dark' 
      ? { 
          background: '#1f2937', 
          text: '#f9fafb',
          grid: '#374151',
          primary: '#60a5fa',
          uncertainty: 'rgba(96, 165, 250, 0.3)',
          confidenceArea: 'rgba(96, 165, 250, 0.2)',
          monteCarloPoints: 'rgba(249, 115, 22, 0.4)',
          sensitivityBars: [
            'rgba(96, 165, 250, 0.8)',
            'rgba(249, 115, 22, 0.8)',
            'rgba(16, 185, 129, 0.8)',
            'rgba(217, 70, 239, 0.8)',
            'rgba(248, 113, 113, 0.8)'
          ]
        }
      : {
          background: '#ffffff',
          text: '#111827',
          grid: '#e5e7eb',
          primary: '#3b82f6',
          uncertainty: 'rgba(59, 130, 246, 0.3)',
          confidenceArea: 'rgba(59, 130, 246, 0.1)',
          monteCarloPoints: 'rgba(234, 88, 12, 0.3)',
          sensitivityBars: [
            'rgba(59, 130, 246, 0.8)',
            'rgba(234, 88, 12, 0.8)',
            'rgba(5, 150, 105, 0.8)',
            'rgba(192, 38, 211, 0.8)',
            'rgba(239, 68, 68, 0.8)'
          ]
        };
  }, [theme]);
  
  // Handle change of visualization method
  const handleMethodChange = (method: UncertaintyQuantificationMethod) => {
    setSelectedMethod(method);
  };
  
  // Handle change of confidence level
  const handleConfidenceLevelChange = (level: number) => {
    setSelectedConfidenceLevel(level);
  };
  
  // Handle sensitivity parameter selection
  const handleSensitivityParameterToggle = (parameter: string) => {
    if (selectedSensitivityParameters.includes(parameter)) {
      setSelectedSensitivityParameters(
        selectedSensitivityParameters.filter(p => p !== parameter)
      );
    } else {
      setSelectedSensitivityParameters([...selectedSensitivityParameters, parameter]);
    }
  };
  
  // Custom tooltip
  const CustomTooltip = useCallback(({ active, payload, label }: any) => {
    if (!active || !payload || !payload.length) return null;
    
    const data = payload[0].payload;
    
    return (
      <div style={{
        backgroundColor: colors.background,
        border: `1px solid ${colors.grid}`,
        padding: '10px',
        borderRadius: '5px'
      }}>
        <p style={{ margin: '0 0 5px', fontWeight: 'bold', color: colors.text }}>
          {data.date} {data.time}
        </p>
        
        <p style={{ margin: '0', color: colors.text }}>
          Value: {data.value.toFixed(3)} {timeSeries?.units || ''}
        </p>
        
        {data.standardDeviation && (
          <p style={{ margin: '0', color: colors.text }}>
            Standard Deviation: ±{data.standardDeviation.toFixed(3)}
          </p>
        )}
        
        {data.confidenceInterval && (
          <p style={{ margin: '0', color: colors.text }}>
            {(selectedConfidenceLevel * 100).toFixed(0)}% CI: [{data.confidenceInterval[0].toFixed(3)}, {data.confidenceInterval[1].toFixed(3)}]
          </p>
        )}
        
        {data.sensitivityScores && selectedSensitivityParameters.length > 0 && (
          <div style={{ marginTop: '5px' }}>
            <p style={{ margin: '0', fontWeight: 'bold', color: colors.text }}>Sensitivity Scores:</p>
            {selectedSensitivityParameters.map(param => (
              data.sensitivityScores![param] !== undefined && (
                <p key={param} style={{ margin: '0', fontSize: '0.9em', color: colors.text }}>
                  {param}: {data.sensitivityScores![param].toFixed(3)}
                </p>
              )
            ))}
          </div>
        )}
        
        {/* Show any additional metadata */}
        {Object.entries(data.metadata || {})
          .filter(([key]) => !['timestamp', 'time', 'date', 'value', 'uncertainty', 'standardDeviation'].includes(key))
          .map(([key, value]) => (
            <p key={key} style={{ margin: '0', fontSize: '0.9em', color: colors.text, opacity: 0.8 }}>
              {key}: {String(value)}
            </p>
          ))
        }
      </div>
    );
  }, [colors, timeSeries, selectedConfidenceLevel, selectedSensitivityParameters]);
  
  // Render sensitivity analysis chart
  const renderSensitivityChart = useCallback(() => {
    if (selectedMethod !== 'sensitivity' || sensitivityChartData.length === 0) {
      return null;
    }
    
    return (
      <div style={{ marginTop: '20px' }}>
        <h4 style={{ margin: '0 0 10px', fontSize: '1em', color: colors.text }}>
          Parameter Sensitivity Analysis
        </h4>
        <ResponsiveContainer width="100%" height={200}>
          <ComposedChart
            data={sensitivityChartData}
            margin={{ top: 10, right: 30, left: 10, bottom: 50 }}
          >
            <CartesianGrid strokeDasharray="3 3" stroke={colors.grid} />
            <XAxis 
              dataKey="parameter" 
              tick={{ fill: colors.text }}
              stroke={colors.text}
              angle={-45}
              textAnchor="end"
              height={60}
            />
            <YAxis 
              label={{ 
                value: 'Sensitivity Score', 
                angle: -90, 
                position: 'insideLeft',
                style: { fill: colors.text }
              }}
              tick={{ fill: colors.text }}
              stroke={colors.text}
            />
            <Tooltip 
              content={({ active, payload }) => {
                if (!active || !payload || !payload.length) return null;
                
                const data = payload[0].payload;
                
                return (
                  <div style={{
                    backgroundColor: colors.background,
                    border: `1px solid ${colors.grid}`,
                    padding: '10px',
                    borderRadius: '5px'
                  }}>
                    <p style={{ margin: '0 0 5px', fontWeight: 'bold', color: colors.text }}>
                      {data.parameter}
                    </p>
                    <p style={{ margin: '0', color: colors.text }}>
                      Sensitivity: {data.score.toFixed(3)}
                    </p>
                    <p style={{ margin: '0', color: colors.text }}>
                      Normalized: {data.normalized.toFixed(2)}
                    </p>
                  </div>
                );
              }}
            />
            <Legend />
            {sensitivityChartData.map((entry, index) => (
              <Area
                key={entry.parameter}
                type="monotone"
                dataKey="score"
                name={entry.parameter}
                fill={colors.sensitivityBars[index % colors.sensitivityBars.length]}
                stroke={colors.sensitivityBars[index % colors.sensitivityBars.length]}
                fillOpacity={0.8}
              />
            ))}
          </ComposedChart>
        </ResponsiveContainer>
      </div>
    );
  }, [colors, selectedMethod, sensitivityChartData]);
  
  // Render Monte Carlo scatter plot
  const renderMonteCarloScatter = useCallback(() => {
    if (selectedMethod !== 'monte_carlo' || chartData.length === 0) {
      return null;
    }
    
    // Flatten all Monte Carlo samples
    const scatterData = chartData.flatMap(point => {
      if (!point.monteCarlo) return [];
      
      return point.monteCarlo.map((value, idx) => ({
        timestamp: point.timestamp,
        value,
        sampleId: idx
      }));
    });
    
    return (
      <ScatterChart
        width={0}
        height={0}
        margin={{ top: 20, right: 20, bottom: 20, left: 20 }}
      >
        <CartesianGrid strokeDasharray="3 3" stroke={colors.grid} />
        <XAxis
          type="number"
          dataKey="timestamp"
          name="Time"
          domain={domain}
          tick={{ fill: colors.text }}
          stroke={colors.text}
          hide
        />
        <YAxis
          type="number"
          dataKey="value"
          name="Value"
          tick={{ fill: colors.text }}
          stroke={colors.text}
          hide
        />
        <Scatter
          name="Monte Carlo Samples"
          data={scatterData}
          fill={colors.monteCarloPoints}
          shape="circle"
          isAnimationActive={false}
        />
      </ScatterChart>
    );
  }, [colors, domain, selectedMethod, chartData]);
  
  // If data is loading
  if (!timeSeries || !rawTimeSeriesData) {
    return (
      <div style={{ 
        height, 
        width, 
        display: 'flex', 
        alignItems: 'center', 
        justifyContent: 'center',
        backgroundColor: colors.background,
        color: colors.text
      }}>
        Loading uncertainty data...
      </div>
    );
  }
  
  // If no data
  if (rawTimeSeriesData.length === 0) {
    return (
      <div style={{ 
        height, 
        width, 
        display: 'flex', 
        alignItems: 'center', 
        justifyContent: 'center',
        backgroundColor: colors.background,
        color: colors.text,
        flexDirection: 'column',
        gap: '10px'
      }}>
        <p>No data available for uncertainty analysis.</p>
        {actions && <div>{actions}</div>}
      </div>
    );
  }
  
  return (
    <div style={{ 
      width, 
      backgroundColor: colors.background,
      color: colors.text,
      padding: '20px',
      borderRadius: '8px'
    }}>
      <div style={{ 
        display: 'flex', 
        justifyContent: 'space-between',
        alignItems: 'center',
        marginBottom: '15px'
      }}>
        <div>
          <h3 style={{ margin: '0 0 5px', fontSize: '1.2em' }}>
            {timeSeries.name} Uncertainty Analysis
            {timeSeries.units && <span style={{ fontSize: '0.8em', opacity: 0.7 }}> ({timeSeries.units})</span>}
          </h3>
          {timeSeries.description && (
            <p style={{ margin: '0', fontSize: '0.9em', opacity: 0.8 }}>
              {timeSeries.description}
            </p>
          )}
        </div>
        
        <div style={{ display: 'flex', gap: '10px', alignItems: 'center' }}>
          {zoomDomain && (
            <button 
              onClick={handleResetZoom}
              style={{
                backgroundColor: colors.primary,
                color: 'white',
                border: 'none',
                padding: '5px 10px',
                borderRadius: '4px',
                cursor: 'pointer',
                fontSize: '0.9em'
              }}
            >
              Reset Zoom
            </button>
          )}
          {actions}
        </div>
      </div>
      
      {/* Method selector */}
      <div style={{ 
        display: 'flex', 
        justifyContent: 'flex-start', 
        gap: '15px', 
        flexWrap: 'wrap',
        marginBottom: '15px',
        padding: '10px',
        backgroundColor: theme === 'dark' ? '#374151' : '#f3f4f6',
        borderRadius: '6px'
      }}>
        <div>
          <label style={{ fontSize: '0.9em', marginRight: '8px', color: colors.text }}>Method:</label>
          <select
            value={selectedMethod}
            onChange={(e) => handleMethodChange(e.target.value as UncertaintyQuantificationMethod)}
            style={{
              padding: '6px',
              backgroundColor: colors.background,
              color: colors.text,
              border: `1px solid ${colors.grid}`,
              borderRadius: '4px'
            }}
          >
            <option value="confidence_interval">Confidence Interval</option>
            <option value="monte_carlo">Monte Carlo Simulation</option>
            <option value="sensitivity">Sensitivity Analysis</option>
            <option value="error_propagation">Error Propagation</option>
          </select>
        </div>
        
        {/* Confidence level selector (only for confidence interval) */}
        {selectedMethod === 'confidence_interval' && (
          <div>
            <label style={{ fontSize: '0.9em', marginRight: '8px', color: colors.text }}>Confidence Level:</label>
            <select
              value={selectedConfidenceLevel}
              onChange={(e) => handleConfidenceLevelChange(parseFloat(e.target.value))}
              style={{
                padding: '6px',
                backgroundColor: colors.background,
                color: colors.text,
                border: `1px solid ${colors.grid}`,
                borderRadius: '4px'
              }}
            >
              <option value="0.99">99%</option>
              <option value="0.95">95%</option>
              <option value="0.90">90%</option>
              <option value="0.80">80%</option>
              <option value="0.50">50%</option>
            </select>
          </div>
        )}
        
        {/* Monte Carlo sample count selector */}
        {selectedMethod === 'monte_carlo' && (
          <div>
            <label style={{ fontSize: '0.9em', marginRight: '8px', color: colors.text }}>Samples:</label>
            <select
              value={selectedMonteCarloSamples}
              onChange={(e) => setSelectedMonteCarloSamples(parseInt(e.target.value))}
              style={{
                padding: '6px',
                backgroundColor: colors.background,
                color: colors.text,
                border: `1px solid ${colors.grid}`,
                borderRadius: '4px'
              }}
            >
              <option value="20">20</option>
              <option value="50">50</option>
              <option value="100">100</option>
              <option value="200">200</option>
              <option value="500">500</option>
            </select>
          </div>
        )}
        
        {/* Sensitivity parameters selector */}
        {selectedMethod === 'sensitivity' && availableSensitivityParameters.length > 0 && (
          <div style={{ display: 'flex', flexWrap: 'wrap', gap: '6px' }}>
            <label style={{ fontSize: '0.9em', marginRight: '8px', color: colors.text }}>Parameters:</label>
            {availableSensitivityParameters.map((param) => (
              <label key={param} style={{ 
                display: 'inline-flex', 
                alignItems: 'center', 
                backgroundColor: selectedSensitivityParameters.includes(param) 
                  ? (theme === 'dark' ? '#4b5563' : '#e5e7eb') 
                  : 'transparent',
                padding: '3px 6px',
                borderRadius: '4px',
                fontSize: '0.85em',
                cursor: 'pointer'
              }}>
                <input 
                  type="checkbox"
                  checked={selectedSensitivityParameters.includes(param)}
                  onChange={() => handleSensitivityParameterToggle(param)}
                  style={{ marginRight: '4px' }}
                />
                {param}
              </label>
            ))}
          </div>
        )}
      </div>
      
      {/* Main visualization */}
      <div style={{ height: 'auto', minHeight: height - 150 }}>
        <ResponsiveContainer width="100%" height={height - 150}>
          <ComposedChart
            data={chartData}
            margin={{ top: 10, right: 30, left: 10, bottom: 30 }}
            onMouseDown={handleMouseDown}
            onMouseUp={handleMouseUp}
            onClick={handleClick}
          >
            <CartesianGrid strokeDasharray="3 3" stroke={colors.grid} />
            <XAxis 
              dataKey="timestamp" 
              type="number"
              domain={domain}
              tickFormatter={(timestamp) => new Date(timestamp).toLocaleTimeString()}
              tick={{ fill: colors.text }}
              stroke={colors.text}
            />
            <YAxis 
              tickFormatter={(value) => value.toFixed(1)}
              tick={{ fill: colors.text }}
              stroke={colors.text}
              label={{ 
                value: timeSeries.units || '', 
                angle: -90, 
                position: 'insideLeft',
                style: { fill: colors.text }
              }}
            />
            <Tooltip content={<CustomTooltip />} />
            <Legend wrapperStyle={{ color: colors.text }} />
            
            {/* Confidence interval area */}
            {selectedMethod === 'confidence_interval' && (
              <Area
                type="monotone"
                dataKey={(data) => data.confidenceInterval?.[0]}
                stroke="transparent"
                fill={colors.confidenceArea}
                activeDot={false}
                name={`${(selectedConfidenceLevel * 100).toFixed(0)}% Confidence Interval`}
                isAnimationActive={false}
              />
            )}
            {selectedMethod === 'confidence_interval' && (
              <Area
                type="monotone"
                dataKey={(data) => data.confidenceInterval?.[1]}
                stroke="transparent"
                fillOpacity={0}
                activeDot={false}
                isAnimationActive={false}
              />
            )}
            
            {/* Standard deviation area */}
            {selectedMethod === 'error_propagation' && (
              <Area
                type="monotone"
                dataKey={(data) => data.value - (data.standardDeviation || 0)}
                stroke="transparent"
                fill={colors.confidenceArea}
                activeDot={false}
                name="Standard Deviation"
                isAnimationActive={false}
              />
            )}
            {selectedMethod === 'error_propagation' && (
              <Area
                type="monotone"
                dataKey={(data) => data.value + (data.standardDeviation || 0)}
                stroke="transparent"
                fillOpacity={0}
                activeDot={false}
                isAnimationActive={false}
              />
            )}
            
            {/* Main data line */}
            <Line
              type="monotone"
              dataKey="value"
              stroke={colors.primary}
              strokeWidth={2}
              dot={{ r: 3, strokeWidth: 1 }}
              activeDot={{ r: 6, strokeWidth: 1 }}
              name="Measured Value"
              isAnimationActive={true}
            />
            
            {/* Monte Carlo scatter points */}
            {selectedMethod === 'monte_carlo' && renderMonteCarloScatter()}
            
            {/* Annotations - reference lines */}
            {annotations
              .filter(annotation => annotation.type !== 'range')
              .map((annotation, index) => (
                <ReferenceLine
                  key={`line-${index}`}
                  x={annotation.time}
                  stroke={annotation.color || '#ff7300'}
                  label={{
                    value: annotation.label,
                    position: 'top',
                    fill: annotation.color || '#ff7300',
                  }}
                />
              ))}
            
            {/* Annotations - reference areas */}
            {annotations
              .filter(annotation => annotation.type === 'range' && annotation.endTime)
              .map((annotation, index) => (
                <ReferenceArea
                  key={`area-${index}`}
                  x1={annotation.time}
                  x2={annotation.endTime}
                  stroke={annotation.color || '#ff7300'}
                  strokeOpacity={0.3}
                  fill={annotation.color || '#ff7300'}
                  fillOpacity={0.1}
                  label={{
                    value: annotation.label,
                    position: 'insideTop',
                    fill: annotation.color || '#ff7300',
                  }}
                />
              ))}
            
            {/* Zoom reference area */}
            {isZooming && zoomStart && (
              <ReferenceArea
                x1={zoomStart}
                x2={zoomStart} // Will be updated on mouse move
                strokeOpacity={0.3}
                fill={colors.primary}
                fillOpacity={0.1}
              />
            )}
            
            {/* Time brush for navigation */}
            <Brush 
              dataKey="timestamp" 
              height={30} 
              stroke={colors.primary}
              tickFormatter={(timestamp) => new Date(timestamp).toLocaleTimeString()}
              travellerWidth={10}
              startIndex={0}
              endIndex={chartData.length - 1}
            />
          </ComposedChart>
        </ResponsiveContainer>
        
        {/* Sensitivity chart */}
        {renderSensitivityChart()}
      </div>
      
      {/* Method explanation */}
      <div style={{ 
        marginTop: '20px',
        padding: '15px',
        backgroundColor: theme === 'dark' ? '#374151' : '#f3f4f6',
        borderRadius: '6px',
        fontSize: '0.9em'
      }}>
        <h4 style={{ margin: '0 0 8px', fontSize: '1em' }}>
          About {getMethodDisplayName(selectedMethod)}
        </h4>
        <p style={{ margin: '0', lineHeight: '1.5' }}>
          {getMethodDescription(selectedMethod)}
        </p>
      </div>
    </div>
  );
}

/**
 * Generate a random number from a normal distribution
 */
function generateNormalRandom(mean: number, stdDev: number): number {
  // Box-Muller transform for normal distribution
  const u1 = Math.random();
  const u2 = Math.random();
  const z0 = Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(2.0 * Math.PI * u2);
  return mean + stdDev * z0;
}

/**
 * Get the Z-score for a given confidence level
 */
function getZScoreForConfidenceLevel(confidenceLevel: number): number {
  // Common Z-scores for confidence levels
  const zScoreMap: Record<number, number> = {
    0.50: 0.674,
    0.80: 1.282,
    0.90: 1.645,
    0.95: 1.960,
    0.99: 2.576
  };
  
  return zScoreMap[confidenceLevel] || 1.960; // Default to 95% confidence
}

/**
 * Get display name for a method
 */
function getMethodDisplayName(method: UncertaintyQuantificationMethod): string {
  switch (method) {
    case 'confidence_interval':
      return 'Confidence Interval Analysis';
    case 'monte_carlo':
      return 'Monte Carlo Simulation';
    case 'sensitivity':
      return 'Sensitivity Analysis';
    case 'error_propagation':
      return 'Error Propagation';
    default:
      return 'Uncertainty Analysis';
  }
}

/**
 * Get description for a method
 */
function getMethodDescription(method: UncertaintyQuantificationMethod): string {
  switch (method) {
    case 'confidence_interval':
      return 'Confidence intervals represent the range of values within which the true value lies with a specified level of confidence. The shaded region shows the confidence interval around the measured values.';
    case 'monte_carlo':
      return 'Monte Carlo simulation generates multiple random samples based on the uncertainty distribution of measurements. Each point represents a possible outcome, providing a visual representation of the probability distribution.';
    case 'sensitivity':
      return 'Sensitivity analysis identifies how different parameters and factors contribute to the overall uncertainty in the experimental results. The bar chart shows the relative importance of each parameter.';
    case 'error_propagation':
      return 'Error propagation analyzes how measurement uncertainties propagate through calculations to affect the final result. The shaded region represents standard deviation bounds around the measured values.';
    default:
      return 'This analysis quantifies and visualizes the uncertainty in experimental measurements to help assess reliability and precision of results.';
  }
}
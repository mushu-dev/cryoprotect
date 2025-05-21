/**
 * Time Series Visualizer Component
 * 
 * A component for visualizing time series data from experiments with
 * uncertainty quantification, interactive controls, and annotations.
 */

import React, { useState, useCallback, useMemo } from 'react';
import { useQuery } from 'convex/react';
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
  ErrorBar,
  Brush
} from 'recharts';

// Types
export interface TimeSeriesData {
  timestamp: number;
  value: number;
  uncertainty?: number;
  metadata?: Record<string, any>;
}

export interface TimeSeriesVisualizerProps {
  timeSeriesId: Id<"timeSeries">;
  showUncertainty?: boolean;
  showRaw?: boolean;
  height?: number;
  width?: number;
  theme?: 'light' | 'dark';
  annotations?: Array<{
    label: string;
    time: number;
    description?: string;
    type?: 'point' | 'range';
    endTime?: number;
    color?: string;
  }>;
  compareWith?: Array<{
    id: Id<"timeSeries">;
    label: string;
    color?: string;
  }>;
  onPointClick?: (data: TimeSeriesData) => void;
  actions?: React.ReactNode;
}

/**
 * TimeSeriesVisualizer Component
 */
export function TimeSeriesVisualizer({
  timeSeriesId,
  showUncertainty = true,
  showRaw = false,
  height = 400,
  width = '100%',
  theme = 'light',
  annotations = [],
  compareWith = [],
  onPointClick,
  actions
}: TimeSeriesVisualizerProps) {
  // State for zoom
  const [zoomDomain, setZoomDomain] = useState<[number, number] | null>(null);
  const [isZooming, setIsZooming] = useState(false);
  const [zoomStart, setZoomStart] = useState<number | null>(null);
  
  // Fetch the time series metadata
  const timeSeries = useQuery(api.experiments.enhanced_index.getTimeSeries, { timeSeriesId });
  
  // Fetch the time series data
  const timeSeriesData = useQuery(api.experiments.enhanced_index.getTimeSeriesData, { timeSeriesId });
  
  // Fetch comparison data
  const comparisonData = useMemo(() => {
    return compareWith.map(comparison => {
      const data = useQuery(api.experiments.enhanced_index.getTimeSeriesData, { timeSeriesId: comparison.id });
      return {
        ...comparison,
        data: data || []
      };
    });
  }, [compareWith]);
  
  // Prepare the chart data
  const chartData = useMemo(() => {
    if (!timeSeriesData) return [];
    
    return timeSeriesData.map(dataPoint => ({
      timestamp: dataPoint.timestamp,
      time: new Date(dataPoint.timestamp).toLocaleTimeString(),
      date: new Date(dataPoint.timestamp).toLocaleDateString(),
      value: dataPoint.value,
      uncertainty: dataPoint.uncertainty,
      ...dataPoint.metadata
    }));
  }, [timeSeriesData]);
  
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
          uncertainty: 'rgba(96, 165, 250, 0.3)'
        }
      : {
          background: '#ffffff',
          text: '#111827',
          grid: '#e5e7eb',
          primary: '#3b82f6',
          uncertainty: 'rgba(59, 130, 246, 0.3)'
        };
  }, [theme]);
  
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
        
        {payload.map((entry: any, index: number) => (
          <p key={index} style={{ margin: '0', color: entry.color || colors.text }}>
            {entry.name}: {entry.value.toFixed(2)} {timeSeries?.units || ''}
            {entry.name === 'value' && data.uncertainty && (
              <span style={{ fontSize: '0.9em', opacity: 0.8 }}>
                {' '}±{data.uncertainty.toFixed(2)}
              </span>
            )}
          </p>
        ))}
        
        {/* Show any metadata */}
        {Object.entries(data)
          .filter(([key]) => !['timestamp', 'time', 'date', 'value', 'uncertainty'].includes(key))
          .map(([key, value]) => (
            <p key={key} style={{ margin: '0', fontSize: '0.9em', color: colors.text, opacity: 0.8 }}>
              {key}: {String(value)}
            </p>
          ))
        }
      </div>
    );
  }, [colors, timeSeries]);
  
  // If data is loading
  if (!timeSeries || !timeSeriesData) {
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
        Loading time series data...
      </div>
    );
  }
  
  // If no data
  if (timeSeriesData.length === 0) {
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
        <p>No data available for this time series.</p>
        {actions && <div>{actions}</div>}
      </div>
    );
  }
  
  return (
    <div style={{ 
      height, 
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
        marginBottom: '10px'
      }}>
        <div>
          <h3 style={{ margin: '0 0 5px', fontSize: '1.2em' }}>
            {timeSeries.name}
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
      
      <ResponsiveContainer width="100%" height={height - 80}>
        <LineChart
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
          
          {/* Main data line */}
          <Line
            name={timeSeries.name}
            type="monotone"
            dataKey="value"
            stroke={colors.primary}
            strokeWidth={2}
            dot={{ r: 3, strokeWidth: 1 }}
            activeDot={{ r: 6, strokeWidth: 1 }}
            isAnimationActive={true}
          >
            {/* Show error bars for uncertainty if enabled */}
            {showUncertainty && (
              <ErrorBar
                dataKey="uncertainty"
                width={4}
                strokeWidth={1}
                stroke={colors.primary}
                direction="y"
              />
            )}
          </Line>
          
          {/* Comparison data */}
          {comparisonData.map((comparison, index) => (
            comparison.data.length > 0 && (
              <Line
                key={comparison.id.toString()}
                name={comparison.label}
                type="monotone"
                data={comparison.data.map(point => ({
                  timestamp: point.timestamp,
                  value: point.value
                }))}
                dataKey="value"
                stroke={comparison.color || getComparisonColor(index, theme)}
                strokeWidth={1.5}
                dot={{ r: 2 }}
                activeDot={{ r: 5 }}
              />
            )
          ))}
          
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
        </LineChart>
      </ResponsiveContainer>
      
      {/* Raw data display if enabled */}
      {showRaw && (
        <div style={{
          marginTop: '20px',
          maxHeight: '200px',
          overflowY: 'auto',
          fontSize: '0.8em',
          backgroundColor: theme === 'dark' ? '#374151' : '#f9fafb',
          padding: '10px',
          borderRadius: '4px'
        }}>
          <table style={{ width: '100%', borderCollapse: 'collapse' }}>
            <thead>
              <tr>
                <th style={{ textAlign: 'left', padding: '4px 8px' }}>Time</th>
                <th style={{ textAlign: 'right', padding: '4px 8px' }}>Value</th>
                {showUncertainty && (
                  <th style={{ textAlign: 'right', padding: '4px 8px' }}>Uncertainty</th>
                )}
              </tr>
            </thead>
            <tbody>
              {chartData.map((data, index) => (
                <tr key={index} style={{ 
                  borderBottom: `1px solid ${colors.grid}`,
                  backgroundColor: index % 2 === 0 
                    ? (theme === 'dark' ? '#4b5563' : '#f3f4f6') 
                    : 'transparent'
                }}>
                  <td style={{ padding: '4px 8px' }}>
                    {new Date(data.timestamp).toLocaleString()}
                  </td>
                  <td style={{ textAlign: 'right', padding: '4px 8px' }}>
                    {data.value.toFixed(3)} {timeSeries.units || ''}
                  </td>
                  {showUncertainty && (
                    <td style={{ textAlign: 'right', padding: '4px 8px' }}>
                      {data.uncertainty ? `±${data.uncertainty.toFixed(3)}` : '—'}
                    </td>
                  )}
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      )}
    </div>
  );
}

/**
 * Get a color for comparison data based on index
 */
function getComparisonColor(index: number, theme: 'light' | 'dark'): string {
  const lightColors = [
    '#ef4444', // Red
    '#14b8a6', // Teal
    '#a855f7', // Purple
    '#f59e0b', // Amber
    '#84cc16', // Lime
    '#06b6d4', // Cyan
    '#ec4899', // Pink
    '#8b5cf6', // Violet
  ];
  
  const darkColors = [
    '#f87171', // Red
    '#2dd4bf', // Teal
    '#c084fc', // Purple
    '#fbbf24', // Amber
    '#a3e635', // Lime
    '#22d3ee', // Cyan
    '#f472b6', // Pink
    '#a78bfa', // Violet
  ];
  
  const colors = theme === 'dark' ? darkColors : lightColors;
  return colors[index % colors.length];
}
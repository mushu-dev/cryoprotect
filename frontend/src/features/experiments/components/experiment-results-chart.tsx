import React from 'react';

interface ExperimentResultsChartProps {
  experiment: any;
  type?: 'bar' | 'line' | 'radar';
}

/**
 * Component for visualizing experiment results
 * This is a placeholder implementation that would be replaced with a real chart library
 * like Chart.js, Recharts, or D3 in a production application
 */
export function ExperimentResultsChart({ experiment, type = 'bar' }: ExperimentResultsChartProps) {
  if (!experiment || !experiment.results) {
    return (
      <div className="p-6 text-center">
        <p className="text-muted-foreground">No results data available</p>
      </div>
    );
  }
  
  // Extract result metrics
  const metrics = {
    viability: parseInt(experiment.results.viability, 10) || 0,
    recovery: parseInt(experiment.results.recovery, 10) || 0,
    functionality: parseInt(experiment.results.functionality, 10) || 0
  };
  
  // Bar chart implementation
  const renderBarChart = () => {
    const maxValue = 100; // Assuming percentage values
    
    return (
      <div className="w-full py-4">
        <div className="space-y-6">
          {Object.entries(metrics).map(([key, value]) => (
            <div key={key} className="space-y-1">
              <div className="flex justify-between text-sm">
                <span className="capitalize">{key}</span>
                <span className="font-medium">{value}%</span>
              </div>
              <div className="w-full bg-muted rounded-full h-2.5">
                <div 
                  className="bg-primary h-2.5 rounded-full" 
                  style={{ width: `${(value / maxValue) * 100}%` }}
                ></div>
              </div>
            </div>
          ))}
        </div>
      </div>
    );
  };
  
  // Radar chart (simplified visual representation)
  const renderRadarChart = () => {
    return (
      <div className="relative w-full aspect-square max-w-[300px] mx-auto p-4">
        <div className="absolute inset-0 flex items-center justify-center">
          <div className="w-2/3 h-2/3 rounded-full border border-dashed border-muted-foreground/30"></div>
          <div className="absolute w-1/3 h-1/3 rounded-full border border-dashed border-muted-foreground/30"></div>
        </div>
        
        {/* Axis lines */}
        <div className="absolute inset-0 flex items-center justify-center">
          <div className="w-full h-[1px] bg-muted-foreground/30"></div>
        </div>
        <div className="absolute inset-0 flex items-center justify-center">
          <div className="h-full w-[1px] bg-muted-foreground/30"></div>
        </div>
        <div className="absolute inset-0 flex items-center justify-center">
          <div className="w-[70%] h-[1px] bg-muted-foreground/30 transform rotate-45"></div>
        </div>
        <div className="absolute inset-0 flex items-center justify-center">
          <div className="w-[70%] h-[1px] bg-muted-foreground/30 transform -rotate-45"></div>
        </div>
        
        {/* Data points */}
        <div className="absolute inset-0">
          <div 
            className="absolute h-2 w-2 rounded-full bg-primary" 
            style={{ 
              top: `${50 - ((metrics.viability / 100) * 40)}%`, 
              left: '50%',
              transform: 'translate(-50%, -50%)' 
            }}
          ></div>
          <div 
            className="absolute h-2 w-2 rounded-full bg-primary" 
            style={{ 
              bottom: `${50 - ((metrics.recovery / 100) * 40)}%`, 
              left: '50%',
              transform: 'translate(-50%, 50%)' 
            }}
          ></div>
          <div 
            className="absolute h-2 w-2 rounded-full bg-primary" 
            style={{ 
              top: '50%', 
              right: `${50 - ((metrics.functionality / 100) * 40)}%`,
              transform: 'translate(50%, -50%)' 
            }}
          ></div>
          
          {/* Connect the dots to form a polygon */}
          <svg className="absolute inset-0 w-full h-full" style={{ zIndex: -1 }}>
            <polygon 
              points={`
                ${50}% ${(50 - ((metrics.viability / 100) * 40))}%, 
                ${(50 + ((metrics.functionality / 100) * 40))}% 50%, 
                ${50}% ${(50 + ((metrics.recovery / 100) * 40))}%
              `}
              fill="rgba(var(--primary) / 0.2)"
              stroke="rgb(var(--primary))"
              strokeWidth="1"
            />
          </svg>
          
          {/* Labels */}
          <div className="absolute text-xs font-medium" style={{ top: '5%', left: '50%', transform: 'translateX(-50%)' }}>
            Viability
          </div>
          <div className="absolute text-xs font-medium" style={{ bottom: '5%', left: '50%', transform: 'translateX(-50%)' }}>
            Recovery
          </div>
          <div className="absolute text-xs font-medium" style={{ top: '50%', right: '5%', transform: 'translateY(-50%)' }}>
            Functionality
          </div>
        </div>
      </div>
    );
  };
  
  // Line chart implementation
  const renderLineChart = () => {
    const points = [
      { name: 'Viability', value: metrics.viability },
      { name: 'Recovery', value: metrics.recovery },
      { name: 'Functionality', value: metrics.functionality }
    ];
    
    const maxValue = 100;
    const chartHeight = 200;
    
    // Calculate point coordinates
    const pointsWithCoords = points.map((point, index) => {
      const x = (index / (points.length - 1)) * 100;
      const y = chartHeight - ((point.value / maxValue) * chartHeight);
      return { ...point, x, y };
    });
    
    // Create the path for the line
    let path = `M ${pointsWithCoords[0].x} ${pointsWithCoords[0].y}`;
    for (let i = 1; i < pointsWithCoords.length; i++) {
      path += ` L ${pointsWithCoords[i].x} ${pointsWithCoords[i].y}`;
    }
    
    return (
      <div className="w-full h-[250px] py-4">
        <svg width="100%" height={chartHeight} className="overflow-visible">
          {/* Grid lines */}
          <line x1="0" y1={chartHeight} x2="100%" y2={chartHeight} stroke="#e5e7eb" />
          <line x1="0" y1={chartHeight / 2} x2="100%" y2={chartHeight / 2} stroke="#e5e7eb" strokeDasharray="4" />
          <line x1="0" y1="0" x2="100%" y2="0" stroke="#e5e7eb" />
          
          {/* Vertical grid lines */}
          {points.map((point, index) => {
            const x = (index / (points.length - 1)) * 100 + '%';
            return (
              <line 
                key={`grid-${index}`}
                x1={x} 
                y1="0" 
                x2={x} 
                y2={chartHeight} 
                stroke="#e5e7eb" 
                strokeDasharray="4"
              />
            );
          })}
          
          {/* Line */}
          <path 
            d={path} 
            fill="none" 
            stroke="rgb(var(--primary))" 
            strokeWidth="2" 
          />
          
          {/* Area under line */}
          <path 
            d={`${path} L ${pointsWithCoords[pointsWithCoords.length-1].x} ${chartHeight} L ${pointsWithCoords[0].x} ${chartHeight} Z`} 
            fill="rgba(var(--primary) / 0.1)" 
          />
          
          {/* Points */}
          {pointsWithCoords.map((point, index) => (
            <circle 
              key={index}
              cx={`${point.x}%`} 
              cy={point.y} 
              r="4" 
              fill="white" 
              stroke="rgb(var(--primary))" 
              strokeWidth="2" 
            />
          ))}
          
          {/* X-axis labels */}
          {points.map((point, index) => (
            <text 
              key={`label-${index}`}
              x={`${(index / (points.length - 1)) * 100}%`} 
              y={chartHeight + 20} 
              textAnchor="middle" 
              className="text-xs font-medium text-muted-foreground"
            >
              {point.name}
            </text>
          ))}
          
          {/* Y-axis labels */}
          <text x="-5" y={chartHeight} textAnchor="end" className="text-xs font-medium text-muted-foreground">0%</text>
          <text x="-5" y={chartHeight / 2} textAnchor="end" className="text-xs font-medium text-muted-foreground">50%</text>
          <text x="-5" y="5" textAnchor="end" className="text-xs font-medium text-muted-foreground">100%</text>
        </svg>
      </div>
    );
  };
  
  return (
    <div className="bg-card rounded-lg border p-6">
      <div className="flex justify-between items-center mb-4">
        <h2 className="text-xl font-semibold">Experiment Results</h2>
        <div className="flex space-x-2">
          <button 
            onClick={() => {}} 
            className={`px-2 py-1 text-xs font-medium rounded-md ${
              type === 'bar' ? 'bg-primary text-primary-foreground' : 'bg-muted text-muted-foreground'
            }`}
          >
            Bar
          </button>
          <button 
            onClick={() => {}} 
            className={`px-2 py-1 text-xs font-medium rounded-md ${
              type === 'line' ? 'bg-primary text-primary-foreground' : 'bg-muted text-muted-foreground'
            }`}
          >
            Line
          </button>
          <button 
            onClick={() => {}} 
            className={`px-2 py-1 text-xs font-medium rounded-md ${
              type === 'radar' ? 'bg-primary text-primary-foreground' : 'bg-muted text-muted-foreground'
            }`}
          >
            Radar
          </button>
        </div>
      </div>
      
      {/* Chart based on type */}
      {type === 'bar' && renderBarChart()}
      {type === 'line' && renderLineChart()}
      {type === 'radar' && renderRadarChart()}
      
      {/* Summary stats */}
      <div className="mt-6 pt-6 border-t">
        <h3 className="text-sm font-medium mb-2">Summary</h3>
        <p className="text-sm text-muted-foreground">
          {experiment.results.notes || 'No additional notes provided for these results.'}
        </p>
      </div>
    </div>
  );
}
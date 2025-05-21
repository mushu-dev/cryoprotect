/**
 * MixtureCompositionChart Component
 * Chart visualization for mixture compositions
 */
import React, { useEffect, useRef, useState } from 'react';
import { Card } from '../ui/card';
import { randomColor } from '../../lib/utils';

export function MixtureCompositionChart({
  components,
  height = 300,
  title = 'Mixture Composition',
  showLegend = true,
}) {
  const chartRef = useRef(null);
  const [isClient, setIsClient] = useState(false);
  const [isLoading, setIsLoading] = useState(true);
  const [error, setError] = useState(null);
  
  // Generate consistent colors for components if not provided
  const componentsWithColors = components.map(component => ({
    ...component,
    color: component.color || randomColor(),
    concentration_unit: component.concentration_unit || '%'
  }));

  useEffect(() => {
    setIsClient(true);
  }, []);

  useEffect(() => {
    if (!isClient || !chartRef.current) return;
    
    const renderPlaceholder = (errorMessage) => {
      if (!chartRef.current) return;
      
      const container = chartRef.current;
      container.innerHTML = '';
      
      // Create placeholder
      const placeholder = document.createElement('div');
      placeholder.className = 'flex flex-col items-center justify-center w-full h-full bg-muted/30 rounded-md';
      placeholder.style.height = `${height}px`;
      
      // Add title
      const titleEl = document.createElement('h4');
      titleEl.className = 'font-medium mb-4';
      titleEl.textContent = title;
      placeholder.appendChild(titleEl);
      
      if (errorMessage) {
        const errorEl = document.createElement('p');
        errorEl.className = 'text-destructive text-sm mb-4';
        errorEl.textContent = errorMessage;
        placeholder.appendChild(errorEl);
      }
      
      // Create a simple visualization of the mixture components
      const componentsContainer = document.createElement('div');
      componentsContainer.className = 'flex flex-wrap justify-center gap-2 max-w-xs';
      
      componentsWithColors.forEach(component => {
        const componentEl = document.createElement('div');
        componentEl.className = 'flex items-center gap-2 bg-background p-2 rounded-md shadow-sm';
        
        const colorBox = document.createElement('div');
        colorBox.className = 'w-4 h-4 rounded-sm';
        colorBox.style.backgroundColor = component.color || randomColor();
        
        const nameEl = document.createElement('span');
        nameEl.className = 'text-xs font-medium';
        nameEl.textContent = component.name;
        
        const valueEl = document.createElement('span');
        valueEl.className = 'text-xs text-muted-foreground';
        valueEl.textContent = `${component.concentration}${component.concentration_unit}`;
        
        componentEl.appendChild(colorBox);
        componentEl.appendChild(nameEl);
        componentEl.appendChild(valueEl);
        componentsContainer.appendChild(componentEl);
      });
      
      placeholder.appendChild(componentsContainer);
      
      // Add loading message or error
      const statusEl = document.createElement('p');
      statusEl.className = 'text-xs text-muted-foreground mt-4';
      statusEl.textContent = errorMessage ? 'Error loading chart' : 'Loading chart...';
      placeholder.appendChild(statusEl);
      
      container.appendChild(placeholder);
    };

    const renderChart = async () => {
      setIsLoading(true);
      
      try {
        // Show placeholder while loading
        renderPlaceholder();
        
        // Import Plotly dynamically to avoid SSR issues
        const Plotly = await import('plotly.js-dist');
        
        // Check if components are provided
        if (!components.length) {
          setError('No mixture components provided');
          renderPlaceholder('No mixture components provided');
          return;
        }
        
        // Prepare data for pie chart
        const data = [{
          values: componentsWithColors.map(c => c.concentration),
          labels: componentsWithColors.map(c => c.name),
          type: 'pie',
          textinfo: 'label+percent',
          textposition: 'inside',
          automargin: true,
          marker: {
            colors: componentsWithColors.map(c => c.color)
          },
          hoverinfo: 'label+percent+value',
          hovertemplate: '%{label}: %{value}%{text}<extra></extra>',
          text: componentsWithColors.map(c => c.concentration_unit)
        }];
        
        // Chart layout
        const layout = {
          height: height,
          margin: { t: title ? 50 : 30, b: 30, l: 30, r: 30 },
          showlegend: showLegend,
          legend: showLegend ? { orientation: 'h', y: -0.2 } : undefined,
          title: title ? {
            text: title,
            font: {
              size: 16
            }
          } : undefined,
          paper_bgcolor: 'transparent',
          plot_bgcolor: 'transparent'
        };
        
        const config = {
          responsive: true,
          displayModeBar: false
        };
        
        // Create the chart
        Plotly.default.newPlot(chartRef.current, data, layout, config);
        setIsLoading(false);
        setError(null);
        
      } catch (err) {
        console.error('Error rendering mixture chart:', err);
        setError(err instanceof Error ? err.message : 'Failed to render chart');
        renderPlaceholder(err instanceof Error ? err.message : 'Failed to render chart');
        setIsLoading(false);
      }
    };
    
    renderChart();
    
    // Cleanup function
    return () => {
      if (chartRef.current) {
        if (typeof window !== 'undefined') {
          try {
            const Plotly = window.Plotly;
            if (Plotly && typeof Plotly.purge === 'function') {
              Plotly.purge(chartRef.current);
            }
          } catch (e) {
            // Silent catch if Plotly is not defined yet
          }
        }
        chartRef.current.innerHTML = '';
      }
    };
  }, [isClient, components, componentsWithColors, height, title, showLegend]);

  if (!isClient) {
    return (
      <div 
        className="flex flex-col items-center justify-center w-full bg-muted/30 rounded-md"
        style={{ height: `${height}px` }}
      >
        <p className="text-muted-foreground text-sm">Loading mixture composition...</p>
      </div>
    );
  }

  return (
    <Card className="overflow-hidden">
      <div 
        ref={chartRef} 
        className="chart-container w-full"
        style={{ height: `${height}px` }}
        role="img"
        aria-label="Pie chart showing mixture composition"
      />
    </Card>
  );
}
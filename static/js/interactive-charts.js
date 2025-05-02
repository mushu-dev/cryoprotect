/**
 * CryoProtect Analyzer - Interactive Charts Module
 * 
 * This module provides functions for creating interactive charts using Plotly.js
 */

const InteractiveCharts = (function() {
    // Chart color palette
    const COLORS = [
        '#2c3e50', // Primary
        '#3498db', // Secondary
        '#e74c3c', // Accent
        '#2ecc71', // Success
        '#f39c12', // Warning
        '#9b59b6', // Purple
        '#1abc9c', // Teal
        '#34495e', // Dark
        '#e67e22', // Orange
        '#95a5a6'  // Gray
    ];
    
    /**
     * Create an interactive bar chart
     * 
     * @param {string} containerId - ID of the container element
     * @param {Array} data - Array of data objects
     * @param {string} xField - Field name for x-axis
     * @param {string} yField - Field name for y-axis
     * @param {Object} options - Chart options
     */
    function createBarChart(containerId, data, xField, yField, options = {}) {
        const container = document.getElementById(containerId);
        if (!container) {
            console.error(`Container element with ID ${containerId} not found`);
            return;
        }
        
        // Default options
        const defaultOptions = {
            title: 'Bar Chart',
            xAxisTitle: xField,
            yAxisTitle: yField,
            color: COLORS[0],
            hoverTemplate: `%{x}: %{y}<extra></extra>`,
            showLegend: false,
            height: 400,
            margin: { t: 50, r: 50, b: 100, l: 80 }
        };
        
        // Merge options
        const chartOptions = { ...defaultOptions, ...options };
        
        // Prepare data for chart
        const xValues = data.map(item => item[xField]);
        const yValues = data.map(item => item[yField]);
        
        // Create trace
        const trace = {
            x: xValues,
            y: yValues,
            type: 'bar',
            marker: {
                color: chartOptions.color
            },
            hovertemplate: chartOptions.hoverTemplate
        };
        
        // Create layout
        const layout = {
            title: chartOptions.title,
            xaxis: {
                title: chartOptions.xAxisTitle,
                tickangle: -45
            },
            yaxis: {
                title: chartOptions.yAxisTitle
            },
            showlegend: chartOptions.showLegend,
            height: chartOptions.height,
            margin: chartOptions.margin,
            hovermode: 'closest'
        };
        
        // Create config
        const config = {
            responsive: true,
            displayModeBar: true,
            modeBarButtonsToRemove: ['lasso2d', 'select2d'],
            displaylogo: false,
            toImageButtonOptions: {
                format: 'png',
                filename: 'chart_export',
                height: 800,
                width: 1200,
                scale: 2
            }
        };
        
        // Create chart
        Plotly.newPlot(containerId, [trace], layout, config);
        
        // Add click event handler if provided
        if (options.onClick) {
            container.on('plotly_click', (data) => {
                const point = data.points[0];
                const pointData = {
                    x: point.x,
                    y: point.y,
                    index: point.pointIndex,
                    data: options.data[point.pointIndex]
                };
                options.onClick(pointData);
            });
        }
    }
    
    /**
     * Create an interactive line chart
     * 
     * @param {string} containerId - ID of the container element
     * @param {Array} data - Array of data objects
     * @param {string} xField - Field name for x-axis
     * @param {string} yField - Field name for y-axis
     * @param {Object} options - Chart options
     */
    function createLineChart(containerId, data, xField, yField, options = {}) {
        const container = document.getElementById(containerId);
        if (!container) {
            console.error(`Container element with ID ${containerId} not found`);
            return;
        }
        
        // Default options
        const defaultOptions = {
            title: 'Line Chart',
            xAxisTitle: xField,
            yAxisTitle: yField,
            color: COLORS[1],
            hoverTemplate: `%{x}: %{y}<extra></extra>`,
            showLegend: false,
            height: 400,
            margin: { t: 50, r: 50, b: 80, l: 80 },
            mode: 'lines+markers'
        };
        
        // Merge options
        const chartOptions = { ...defaultOptions, ...options };
        
        // Prepare data for chart
        const xValues = data.map(item => item[xField]);
        const yValues = data.map(item => item[yField]);
        
        // Create trace
        const trace = {
            x: xValues,
            y: yValues,
            type: 'scatter',
            mode: chartOptions.mode,
            marker: {
                color: chartOptions.color,
                size: 8
            },
            line: {
                color: chartOptions.color,
                width: 2
            },
            hovertemplate: chartOptions.hoverTemplate
        };
        
        // Create layout
        const layout = {
            title: chartOptions.title,
            xaxis: {
                title: chartOptions.xAxisTitle
            },
            yaxis: {
                title: chartOptions.yAxisTitle
            },
            showlegend: chartOptions.showLegend,
            height: chartOptions.height,
            margin: chartOptions.margin,
            hovermode: 'closest'
        };
        
        // Create config
        const config = {
            responsive: true,
            displayModeBar: true,
            modeBarButtonsToRemove: ['lasso2d', 'select2d'],
            displaylogo: false,
            toImageButtonOptions: {
                format: 'png',
                filename: 'chart_export',
                height: 800,
                width: 1200,
                scale: 2
            }
        };
        
        // Create chart
        Plotly.newPlot(containerId, [trace], layout, config);
        
        // Add click event handler if provided
        if (options.onClick) {
            container.on('plotly_click', (data) => {
                const point = data.points[0];
                const pointData = {
                    x: point.x,
                    y: point.y,
                    index: point.pointIndex,
                    data: options.data[point.pointIndex]
                };
                options.onClick(pointData);
            });
        }
    }
    
    /**
     * Create an interactive scatter plot
     * 
     * @param {string} containerId - ID of the container element
     * @param {Array} data - Array of data objects
     * @param {string} xField - Field name for x-axis
     * @param {string} yField - Field name for y-axis
     * @param {string} labelField - Field name for point labels
     * @param {Object} options - Chart options
     */
    function createScatterPlot(containerId, data, xField, yField, labelField, options = {}) {
        const container = document.getElementById(containerId);
        if (!container) {
            console.error(`Container element with ID ${containerId} not found`);
            return;
        }
        
        // Default options
        const defaultOptions = {
            title: 'Scatter Plot',
            xAxisTitle: xField,
            yAxisTitle: yField,
            color: COLORS[2],
            hoverTemplate: `%{text}<br>${xField}: %{x}<br>${yField}: %{y}<extra></extra>`,
            showLegend: false,
            height: 500,
            margin: { t: 50, r: 50, b: 80, l: 80 },
            colorField: null,
            sizeField: null,
            minSize: 6,
            maxSize: 20
        };
        
        // Merge options
        const chartOptions = { ...defaultOptions, ...options };
        
        // Prepare data for chart
        const xValues = data.map(item => item[xField]);
        const yValues = data.map(item => item[yField]);
        const labels = data.map(item => item[labelField] || '');
        
        // Create trace
        const trace = {
            x: xValues,
            y: yValues,
            text: labels,
            type: 'scatter',
            mode: 'markers',
            marker: {
                color: chartOptions.color,
                size: 10
            },
            hovertemplate: chartOptions.hoverTemplate
        };
        
        // Apply color mapping if colorField is provided
        if (chartOptions.colorField && data[0] && data[0][chartOptions.colorField] !== undefined) {
            const colorValues = data.map(item => item[chartOptions.colorField]);
            trace.marker.color = colorValues;
            trace.marker.colorscale = 'Viridis';
            trace.marker.showscale = true;
            trace.marker.colorbar = {
                title: chartOptions.colorField,
                thickness: 20
            };
        }
        
        // Apply size mapping if sizeField is provided
        if (chartOptions.sizeField && data[0] && data[0][chartOptions.sizeField] !== undefined) {
            const sizeValues = data.map(item => item[chartOptions.sizeField]);
            const minValue = Math.min(...sizeValues);
            const maxValue = Math.max(...sizeValues);
            const range = maxValue - minValue;
            
            trace.marker.size = sizeValues.map(value => {
                if (range === 0) return chartOptions.minSize;
                const normalized = (value - minValue) / range;
                return chartOptions.minSize + normalized * (chartOptions.maxSize - chartOptions.minSize);
            });
        }
        
        // Create layout
        const layout = {
            title: chartOptions.title,
            xaxis: {
                title: chartOptions.xAxisTitle
            },
            yaxis: {
                title: chartOptions.yAxisTitle
            },
            showlegend: chartOptions.showLegend,
            height: chartOptions.height,
            margin: chartOptions.margin,
            hovermode: 'closest'
        };
        
        // Create config
        const config = {
            responsive: true,
            displayModeBar: true,
            modeBarButtonsToRemove: [],
            displaylogo: false,
            toImageButtonOptions: {
                format: 'png',
                filename: 'chart_export',
                height: 800,
                width: 1200,
                scale: 2
            }
        };
        
        // Create chart
        Plotly.newPlot(containerId, [trace], layout, config);
        
        // Add click event handler if provided
        if (options.onClick) {
            container.on('plotly_click', (data) => {
                const point = data.points[0];
                const pointData = {
                    x: point.x,
                    y: point.y,
                    text: point.text,
                    index: point.pointIndex,
                    data: options.data[point.pointIndex]
                };
                options.onClick(pointData);
            });
        }
    }
    
    /**
     * Create an interactive pie chart
     * 
     * @param {string} containerId - ID of the container element
     * @param {Array} labels - Array of labels
     * @param {Array} values - Array of values
     * @param {Object} options - Chart options
     */
    function createPieChart(containerId, labels, values, options = {}) {
        const container = document.getElementById(containerId);
        if (!container) {
            console.error(`Container element with ID ${containerId} not found`);
            return;
        }
        
        // Default options
        const defaultOptions = {
            title: 'Pie Chart',
            colors: COLORS,
            hoverTemplate: '%{label}: %{value} (%{percent})<extra></extra>',
            showLegend: true,
            height: 400,
            margin: { t: 50, r: 50, b: 50, l: 50 },
            hole: 0
        };
        
        // Merge options
        const chartOptions = { ...defaultOptions, ...options };
        
        // Create trace
        const trace = {
            labels: labels,
            values: values,
            type: 'pie',
            marker: {
                colors: chartOptions.colors
            },
            hovertemplate: chartOptions.hoverTemplate,
            hole: chartOptions.hole,
            textinfo: 'label+percent',
            insidetextorientation: 'radial'
        };
        
        // Create layout
        const layout = {
            title: chartOptions.title,
            showlegend: chartOptions.showLegend,
            height: chartOptions.height,
            margin: chartOptions.margin,
            legend: {
                orientation: 'h',
                y: -0.2
            }
        };
        
        // Create config
        const config = {
            responsive: true,
            displayModeBar: true,
            modeBarButtonsToRemove: ['lasso2d', 'select2d'],
            displaylogo: false,
            toImageButtonOptions: {
                format: 'png',
                filename: 'chart_export',
                height: 800,
                width: 1200,
                scale: 2
            }
        };
        
        // Create chart
        Plotly.newPlot(containerId, [trace], layout, config);
        
        // Add click event handler if provided
        if (options.onClick) {
            container.on('plotly_click', (data) => {
                const point = data.points[0];
                const pointData = {
                    label: point.label,
                    value: point.value,
                    percent: point.percent,
                    index: point.pointIndex
                };
                options.onClick(pointData);
            });
        }
    }
    
    /**
     * Create an interactive radar chart
     * 
     * @param {string} containerId - ID of the container element
     * @param {Array} categories - Array of category names
     * @param {Array} values - Array of values
     * @param {Object} options - Chart options
     */
    function createRadarChart(containerId, categories, values, options = {}) {
        const container = document.getElementById(containerId);
        if (!container) {
            console.error(`Container element with ID ${containerId} not found`);
            return;
        }
        
        // Default options
        const defaultOptions = {
            title: 'Radar Chart',
            color: COLORS[0],
            fillColor: `${COLORS[0]}40`,
            hoverTemplate: '%{theta}: %{r}<extra></extra>',
            showLegend: false,
            height: 400,
            margin: { t: 50, r: 50, b: 50, l: 50 },
            maxValue: null
        };
        
        // Merge options
        const chartOptions = { ...defaultOptions, ...options };
        
        // Create trace
        const trace = {
            type: 'scatterpolar',
            r: values,
            theta: categories,
            fill: 'toself',
            fillcolor: chartOptions.fillColor,
            line: {
                color: chartOptions.color,
                width: 2
            },
            hovertemplate: chartOptions.hoverTemplate
        };
        
        // Create layout
        const layout = {
            title: chartOptions.title,
            showlegend: chartOptions.showLegend,
            height: chartOptions.height,
            margin: chartOptions.margin,
            polar: {
                radialaxis: {
                    visible: true,
                    range: [0, chartOptions.maxValue || Math.max(...values) * 1.1]
                }
            }
        };
        
        // Create config
        const config = {
            responsive: true,
            displayModeBar: true,
            modeBarButtonsToRemove: ['lasso2d', 'select2d'],
            displaylogo: false,
            toImageButtonOptions: {
                format: 'png',
                filename: 'chart_export',
                height: 800,
                width: 1200,
                scale: 2
            }
        };
        
        // Create chart
        Plotly.newPlot(containerId, [trace], layout, config);
        
        // Add click event handler if provided
        if (options.onClick) {
            container.on('plotly_click', (data) => {
                const point = data.points[0];
                const pointData = {
                    category: point.theta,
                    value: point.r,
                    index: point.pointIndex
                };
                options.onClick(pointData);
            });
        }
    }
    
    /**
     * Create an interactive heatmap
     * 
     * @param {string} containerId - ID of the container element
     * @param {Array} xLabels - Array of x-axis labels
     * @param {Array} yLabels - Array of y-axis labels
     * @param {Array} zValues - 2D array of z values
     * @param {Object} options - Chart options
     */
    function createHeatmap(containerId, xLabels, yLabels, zValues, options = {}) {
        const container = document.getElementById(containerId);
        if (!container) {
            console.error(`Container element with ID ${containerId} not found`);
            return;
        }
        
        // Default options
        const defaultOptions = {
            title: 'Heatmap',
            xAxisTitle: '',
            yAxisTitle: '',
            colorscale: 'Viridis',
            hoverTemplate: '%{x}, %{y}: %{z}<extra></extra>',
            showScale: true,
            height: 500,
            margin: { t: 50, r: 80, b: 80, l: 80 },
            annotate: false,
            annotationFormat: '.2f'
        };
        
        // Merge options
        const chartOptions = { ...defaultOptions, ...options };
        
        // Create trace
        const trace = {
            x: xLabels,
            y: yLabels,
            z: zValues,
            type: 'heatmap',
            colorscale: chartOptions.colorscale,
            showscale: chartOptions.showScale,
            hovertemplate: chartOptions.hoverTemplate
        };
        
        // Create layout
        const layout = {
            title: chartOptions.title,
            xaxis: {
                title: chartOptions.xAxisTitle
            },
            yaxis: {
                title: chartOptions.yAxisTitle
            },
            height: chartOptions.height,
            margin: chartOptions.margin
        };
        
        // Add annotations if requested
        if (chartOptions.annotate) {
            const annotations = [];
            
            for (let i = 0; i < yLabels.length; i++) {
                for (let j = 0; j < xLabels.length; j++) {
                    annotations.push({
                        x: xLabels[j],
                        y: yLabels[i],
                        text: zValues[i][j].toFixed(2),
                        showarrow: false,
                        font: {
                            color: getContrastColor(zValues[i][j], zValues)
                        }
                    });
                }
            }
            
            layout.annotations = annotations;
        }
        
        // Create config
        const config = {
            responsive: true,
            displayModeBar: true,
            modeBarButtonsToRemove: ['lasso2d', 'select2d'],
            displaylogo: false,
            toImageButtonOptions: {
                format: 'png',
                filename: 'chart_export',
                height: 800,
                width: 1200,
                scale: 2
            }
        };
        
        // Create chart
        Plotly.newPlot(containerId, [trace], layout, config);
        
        // Add click event handler if provided
        if (options.onClick) {
            container.on('plotly_click', (data) => {
                const point = data.points[0];
                const pointData = {
                    x: point.x,
                    y: point.y,
                    z: point.z,
                    xIndex: point.pointIndex[1],
                    yIndex: point.pointIndex[0]
                };
                options.onClick(pointData);
            });
        }
    }
    
    /**
     * Create an interactive 3D scatter plot
     * 
     * @param {string} containerId - ID of the container element
     * @param {Array} data - Array of data objects
     * @param {string} xField - Field name for x-axis
     * @param {string} yField - Field name for y-axis
     * @param {string} zField - Field name for z-axis
     * @param {string} labelField - Field name for point labels
     * @param {Object} options - Chart options
     */
    function create3DScatterPlot(containerId, data, xField, yField, zField, labelField, options = {}) {
        const container = document.getElementById(containerId);
        if (!container) {
            console.error(`Container element with ID ${containerId} not found`);
            return;
        }
        
        // Default options
        const defaultOptions = {
            title: '3D Scatter Plot',
            xAxisTitle: xField,
            yAxisTitle: yField,
            zAxisTitle: zField,
            color: COLORS[3],
            hoverTemplate: `%{text}<br>${xField}: %{x}<br>${yField}: %{y}<br>${zField}: %{z}<extra></extra>`,
            showLegend: false,
            height: 600,
            margin: { t: 50, r: 50, b: 50, l: 50 },
            colorField: null,
            sizeField: null,
            minSize: 4,
            maxSize: 12
        };
        
        // Merge options
        const chartOptions = { ...defaultOptions, ...options };
        
        // Prepare data for chart
        const xValues = data.map(item => item[xField]);
        const yValues = data.map(item => item[yField]);
        const zValues = data.map(item => item[zField]);
        const labels = data.map(item => item[labelField] || '');
        
        // Create trace
        const trace = {
            x: xValues,
            y: yValues,
            z: zValues,
            text: labels,
            type: 'scatter3d',
            mode: 'markers',
            marker: {
                color: chartOptions.color,
                size: 6
            },
            hovertemplate: chartOptions.hoverTemplate
        };
        
        // Apply color mapping if colorField is provided
        if (chartOptions.colorField && data[0] && data[0][chartOptions.colorField] !== undefined) {
            const colorValues = data.map(item => item[chartOptions.colorField]);
            trace.marker.color = colorValues;
            trace.marker.colorscale = 'Viridis';
            trace.marker.showscale = true;
            trace.marker.colorbar = {
                title: chartOptions.colorField,
                thickness: 20
            };
        }
        
        // Apply size mapping if sizeField is provided
        if (chartOptions.sizeField && data[0] && data[0][chartOptions.sizeField] !== undefined) {
            const sizeValues = data.map(item => item[chartOptions.sizeField]);
            const minValue = Math.min(...sizeValues);
            const maxValue = Math.max(...sizeValues);
            const range = maxValue - minValue;
            
            trace.marker.size = sizeValues.map(value => {
                if (range === 0) return chartOptions.minSize;
                const normalized = (value - minValue) / range;
                return chartOptions.minSize + normalized * (chartOptions.maxSize - chartOptions.minSize);
            });
        }
        
        // Create layout
        const layout = {
            title: chartOptions.title,
            scene: {
                xaxis: {
                    title: chartOptions.xAxisTitle
                },
                yaxis: {
                    title: chartOptions.yAxisTitle
                },
                zaxis: {
                    title: chartOptions.zAxisTitle
                }
            },
            showlegend: chartOptions.showLegend,
            height: chartOptions.height,
            margin: chartOptions.margin
        };
        
        // Create config
        const config = {
            responsive: true,
            displayModeBar: true,
            modeBarButtonsToRemove: [],
            displaylogo: false,
            toImageButtonOptions: {
                format: 'png',
                filename: 'chart_export',
                height: 800,
                width: 1200,
                scale: 2
            }
        };
        
        // Create chart
        Plotly.newPlot(containerId, [trace], layout, config);
        
        // Add click event handler if provided
        if (options.onClick) {
            container.on('plotly_click', (data) => {
                const point = data.points[0];
                const pointData = {
                    x: point.x,
                    y: point.y,
                    z: point.z,
                    text: point.text,
                    index: point.pointIndex,
                    data: options.data[point.pointIndex]
                };
                options.onClick(pointData);
            });
        }
    }
    
    /**
     * Helper function to get contrast color for heatmap annotations
     * 
     * @param {number} value - The value at the current position
     * @param {Array} zValues - 2D array of all z values
     * @returns {string} - Color for the annotation text
     */
    function getContrastColor(value, zValues) {
        // Flatten the 2D array
        const flatValues = zValues.flat();
        
        // Get min and max values
        const minValue = Math.min(...flatValues);
        const maxValue = Math.max(...flatValues);
        
        // Calculate normalized value (0-1)
        const range = maxValue - minValue;
        const normalized = range === 0 ? 0.5 : (value - minValue) / range;
        
        // Return white for dark backgrounds, black for light backgrounds
        return normalized > 0.6 ? 'white' : 'black';
    }
    
    // Public API
    return {
        createBarChart,
        createLineChart,
        createScatterPlot,
        createPieChart,
        createRadarChart,
        createHeatmap,
        create3DScatterPlot
    };
})();

// Export the module
if (typeof module !== 'undefined' && module.exports) {
    module.exports = InteractiveCharts;
} else {
    window.InteractiveCharts = InteractiveCharts;
}
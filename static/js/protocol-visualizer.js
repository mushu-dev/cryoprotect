/**
 * CryoProtect Analyzer - Protocol Visualizer
 * 
 * This module provides functions for interactive visualization of cryopreservation protocols
 */

const ProtocolVisualizer = (function() {
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
     * Create an interactive timeline visualization for a protocol
     * 
     * @param {string} containerId - ID of the container element
     * @param {Object} protocol - Protocol object with steps
     * @param {Object} options - Visualization options
     */
    function createProtocolTimeline(containerId, protocol, options = {}) {
        const container = document.getElementById(containerId);
        if (!container) {
            console.error(`Container element with ID ${containerId} not found`);
            return;
        }
        
        if (!protocol || !protocol.steps || protocol.steps.length === 0) {
            container.innerHTML = '<div class="alert alert-warning">No protocol steps found</div>';
            return;
        }
        
        // Default options
        const defaultOptions = {
            title: `Protocol Timeline: ${protocol.name || 'Unnamed Protocol'}`,
            height: 400,
            showConcentration: true,
            showTemperature: true,
            onClick: null
        };
        
        // Merge options
        const chartOptions = { ...defaultOptions, ...options };
        
        // Prepare data for chart
        const steps = protocol.steps.sort((a, b) => a.step_number - b.step_number);
        
        // Calculate cumulative time for x-axis
        let cumulativeTime = 0;
        const timePoints = steps.map(step => {
            const startTime = cumulativeTime;
            cumulativeTime += step.duration;
            return {
                step: step,
                startTime: startTime,
                endTime: cumulativeTime
            };
        });
        
        // Create traces
        const traces = [];
        
        // Temperature trace
        if (chartOptions.showTemperature) {
            const temperatureTrace = {
                x: [],
                y: [],
                type: 'scatter',
                mode: 'lines+markers',
                name: 'Temperature (°C)',
                line: {
                    color: COLORS[2],
                    width: 3
                },
                marker: {
                    size: 8,
                    color: COLORS[2]
                },
                yaxis: 'y2',
                hovertemplate: 'Time: %{x} min<br>Temperature: %{y}°C<extra></extra>'
            };
            
            // Add points for each step (start and end)
            timePoints.forEach((point, index) => {
                const step = point.step;
                
                // Add start point
                temperatureTrace.x.push(point.startTime);
                temperatureTrace.y.push(step.temperature);
                
                // Add end point
                temperatureTrace.x.push(point.endTime);
                temperatureTrace.y.push(step.temperature);
            });
            
            traces.push(temperatureTrace);
        }
        
        // Concentration trace
        if (chartOptions.showConcentration) {
            const concentrationTrace = {
                x: [],
                y: [],
                type: 'scatter',
                mode: 'lines+markers',
                name: 'Concentration (%)',
                line: {
                    color: COLORS[1],
                    width: 3
                },
                marker: {
                    size: 8,
                    color: COLORS[1]
                },
                hovertemplate: 'Time: %{x} min<br>Concentration: %{y}%<extra></extra>'
            };
            
            // Add points for each step (start and end)
            timePoints.forEach((point, index) => {
                const step = point.step;
                
                // Add start point
                concentrationTrace.x.push(point.startTime);
                concentrationTrace.y.push(step.concentration);
                
                // Add end point
                concentrationTrace.x.push(point.endTime);
                concentrationTrace.y.push(step.concentration);
            });
            
            traces.push(concentrationTrace);
        }
        
        // Create step regions
        timePoints.forEach((point, index) => {
            const step = point.step;
            
            // Create a shape for the step
            const shape = {
                type: 'rect',
                xref: 'x',
                yref: 'paper',
                x0: point.startTime,
                x1: point.endTime,
                y0: 0,
                y1: 1,
                fillcolor: `${COLORS[index % COLORS.length]}20`,
                line: {
                    width: 0
                }
            };
            
            // Add step number annotation
            const annotation = {
                x: (point.startTime + point.endTime) / 2,
                y: 1.05,
                xref: 'x',
                yref: 'paper',
                text: `Step ${step.step_number}`,
                showarrow: false,
                font: {
                    size: 12,
                    color: COLORS[index % COLORS.length]
                }
            };
            
            // Add to layout
            if (!layout.shapes) layout.shapes = [];
            layout.shapes.push(shape);
            
            if (!layout.annotations) layout.annotations = [];
            layout.annotations.push(annotation);
        });
        
        // Create layout
        const layout = {
            title: chartOptions.title,
            height: chartOptions.height,
            xaxis: {
                title: 'Time (minutes)',
                range: [0, cumulativeTime]
            },
            yaxis: {
                title: 'Concentration (%)',
                range: [0, Math.max(...protocol.steps.map(step => step.concentration)) * 1.1]
            },
            yaxis2: {
                title: 'Temperature (°C)',
                overlaying: 'y',
                side: 'right',
                range: [
                    Math.min(...protocol.steps.map(step => step.temperature)) - 5,
                    Math.max(...protocol.steps.map(step => step.temperature)) + 5
                ]
            },
            hovermode: 'closest',
            legend: {
                orientation: 'h',
                y: -0.2
            },
            shapes: [],
            annotations: []
        };
        
        // Create config
        const config = {
            responsive: true,
            displayModeBar: true,
            modeBarButtonsToRemove: ['lasso2d', 'select2d'],
            displaylogo: false,
            toImageButtonOptions: {
                format: 'png',
                filename: 'protocol_timeline',
                height: 800,
                width: 1200,
                scale: 2
            }
        };
        
        // Create chart
        Plotly.newPlot(containerId, traces, layout, config);
        
        // Add click event handler if provided
        if (chartOptions.onClick) {
            container.on('plotly_click', (data) => {
                const point = data.points[0];
                const time = point.x;
                
                // Find which step this time belongs to
                const stepPoint = timePoints.find(p => time >= p.startTime && time <= p.endTime);
                if (stepPoint) {
                    chartOptions.onClick(stepPoint.step, stepPoint.step.step_number - 1);
                }
            });
        }
    }
    
    /**
     * Create an interactive heatmap visualization for protocol comparison
     * 
     * @param {string} containerId - ID of the container element
     * @param {Array} protocols - Array of protocol objects
     * @param {Object} options - Visualization options
     */
    function createProtocolComparisonHeatmap(containerId, protocols, options = {}) {
        const container = document.getElementById(containerId);
        if (!container) {
            console.error(`Container element with ID ${containerId} not found`);
            return;
        }
        
        if (!protocols || protocols.length < 2) {
            container.innerHTML = '<div class="alert alert-warning">At least two protocols are required for comparison</div>';
            return;
        }
        
        // Default options
        const defaultOptions = {
            title: 'Protocol Comparison',
            height: 500,
            parameter: 'concentration', // 'concentration' or 'temperature'
            colorScale: 'Viridis',
            annotate: true,
            onClick: null
        };
        
        // Merge options
        const chartOptions = { ...defaultOptions, ...options };
        
        // Get max step count
        const maxSteps = Math.max(...protocols.map(p => p.steps ? p.steps.length : 0));
        
        // Prepare data for heatmap
        const zValues = [];
        const yLabels = [];
        
        protocols.forEach(protocol => {
            if (!protocol.steps || protocol.steps.length === 0) return;
            
            // Sort steps by step number
            const sortedSteps = protocol.steps.sort((a, b) => a.step_number - b.step_number);
            
            // Create row for this protocol
            const row = Array(maxSteps).fill(null);
            
            // Fill in values
            sortedSteps.forEach((step, index) => {
                if (index < maxSteps) {
                    row[index] = chartOptions.parameter === 'temperature' ? 
                        step.temperature : step.concentration;
                }
            });
            
            // Add to matrix
            zValues.push(row);
            yLabels.push(protocol.name || `Protocol ${protocol.id}`);
        });
        
        // Create x-axis labels
        const xLabels = Array.from({length: maxSteps}, (_, i) => `Step ${i + 1}`);
        
        // Create trace
        const trace = {
            z: zValues,
            x: xLabels,
            y: yLabels,
            type: 'heatmap',
            colorscale: chartOptions.colorScale,
            hovertemplate: '%{y}, %{x}: %{z}<extra></extra>'
        };
        
        // Create layout
        const layout = {
            title: chartOptions.title,
            height: chartOptions.height,
            xaxis: {
                title: 'Protocol Steps'
            },
            yaxis: {
                title: 'Protocols'
            }
        };
        
        // Add annotations if requested
        if (chartOptions.annotate) {
            const annotations = [];
            
            for (let i = 0; i < yLabels.length; i++) {
                for (let j = 0; j < xLabels.length; j++) {
                    if (zValues[i][j] !== null) {
                        annotations.push({
                            x: xLabels[j],
                            y: yLabels[i],
                            text: zValues[i][j].toFixed(1),
                            showarrow: false,
                            font: {
                                color: getContrastColor(zValues[i][j], zValues)
                            }
                        });
                    }
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
                filename: 'protocol_comparison',
                height: 800,
                width: 1200,
                scale: 2
            }
        };
        
        // Create chart
        Plotly.newPlot(containerId, [trace], layout, config);
        
        // Add click event handler if provided
        if (chartOptions.onClick) {
            container.on('plotly_click', (data) => {
                const point = data.points[0];
                const protocolIndex = point.pointIndex[0];
                const stepIndex = point.pointIndex[1];
                const protocol = protocols[protocolIndex];
                
                if (protocol && protocol.steps && protocol.steps.length > stepIndex) {
                    const step = protocol.steps.find(s => s.step_number === stepIndex + 1);
                    if (step) {
                        chartOptions.onClick(protocol, step, protocolIndex, stepIndex);
                    }
                }
            });
        }
    }
    
    /**
     * Create an interactive parallel coordinates plot for protocol comparison
     * 
     * @param {string} containerId - ID of the container element
     * @param {Array} protocols - Array of protocol objects
     * @param {Object} options - Visualization options
     */
    function createProtocolParallelCoordinates(containerId, protocols, options = {}) {
        const container = document.getElementById(containerId);
        if (!container) {
            console.error(`Container element with ID ${containerId} not found`);
            return;
        }
        
        if (!protocols || protocols.length === 0) {
            container.innerHTML = '<div class="alert alert-warning">No protocols found</div>';
            return;
        }
        
        // Default options
        const defaultOptions = {
            title: 'Protocol Parameter Comparison',
            height: 500,
            colorBy: 'protocol', // 'protocol' or 'score'
            onClick: null
        };
        
        // Merge options
        const chartOptions = { ...defaultOptions, ...options };
        
        // Define dimensions for parallel coordinates
        const dimensions = [
            {
                label: 'Target Concentration (%)',
                values: protocols.map(p => p.target_concentration || 0)
            },
            {
                label: 'Starting Temperature (°C)',
                values: protocols.map(p => p.starting_temperature || 0)
            },
            {
                label: 'Target Temperature (°C)',
                values: protocols.map(p => p.target_temperature || 0)
            },
            {
                label: 'Step Count',
                values: protocols.map(p => p.steps ? p.steps.length : 0)
            },
            {
                label: 'Total Duration (min)',
                values: protocols.map(p => p.steps ? 
                    p.steps.reduce((total, step) => total + step.duration, 0) : 0
                )
            }
        ];
        
        // Add score dimension if available
        if (protocols.some(p => p.score !== undefined)) {
            dimensions.push({
                label: 'Protocol Score',
                values: protocols.map(p => p.score !== undefined ? p.score : 0)
            });
        }
        
        // Prepare data for parallel coordinates
        const trace = {
            type: 'parcoords',
            line: {
                color: chartOptions.colorBy === 'score' && protocols.some(p => p.score !== undefined) ?
                    protocols.map(p => p.score !== undefined ? p.score : 0) :
                    protocols.map((_, i) => i),
                colorscale: 'Viridis'
            },
            dimensions: dimensions.map(dim => ({
                label: dim.label,
                values: dim.values
            }))
        };
        
        // Create layout
        const layout = {
            title: chartOptions.title,
            height: chartOptions.height
        };
        
        // Create config
        const config = {
            responsive: true,
            displayModeBar: true,
            modeBarButtonsToRemove: ['lasso2d', 'select2d'],
            displaylogo: false,
            toImageButtonOptions: {
                format: 'png',
                filename: 'protocol_parallel_coordinates',
                height: 800,
                width: 1200,
                scale: 2
            }
        };
        
        // Create chart
        Plotly.newPlot(containerId, [trace], layout, config);
    }
    
    /**
     * Create an interactive 3D surface plot for cooling rate visualization
     * 
     * @param {string} containerId - ID of the container element
     * @param {Object} protocol - Protocol object with steps
     * @param {Object} options - Visualization options
     */
    function createCoolingRateSurface(containerId, protocol, options = {}) {
        const container = document.getElementById(containerId);
        if (!container) {
            console.error(`Container element with ID ${containerId} not found`);
            return;
        }
        
        if (!protocol || !protocol.steps || protocol.steps.length === 0) {
            container.innerHTML = '<div class="alert alert-warning">No protocol steps found</div>';
            return;
        }
        
        // Default options
        const defaultOptions = {
            title: `Cooling Rate Visualization: ${protocol.name || 'Unnamed Protocol'}`,
            height: 600,
            resolution: 20, // Number of points in each dimension
            onClick: null
        };
        
        // Merge options
        const chartOptions = { ...defaultOptions, ...options };
        
        // Sort steps by step number
        const steps = protocol.steps.sort((a, b) => a.step_number - b.step_number);
        
        // Calculate cumulative time for each step
        let cumulativeTime = 0;
        const timePoints = steps.map(step => {
            const startTime = cumulativeTime;
            cumulativeTime += step.duration;
            return {
                step: step,
                startTime: startTime,
                endTime: cumulativeTime
            };
        });
        
        // Generate grid for surface
        const resolution = chartOptions.resolution;
        const x = []; // Time
        const y = []; // Position (depth)
        const z = []; // Temperature
        
        // Create time points
        for (let i = 0; i < resolution; i++) {
            const time = (cumulativeTime * i) / (resolution - 1);
            x.push(time);
        }
        
        // Create position points (depth from surface)
        for (let i = 0; i < resolution; i++) {
            const position = i / (resolution - 1);
            y.push(position);
        }
        
        // Create temperature grid
        for (let i = 0; i < resolution; i++) {
            const row = [];
            const time = x[i];
            
            // Find which step this time belongs to
            const stepPoint = timePoints.find(p => time >= p.startTime && time <= p.endTime) || timePoints[timePoints.length - 1];
            const step = stepPoint.step;
            
            for (let j = 0; j < resolution; j++) {
                const position = y[j];
                
                // Calculate temperature at this position and time
                // This is a simplified model - in reality, heat transfer would be more complex
                // We're using a simple exponential decay from surface to center
                const surfaceTemp = step.temperature;
                const centerTemp = steps[0].temperature; // Starting temperature at center
                
                // Calculate temperature gradient based on position and time in step
                const stepProgress = (time - stepPoint.startTime) / (stepPoint.endTime - stepPoint.startTime);
                const depthFactor = Math.exp(-5 * position); // Exponential decay with depth
                const tempDiff = surfaceTemp - centerTemp;
                
                // Temperature at this point
                const temperature = centerTemp + tempDiff * depthFactor * Math.min(1, stepProgress * 2);
                
                row.push(temperature);
            }
            
            z.push(row);
        }
        
        // Create trace
        const trace = {
            type: 'surface',
            x: x,
            y: y,
            z: z,
            colorscale: 'Viridis',
            contours: {
                z: {
                    show: true,
                    usecolormap: true,
                    highlightcolor: "#42f462",
                    project: {z: true}
                }
            },
            hovertemplate: 'Time: %{x} min<br>Depth: %{y}<br>Temperature: %{z}°C<extra></extra>'
        };
        
        // Create layout
        const layout = {
            title: chartOptions.title,
            height: chartOptions.height,
            scene: {
                xaxis: {
                    title: 'Time (minutes)'
                },
                yaxis: {
                    title: 'Depth (0=surface, 1=center)'
                },
                zaxis: {
                    title: 'Temperature (°C)'
                },
                camera: {
                    eye: {x: 1.5, y: 1.5, z: 1}
                }
            }
        };
        
        // Create config
        const config = {
            responsive: true,
            displayModeBar: true,
            modeBarButtonsToRemove: [],
            displaylogo: false,
            toImageButtonOptions: {
                format: 'png',
                filename: 'cooling_rate_surface',
                height: 800,
                width: 1200,
                scale: 2
            }
        };
        
        // Create chart
        Plotly.newPlot(containerId, [trace], layout, config);
    }
    
    /**
     * Helper function to get contrast color for heatmap annotations
     * 
     * @param {number} value - The value at the current position
     * @param {Array} zValues - 2D array of all z values
     * @returns {string} - Color for the annotation text
     */
    function getContrastColor(value, zValues) {
        // Flatten the 2D array and filter out null values
        const flatValues = zValues.flat().filter(v => v !== null);
        
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
        createProtocolTimeline,
        createProtocolComparisonHeatmap,
        createProtocolParallelCoordinates,
        createCoolingRateSurface
    };
})();

// Export the module
if (typeof module !== 'undefined' && module.exports) {
    module.exports = ProtocolVisualizer;
} else {
    window.ProtocolVisualizer = ProtocolVisualizer;
}
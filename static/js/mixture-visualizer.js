/**
 * CryoProtect Analyzer - Mixture Visualizer
 * 
 * This module provides functions for interactive visualization of mixture compositions
 */

const MixtureVisualizer = (function() {
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
     * Create an interactive pie chart for mixture composition
     * 
     * @param {string} containerId - ID of the container element
     * @param {Object} mixture - Mixture object with components
     * @param {Object} options - Visualization options
     */
    function createCompositionPieChart(containerId, mixture, options = {}) {
        const container = document.getElementById(containerId);
        if (!container) {
            console.error(`Container element with ID ${containerId} not found`);
            return;
        }
        
        if (!mixture || !mixture.components || mixture.components.length === 0) {
            container.innerHTML = '<div class="alert alert-warning">No mixture components found</div>';
            return;
        }
        
        // Default options
        const defaultOptions = {
            title: `Composition of ${mixture.name || 'Mixture'}`,
            showLegend: true,
            height: 400,
            donut: false,
            onClick: null
        };
        
        // Merge options
        const chartOptions = { ...defaultOptions, ...options };
        
        // Prepare data for chart
        const labels = mixture.components.map(component => component.name || `CID: ${component.cid}`);
        const values = mixture.components.map(component => component.concentration);
        const concentrationUnits = mixture.components.map(component => component.concentration_unit || '%');
        
        // Create hover template
        const hoverTemplate = '%{label}: %{value}%{text} (%{percent})<extra></extra>';
        
        // Create trace
        const trace = {
            labels: labels,
            values: values,
            type: 'pie',
            textinfo: 'percent',
            hovertemplate: hoverTemplate,
            text: concentrationUnits,
            hole: chartOptions.donut ? 0.4 : 0
        };
        
        // Create layout
        const layout = {
            title: chartOptions.title,
            height: chartOptions.height,
            showlegend: chartOptions.showLegend,
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
                filename: 'mixture_composition',
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
                const component = mixture.components[point.pointIndex];
                chartOptions.onClick(component, point.pointIndex);
            });
        }
    }
    
    /**
     * Create an interactive stacked bar chart for mixture composition
     * 
     * @param {string} containerId - ID of the container element
     * @param {Array} mixtures - Array of mixture objects
     * @param {Object} options - Visualization options
     */
    function createCompositionBarChart(containerId, mixtures, options = {}) {
        const container = document.getElementById(containerId);
        if (!container) {
            console.error(`Container element with ID ${containerId} not found`);
            return;
        }
        
        if (!mixtures || mixtures.length === 0) {
            container.innerHTML = '<div class="alert alert-warning">No mixtures found</div>';
            return;
        }
        
        // Default options
        const defaultOptions = {
            title: 'Mixture Compositions',
            height: 500,
            barMode: 'stack',
            onClick: null
        };
        
        // Merge options
        const chartOptions = { ...defaultOptions, ...options };
        
        // Get all unique components across all mixtures
        const allComponents = new Set();
        mixtures.forEach(mixture => {
            if (mixture.components) {
                mixture.components.forEach(component => {
                    allComponents.add(component.name || `CID: ${component.cid}`);
                });
            }
        });
        
        // Convert to array
        const componentNames = Array.from(allComponents);
        
        // Create traces for each component
        const traces = componentNames.map((componentName, index) => {
            const trace = {
                x: mixtures.map(mixture => mixture.name || `Mixture ${mixture.id}`),
                y: mixtures.map(mixture => {
                    if (!mixture.components) return 0;
                    const component = mixture.components.find(c => (c.name || `CID: ${c.cid}`) === componentName);
                    return component ? component.concentration : 0;
                }),
                name: componentName,
                type: 'bar',
                marker: {
                    color: COLORS[index % COLORS.length]
                },
                hovertemplate: '%{x}<br>%{y}%<br>%{fullData.name}<extra></extra>'
            };
            return trace;
        });
        
        // Create layout
        const layout = {
            title: chartOptions.title,
            barmode: chartOptions.barMode,
            height: chartOptions.height,
            xaxis: {
                title: 'Mixtures'
            },
            yaxis: {
                title: 'Concentration (%)',
                range: [0, 100]
            },
            legend: {
                title: {
                    text: 'Components'
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
                filename: 'mixture_compositions',
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
                const mixtureIndex = point.pointIndex;
                const componentName = point.fullData.name;
                const mixture = mixtures[mixtureIndex];
                const component = mixture.components.find(c => (c.name || `CID: ${c.cid}`) === componentName);
                
                chartOptions.onClick(mixture, component, mixtureIndex);
            });
        }
    }
    
    /**
     * Create an interactive treemap for mixture composition
     * 
     * @param {string} containerId - ID of the container element
     * @param {Object} mixture - Mixture object with components
     * @param {Object} options - Visualization options
     */
    function createCompositionTreemap(containerId, mixture, options = {}) {
        const container = document.getElementById(containerId);
        if (!container) {
            console.error(`Container element with ID ${containerId} not found`);
            return;
        }
        
        if (!mixture || !mixture.components || mixture.components.length === 0) {
            container.innerHTML = '<div class="alert alert-warning">No mixture components found</div>';
            return;
        }
        
        // Default options
        const defaultOptions = {
            title: `Composition of ${mixture.name || 'Mixture'}`,
            height: 500,
            colorScale: 'Viridis',
            onClick: null
        };
        
        // Merge options
        const chartOptions = { ...defaultOptions, ...options };
        
        // Prepare data for treemap
        const labels = mixture.components.map(component => component.name || `CID: ${component.cid}`);
        const values = mixture.components.map(component => component.concentration);
        const parents = mixture.components.map(() => mixture.name || 'Mixture');
        
        // Add root node
        labels.unshift(mixture.name || 'Mixture');
        values.unshift(0);
        parents.unshift('');
        
        // Create trace
        const trace = {
            type: 'treemap',
            labels: labels,
            parents: parents,
            values: values,
            textinfo: 'label+value+percent parent',
            hovertemplate: '<b>%{label}</b><br>Concentration: %{value}%<br>Percentage: %{percentParent:.1%}<extra></extra>',
            marker: {
                colorscale: chartOptions.colorScale
            }
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
                filename: 'mixture_composition_treemap',
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
                // Skip the root node
                if (point.pointIndex > 0) {
                    const component = mixture.components[point.pointIndex - 1];
                    chartOptions.onClick(component, point.pointIndex - 1);
                }
            });
        }
    }
    
    /**
     * Create an interactive radar chart for mixture properties
     * 
     * @param {string} containerId - ID of the container element
     * @param {Array} mixtures - Array of mixture objects
     * @param {Array} properties - Array of property names to display
     * @param {Object} options - Visualization options
     */
    function createMixturePropertiesRadar(containerId, mixtures, properties, options = {}) {
        const container = document.getElementById(containerId);
        if (!container) {
            console.error(`Container element with ID ${containerId} not found`);
            return;
        }
        
        if (!mixtures || mixtures.length === 0) {
            container.innerHTML = '<div class="alert alert-warning">No mixtures found</div>';
            return;
        }
        
        if (!properties || properties.length < 3) {
            container.innerHTML = '<div class="alert alert-warning">At least 3 properties are required for a radar chart</div>';
            return;
        }
        
        // Default options
        const defaultOptions = {
            title: 'Mixture Properties Comparison',
            height: 500,
            maxValue: 100,
            onClick: null
        };
        
        // Merge options
        const chartOptions = { ...defaultOptions, ...options };
        
        // Create traces for each mixture
        const traces = mixtures.map((mixture, index) => {
            const trace = {
                type: 'scatterpolar',
                r: properties.map(prop => mixture.properties && mixture.properties[prop] ? mixture.properties[prop] : 0),
                theta: properties,
                fill: 'toself',
                name: mixture.name || `Mixture ${mixture.id}`,
                line: {
                    color: COLORS[index % COLORS.length],
                    width: 2
                },
                fillcolor: `${COLORS[index % COLORS.length]}40`,
                hovertemplate: '%{theta}: %{r}<extra>%{fullData.name}</extra>'
            };
            return trace;
        });
        
        // Create layout
        const layout = {
            title: chartOptions.title,
            height: chartOptions.height,
            polar: {
                radialaxis: {
                    visible: true,
                    range: [0, chartOptions.maxValue]
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
                filename: 'mixture_properties_radar',
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
                const mixtureIndex = point.curveNumber;
                const propertyName = point.theta;
                const propertyValue = point.r;
                const mixture = mixtures[mixtureIndex];
                
                chartOptions.onClick(mixture, propertyName, propertyValue, mixtureIndex);
            });
        }
    }
    
    /**
     * Create an interactive heatmap for mixture compatibility
     * 
     * @param {string} containerId - ID of the container element
     * @param {Object} compatibilityData - Compatibility data with component pairs
     * @param {Object} options - Visualization options
     */
    function createCompatibilityHeatmap(containerId, compatibilityData, options = {}) {
        const container = document.getElementById(containerId);
        if (!container) {
            console.error(`Container element with ID ${containerId} not found`);
            return;
        }
        
        if (!compatibilityData || !compatibilityData.component_pairs || compatibilityData.component_pairs.length === 0) {
            container.innerHTML = '<div class="alert alert-warning">No compatibility data found</div>';
            return;
        }
        
        // Default options
        const defaultOptions = {
            title: 'Component Compatibility Matrix',
            height: 600,
            colorScale: [
                [0, '#e74c3c'],    // Red for incompatible
                [0.4, '#f39c12'],  // Orange for low compatibility
                [0.6, '#f1c40f'],  // Yellow for medium compatibility
                [0.8, '#2ecc71'],  // Green for high compatibility
                [1, '#27ae60']     // Dark green for fully compatible
            ],
            annotate: true,
            onClick: null
        };
        
        // Merge options
        const chartOptions = { ...defaultOptions, ...options };
        
        // Extract unique components
        const components = new Set();
        compatibilityData.component_pairs.forEach(pair => {
            components.add(pair.component1.name);
            components.add(pair.component2.name);
        });
        
        // Convert to array
        const componentNames = Array.from(components);
        
        // Create compatibility matrix
        const matrix = Array(componentNames.length).fill().map(() => Array(componentNames.length).fill(null));
        
        // Fill matrix with compatibility scores
        compatibilityData.component_pairs.forEach(pair => {
            const i = componentNames.indexOf(pair.component1.name);
            const j = componentNames.indexOf(pair.component2.name);
            
            if (i !== -1 && j !== -1) {
                matrix[i][j] = pair.compatibility_score;
                // Make the matrix symmetric
                matrix[j][i] = pair.compatibility_score;
            }
        });
        
        // Set diagonal to 1 (self-compatibility)
        for (let i = 0; i < componentNames.length; i++) {
            matrix[i][i] = 1;
        }
        
        // Create trace
        const trace = {
            z: matrix,
            x: componentNames,
            y: componentNames,
            type: 'heatmap',
            colorscale: chartOptions.colorScale,
            zmin: 0,
            zmax: 1,
            hovertemplate: '%{y} + %{x}<br>Compatibility: %{z:.2f}<extra></extra>'
        };
        
        // Create layout
        const layout = {
            title: chartOptions.title,
            height: chartOptions.height,
            xaxis: {
                title: 'Component',
                tickangle: -45
            },
            yaxis: {
                title: 'Component'
            }
        };
        
        // Add annotations if requested
        if (chartOptions.annotate) {
            const annotations = [];
            
            for (let i = 0; i < componentNames.length; i++) {
                for (let j = 0; j < componentNames.length; j++) {
                    if (matrix[i][j] !== null) {
                        annotations.push({
                            x: componentNames[j],
                            y: componentNames[i],
                            text: matrix[i][j].toFixed(2),
                            showarrow: false,
                            font: {
                                color: matrix[i][j] > 0.5 ? 'white' : 'black'
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
                filename: 'compatibility_matrix',
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
                const component1 = componentNames[point.pointIndex[0]];
                const component2 = componentNames[point.pointIndex[1]];
                const compatibilityScore = point.z;
                
                // Find the pair in the original data
                const pair = compatibilityData.component_pairs.find(p => 
                    (p.component1.name === component1 && p.component2.name === component2) ||
                    (p.component1.name === component2 && p.component2.name === component1)
                );
                
                chartOptions.onClick(pair, compatibilityScore);
            });
        }
    }
    
    /**
     * Create an interactive 3D scatter plot for mixture property space
     * 
     * @param {string} containerId - ID of the container element
     * @param {Array} mixtures - Array of mixture objects
     * @param {string} xProperty - Property name for x-axis
     * @param {string} yProperty - Property name for y-axis
     * @param {string} zProperty - Property name for z-axis
     * @param {Object} options - Visualization options
     */
    function createMixturePropertySpace(containerId, mixtures, xProperty, yProperty, zProperty, options = {}) {
        const container = document.getElementById(containerId);
        if (!container) {
            console.error(`Container element with ID ${containerId} not found`);
            return;
        }
        
        if (!mixtures || mixtures.length === 0) {
            container.innerHTML = '<div class="alert alert-warning">No mixtures found</div>';
            return;
        }
        
        // Default options
        const defaultOptions = {
            title: 'Mixture Property Space',
            height: 700,
            colorProperty: 'Total Score',
            sizeProperty: null,
            minSize: 6,
            maxSize: 15,
            onClick: null
        };
        
        // Merge options
        const chartOptions = { ...defaultOptions, ...options };
        
        // Prepare data for chart
        const xValues = mixtures.map(mixture => 
            mixture.properties && mixture.properties[xProperty] ? mixture.properties[xProperty] : 0
        );
        
        const yValues = mixtures.map(mixture => 
            mixture.properties && mixture.properties[yProperty] ? mixture.properties[yProperty] : 0
        );
        
        const zValues = mixtures.map(mixture => 
            mixture.properties && mixture.properties[zProperty] ? mixture.properties[zProperty] : 0
        );
        
        const labels = mixtures.map(mixture => mixture.name || `Mixture ${mixture.id}`);
        
        // Create trace
        const trace = {
            x: xValues,
            y: yValues,
            z: zValues,
            text: labels,
            type: 'scatter3d',
            mode: 'markers',
            marker: {
                size: 8,
                color: COLORS[0]
            },
            hovertemplate: `%{text}<br>${xProperty}: %{x}<br>${yProperty}: %{y}<br>${zProperty}: %{z}<extra></extra>`
        };
        
        // Apply color mapping if colorProperty is provided
        if (chartOptions.colorProperty && mixtures[0] && 
            mixtures[0].properties && mixtures[0].properties[chartOptions.colorProperty] !== undefined) {
            
            const colorValues = mixtures.map(mixture => 
                mixture.properties && mixture.properties[chartOptions.colorProperty] ? 
                mixture.properties[chartOptions.colorProperty] : 0
            );
            
            trace.marker.color = colorValues;
            trace.marker.colorscale = 'Viridis';
            trace.marker.colorbar = {
                title: chartOptions.colorProperty
            };
            
            // Update hover template
            trace.hovertemplate = `%{text}<br>${xProperty}: %{x}<br>${yProperty}: %{y}<br>${zProperty}: %{z}<br>${chartOptions.colorProperty}: %{marker.color}<extra></extra>`;
        }
        
        // Apply size mapping if sizeProperty is provided
        if (chartOptions.sizeProperty && mixtures[0] && 
            mixtures[0].properties && mixtures[0].properties[chartOptions.sizeProperty] !== undefined) {
            
            const sizeValues = mixtures.map(mixture => 
                mixture.properties && mixture.properties[chartOptions.sizeProperty] ? 
                mixture.properties[chartOptions.sizeProperty] : 0
            );
            
            const minValue = Math.min(...sizeValues);
            const maxValue = Math.max(...sizeValues);
            const range = maxValue - minValue;
            
            trace.marker.size = sizeValues.map(value => {
                if (range === 0) return chartOptions.minSize;
                const normalized = (value - minValue) / range;
                return chartOptions.minSize + normalized * (chartOptions.maxSize - chartOptions.minSize);
            });
            
            // Update hover template
            if (chartOptions.colorProperty) {
                trace.hovertemplate = `%{text}<br>${xProperty}: %{x}<br>${yProperty}: %{y}<br>${zProperty}: %{z}<br>${chartOptions.colorProperty}: %{marker.color}<br>${chartOptions.sizeProperty}: %{marker.size}<extra></extra>`;
            } else {
                trace.hovertemplate = `%{text}<br>${xProperty}: %{x}<br>${yProperty}: %{y}<br>${zProperty}: %{z}<br>${chartOptions.sizeProperty}: %{marker.size}<extra></extra>`;
            }
        }
        
        // Create layout
        const layout = {
            title: chartOptions.title,
            height: chartOptions.height,
            scene: {
                xaxis: {
                    title: xProperty
                },
                yaxis: {
                    title: yProperty
                },
                zaxis: {
                    title: zProperty
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
                filename: 'mixture_property_space',
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
                const mixtureIndex = point.pointIndex;
                const mixture = mixtures[mixtureIndex];
                
                chartOptions.onClick(mixture, mixtureIndex);
            });
        }
    }
    
    // Public API
    return {
        createCompositionPieChart,
        createCompositionBarChart,
        createCompositionTreemap,
        createMixturePropertiesRadar,
        createCompatibilityHeatmap,
        createMixturePropertySpace
    };
})();

// Export the module
if (typeof module !== 'undefined' && module.exports) {
    module.exports = MixtureVisualizer;
} else {
    window.MixtureVisualizer = MixtureVisualizer;
}
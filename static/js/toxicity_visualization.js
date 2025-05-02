/**
 * CryoProtect Analyzer - Toxicity Visualization Module
 * 
 * This module provides functions for visualizing toxicity data and unified scores
 * in the CryoProtect Analyzer UI. It includes:
 * - Toxicity endpoint visualization
 * - Unified score visualization with component breakdown
 * - Comparative toxicity visualization for mixtures
 * - Interactive charts and tooltips
 * 
 * Dependencies:
 * - charts.js - For chart rendering
 * - api.js - For API communication
 */

// Toxicity visualization namespace
const ToxicityViz = (function() {
    'use strict';

    // Color schemes for toxicity visualization
    const COLORS = {
        toxicity: {
            safe: '#4CAF50',      // Green for low toxicity
            moderate: '#FFC107',   // Amber for moderate toxicity
            high: '#F44336',       // Red for high toxicity
            unknown: '#9E9E9E'     // Gray for unknown/no data
        },
        endpoints: {
            'Acute Toxicity': '#E57373',
            'Carcinogenicity': '#F06292',
            'Developmental Toxicity': '#BA68C8',
            'Endocrine Disruption': '#9575CD',
            'Genotoxicity': '#7986CB',
            'Hepatotoxicity': '#64B5F6',
            'Neurotoxicity': '#4FC3F7',
            'Reproductive Toxicity': '#4DD0E1',
            'Skin Sensitization': '#4DB6AC',
            'Cytotoxicity': '#81C784',
            'default': '#A1887F'
        },
        unified: {
            efficacy: '#42A5F5',    // Blue for efficacy
            toxicity: '#EF5350',    // Red for toxicity
            tg: '#66BB6A'           // Green for glass transition
        }
    };

    // Toxicity score thresholds
    const TOXICITY_THRESHOLDS = {
        safe: 30,       // Scores below this are considered safe
        moderate: 70    // Scores below this but above safe are moderate, above is high
    };

    /**
     * Get color for toxicity score
     * @param {number} score - Toxicity score (0-100)
     * @returns {string} - Color hex code
     */
    function getToxicityColor(score) {
        if (score === null || score === undefined) {
            return COLORS.toxicity.unknown;
        }
        
        if (score <= TOXICITY_THRESHOLDS.safe) {
            return COLORS.toxicity.safe;
        } else if (score <= TOXICITY_THRESHOLDS.moderate) {
            return COLORS.toxicity.moderate;
        } else {
            return COLORS.toxicity.high;
        }
    }

    /**
     * Get color for toxicity endpoint
     * @param {string} endpoint - Toxicity endpoint name
     * @returns {string} - Color hex code
     */
    function getEndpointColor(endpoint) {
        return COLORS.endpoints[endpoint] || COLORS.endpoints.default;
    }

    /**
     * Format toxicity score with appropriate label
     * @param {number} score - Toxicity score (0-100)
     * @returns {object} - Formatted score with label and color
     */
    function formatToxicityScore(score) {
        if (score === null || score === undefined) {
            return { value: 'Unknown', label: 'Unknown', color: COLORS.toxicity.unknown };
        }
        
        const formattedScore = Math.round(score);
        let label, color;
        
        if (score <= TOXICITY_THRESHOLDS.safe) {
            label = 'Low Risk';
            color = COLORS.toxicity.safe;
        } else if (score <= TOXICITY_THRESHOLDS.moderate) {
            label = 'Moderate Risk';
            color = COLORS.toxicity.moderate;
        } else {
            label = 'High Risk';
            color = COLORS.toxicity.high;
        }
        
        return { value: formattedScore, label, color };
    }

    /**
     * Create toxicity score badge
     * @param {number} score - Toxicity score (0-100)
     * @param {string} size - Badge size ('sm', 'md', 'lg')
     * @returns {string} - HTML for score badge
     */
    function createScoreBadge(score, size = 'md') {
        const formatted = formatToxicityScore(score);
        const sizeClass = `badge-${size}`;
        
        return `
            <span class="badge toxicity-badge ${sizeClass}" 
                  style="background-color: ${formatted.color};" 
                  title="Toxicity Score: ${formatted.value}">
                ${formatted.label}: ${formatted.value}
            </span>
        `;
    }

    /**
     * Create toxicity endpoint visualization
     * @param {Array} endpointScores - Array of endpoint scores
     * @param {string} containerId - ID of container element
     */
    function createEndpointVisualization(endpointScores, containerId) {
        const container = document.getElementById(containerId);
        if (!container) return;
        
        // Clear container
        container.innerHTML = '';
        
        if (!endpointScores || endpointScores.length === 0) {
            container.innerHTML = '<div class="alert alert-info">No endpoint data available</div>';
            return;
        }
        
        // Create endpoint visualization
        const endpointContainer = document.createElement('div');
        endpointContainer.className = 'toxicity-endpoints';
        
        // Sort endpoints by score (highest first)
        const sortedEndpoints = [...endpointScores].sort((a, b) => b.score_value - a.score_value);
        
        sortedEndpoints.forEach(endpoint => {
            const formatted = formatToxicityScore(endpoint.score_value);
            const endpointColor = getEndpointColor(endpoint.endpoint_name);
            
            const endpointEl = document.createElement('div');
            endpointEl.className = 'toxicity-endpoint-item';
            endpointEl.innerHTML = `
                <div class="endpoint-header">
                    <span class="endpoint-name" style="border-left: 4px solid ${endpointColor}; padding-left: 8px;">
                        ${endpoint.endpoint_name}
                    </span>
                    <span class="endpoint-score" style="color: ${formatted.color};">
                        ${formatted.value}
                    </span>
                </div>
                <div class="endpoint-bar-container">
                    <div class="endpoint-bar" style="width: ${endpoint.score_value}%; background-color: ${formatted.color};"></div>
                </div>
                <div class="endpoint-details">
                    <small>Confidence: ${Math.round(endpoint.confidence * 100)}%</small>
                </div>
            `;
            
            endpointContainer.appendChild(endpointEl);
        });
        
        container.appendChild(endpointContainer);
    }

    /**
     * Create unified score visualization
     * @param {object} unifiedScore - Unified score data
     * @param {string} containerId - ID of container element
     */
    function createUnifiedScoreVisualization(unifiedScore, containerId) {
        const container = document.getElementById(containerId);
        if (!container) return;
        
        // Clear container
        container.innerHTML = '';
        
        if (!unifiedScore || !unifiedScore.score) {
            container.innerHTML = '<div class="alert alert-info">No unified score data available</div>';
            return;
        }
        
        // Create unified score visualization
        const scoreContainer = document.createElement('div');
        scoreContainer.className = 'unified-score-container';
        
        // Main score display
        const scoreDisplay = document.createElement('div');
        scoreDisplay.className = 'unified-score-display';
        
        const score = Math.round(unifiedScore.score);
        const scoreColor = getScoreColor(score);
        
        scoreDisplay.innerHTML = `
            <div class="score-circle" style="background: conic-gradient(${scoreColor} ${score}%, #e0e0e0 0);">
                <div class="score-value">${score}</div>
            </div>
            <div class="score-label">
                <h4>Unified Score</h4>
                <div class="score-context">${unifiedScore.application_context}</div>
            </div>
        `;
        
        scoreContainer.appendChild(scoreDisplay);
        
        // Component scores
        if (unifiedScore.component_scores) {
            const componentContainer = document.createElement('div');
            componentContainer.className = 'component-scores';
            
            // Create component score bars
            const components = [
                {
                    name: 'Efficacy',
                    key: 'efficacy',
                    color: COLORS.unified.efficacy,
                    score: unifiedScore.component_scores.efficacy.score,
                    weight: unifiedScore.component_scores.efficacy.weight
                },
                {
                    name: 'Toxicity',
                    key: 'toxicity',
                    color: COLORS.unified.toxicity,
                    score: unifiedScore.component_scores.toxicity.score,
                    weight: unifiedScore.component_scores.toxicity.weight
                },
                {
                    name: 'Glass Transition',
                    key: 'glass_transition',
                    color: COLORS.unified.tg,
                    score: unifiedScore.component_scores.glass_transition.score,
                    weight: unifiedScore.component_scores.glass_transition.weight
                }
            ];
            
            components.forEach(component => {
                const componentEl = document.createElement('div');
                componentEl.className = 'component-score-item';
                
                componentEl.innerHTML = `
                    <div class="component-header">
                        <span class="component-name">${component.name}</span>
                        <span class="component-score">${Math.round(component.score)}</span>
                    </div>
                    <div class="component-bar-container">
                        <div class="component-bar" style="width: ${component.score}%; background-color: ${component.color};"></div>
                    </div>
                    <div class="component-details">
                        <small>Weight: ${Math.round(component.weight * 100)}%</small>
                    </div>
                `;
                
                componentContainer.appendChild(componentEl);
            });
            
            scoreContainer.appendChild(componentContainer);
        }
        
        // Add application context selector
        const contextSelector = document.createElement('div');
        contextSelector.className = 'context-selector';
        contextSelector.innerHTML = `
            <label for="application-context">Application Context:</label>
            <select id="application-context" class="form-control">
                <option value="general" ${unifiedScore.application_context === 'general' ? 'selected' : ''}>General</option>
                <option value="cell_preservation" ${unifiedScore.application_context === 'cell_preservation' ? 'selected' : ''}>Cell Preservation</option>
                <option value="organ_preservation" ${unifiedScore.application_context === 'organ_preservation' ? 'selected' : ''}>Organ Preservation</option>
                <option value="long_term_storage" ${unifiedScore.application_context === 'long_term_storage' ? 'selected' : ''}>Long-term Storage</option>
                <option value="sensitive_tissues" ${unifiedScore.application_context === 'sensitive_tissues' ? 'selected' : ''}>Sensitive Tissues</option>
            </select>
        `;
        
        scoreContainer.appendChild(contextSelector);
        
        // Add event listener for context change
        setTimeout(() => {
            const select = document.getElementById('application-context');
            if (select) {
                select.addEventListener('change', function() {
                    const entityId = container.dataset.entityId;
                    const entityType = container.dataset.entityType;
                    
                    if (entityId && entityType) {
                        loadUnifiedScore(entityId, entityType, this.value, containerId);
                    }
                });
            }
        }, 0);
        
        container.appendChild(scoreContainer);
    }

    /**
     * Get color for unified score
     * @param {number} score - Unified score (0-100)
     * @returns {string} - Color hex code
     */
    function getScoreColor(score) {
        // Create gradient from red (0) to yellow (50) to green (100)
        if (score <= 50) {
            // Red to yellow gradient
            const r = 255;
            const g = Math.round((score / 50) * 255);
            return `rgb(${r}, ${g}, 0)`;
        } else {
            // Yellow to green gradient
            const r = Math.round(255 - ((score - 50) / 50) * 255);
            const g = 255;
            return `rgb(${r}, ${g}, 0)`;
        }
    }

    /**
     * Create comparative toxicity visualization for mixtures
     * @param {object} mixtureData - Mixture toxicity data
     * @param {string} containerId - ID of container element
     */
    function createMixtureToxicityVisualization(mixtureData, containerId) {
        const container = document.getElementById(containerId);
        if (!container) return;
        
        // Clear container
        container.innerHTML = '';
        
        if (!mixtureData || !mixtureData.component_scores || Object.keys(mixtureData.component_scores).length === 0) {
            container.innerHTML = '<div class="alert alert-info">No mixture toxicity data available</div>';
            return;
        }
        
        // Create mixture visualization
        const mixtureContainer = document.createElement('div');
        mixtureContainer.className = 'mixture-toxicity-container';
        
        // Aggregate score
        if (mixtureData.aggregate_score) {
            const aggregateEl = document.createElement('div');
            aggregateEl.className = 'aggregate-score';
            
            const formatted = formatToxicityScore(mixtureData.aggregate_score.aggregate_toxicity_score);
            
            aggregateEl.innerHTML = `
                <h4>Mixture Toxicity Score</h4>
                <div class="score-badge" style="background-color: ${formatted.color};">
                    ${Math.round(mixtureData.aggregate_score.aggregate_toxicity_score)}
                </div>
                <div class="score-label">${formatted.label}</div>
                <div class="calculation-method">
                    <small>Method: ${mixtureData.aggregate_score.calculation_method}</small>
                </div>
            `;
            
            mixtureContainer.appendChild(aggregateEl);
        }
        
        // Component scores
        const componentsEl = document.createElement('div');
        componentsEl.className = 'component-toxicity-scores';
        componentsEl.innerHTML = '<h4>Component Toxicity Scores</h4>';
        
        // Get component details
        API.getMixtureComponents(mixtureData.mixture_id)
            .then(components => {
                const componentMap = {};
                components.forEach(component => {
                    componentMap[component.molecule_id] = component;
                });
                
                // Create component list
                const componentList = document.createElement('ul');
                componentList.className = 'component-list';
                
                Object.entries(mixtureData.component_scores).forEach(([moleculeId, scores]) => {
                    const component = componentMap[moleculeId] || { name: 'Unknown Component' };
                    const componentScore = scores.find(s => s.score_type === 'overall_toxicity');
                    
                    if (componentScore) {
                        const formatted = formatToxicityScore(componentScore.score_value);
                        
                        const componentEl = document.createElement('li');
                        componentEl.className = 'component-item';
                        componentEl.innerHTML = `
                            <div class="component-name">${component.name || 'Unknown'}</div>
                            <div class="component-concentration">${component.concentration} ${component.concentration_unit}</div>
                            <div class="component-score" style="color: ${formatted.color};">${Math.round(componentScore.score_value)}</div>
                        `;
                        
                        componentList.appendChild(componentEl);
                    }
                });
                
                componentsEl.appendChild(componentList);
            })
            .catch(error => {
                console.error('Error loading mixture components:', error);
                componentsEl.innerHTML += '<div class="alert alert-danger">Error loading component details</div>';
            });
        
        mixtureContainer.appendChild(componentsEl);
        container.appendChild(mixtureContainer);
    }

    /**
     * Create toxicity radar chart
     * @param {Array} endpointScores - Array of endpoint scores
     * @param {string} canvasId - ID of canvas element
     */
    function createToxicityRadarChart(endpointScores, canvasId) {
        const canvas = document.getElementById(canvasId);
        if (!canvas || !endpointScores || endpointScores.length === 0) return;
        
        // Prepare data for radar chart
        const labels = endpointScores.map(endpoint => endpoint.endpoint_name);
        const data = endpointScores.map(endpoint => endpoint.score_value);
        const backgroundColor = 'rgba(255, 99, 132, 0.2)';
        const borderColor = 'rgb(255, 99, 132)';
        
        // Create chart
        new Chart(canvas, {
            type: 'radar',
            data: {
                labels: labels,
                datasets: [{
                    label: 'Toxicity Profile',
                    data: data,
                    backgroundColor: backgroundColor,
                    borderColor: borderColor,
                    borderWidth: 1
                }]
            },
            options: {
                scales: {
                    r: {
                        angleLines: {
                            display: true
                        },
                        suggestedMin: 0,
                        suggestedMax: 100
                    }
                }
            }
        });
    }

    /**
     * Create unified score comparison chart
     * @param {Array} scores - Array of unified scores for different entities
     * @param {string} canvasId - ID of canvas element
     */
    function createUnifiedScoreComparisonChart(scores, canvasId) {
        const canvas = document.getElementById(canvasId);
        if (!canvas || !scores || scores.length === 0) return;
        
        // Prepare data for bar chart
        const labels = scores.map(item => item.name || item.id);
        
        const datasets = [
            {
                label: 'Unified Score',
                data: scores.map(item => item.score),
                backgroundColor: scores.map(item => getScoreColor(item.score)),
                borderColor: 'rgba(0, 0, 0, 0.1)',
                borderWidth: 1
            },
            {
                label: 'Efficacy',
                data: scores.map(item => item.component_scores?.efficacy?.score || 0),
                backgroundColor: COLORS.unified.efficacy,
                borderColor: 'rgba(0, 0, 0, 0.1)',
                borderWidth: 1
            },
            {
                label: 'Toxicity',
                data: scores.map(item => item.component_scores?.toxicity?.score || 0),
                backgroundColor: COLORS.unified.toxicity,
                borderColor: 'rgba(0, 0, 0, 0.1)',
                borderWidth: 1
            },
            {
                label: 'Glass Transition',
                data: scores.map(item => item.component_scores?.glass_transition?.score || 0),
                backgroundColor: COLORS.unified.tg,
                borderColor: 'rgba(0, 0, 0, 0.1)',
                borderWidth: 1
            }
        ];
        
        // Create chart
        new Chart(canvas, {
            type: 'bar',
            data: {
                labels: labels,
                datasets: datasets
            },
            options: {
                scales: {
                    y: {
                        beginAtZero: true,
                        max: 100
                    }
                }
            }
        });
    }

    /**
     * Load toxicity data for a molecule
     * @param {string} moleculeId - Molecule ID
     * @param {string} containerId - ID of container element
     */
    function loadMoleculeToxicityData(moleculeId, containerId) {
        const container = document.getElementById(containerId);
        if (!container) return;
        
        // Show loading indicator
        container.innerHTML = '<div class="loading">Loading toxicity data...</div>';
        
        // Load toxicity scores
        API.getToxicityScores(moleculeId)
            .then(data => {
                if (!data || !data.scores || data.scores.length === 0) {
                    container.innerHTML = '<div class="alert alert-info">No toxicity data available for this molecule</div>';
                    return;
                }
                
                // Prepare endpoint scores for visualization
                const endpointScores = data.scores
                    .filter(score => score.endpoint_id && score.endpoint)
                    .map(score => ({
                        endpoint_id: score.endpoint_id,
                        endpoint_name: score.endpoint.name,
                        score_value: score.score_value,
                        confidence: score.confidence
                    }));
                
                // Create visualization
                createEndpointVisualization(endpointScores, containerId);
                
                // If canvas exists, create radar chart
                const canvasId = `${containerId}-chart`;
                const canvas = document.getElementById(canvasId);
                if (canvas) {
                    createToxicityRadarChart(endpointScores, canvasId);
                }
            })
            .catch(error => {
                console.error('Error loading toxicity data:', error);
                container.innerHTML = '<div class="alert alert-danger">Error loading toxicity data</div>';
            });
    }

    /**
     * Load unified score for an entity
     * @param {string} entityId - Entity ID (molecule or mixture)
     * @param {string} entityType - Entity type ('molecule' or 'mixture')
     * @param {string} context - Application context
     * @param {string} containerId - ID of container element
     */
    function loadUnifiedScore(entityId, entityType, context, containerId) {
        const container = document.getElementById(containerId);
        if (!container) return;
        
        // Set data attributes for context selector
        container.dataset.entityId = entityId;
        container.dataset.entityType = entityType;
        
        // Show loading indicator
        container.innerHTML = '<div class="loading">Loading unified score...</div>';
        
        // Load unified score
        const apiMethod = entityType === 'molecule' ? API.getUnifiedMoleculeScore : API.getUnifiedMixtureScore;
        
        apiMethod(entityId, context)
            .then(data => {
                if (!data || !data.score) {
                    container.innerHTML = '<div class="alert alert-info">No unified score available</div>';
                    return;
                }
                
                // Create visualization
                createUnifiedScoreVisualization(data, containerId);
            })
            .catch(error => {
                console.error('Error loading unified score:', error);
                container.innerHTML = '<div class="alert alert-danger">Error loading unified score</div>';
            });
    }

    /**
     * Load mixture toxicity data
     * @param {string} mixtureId - Mixture ID
     * @param {string} containerId - ID of container element
     */
    function loadMixtureToxicityData(mixtureId, containerId) {
        const container = document.getElementById(containerId);
        if (!container) return;
        
        // Show loading indicator
        container.innerHTML = '<div class="loading">Loading mixture toxicity data...</div>';
        
        // Load mixture toxicity scores
        API.getMixtureToxicityScores(mixtureId)
            .then(data => {
                if (!data || !data.component_scores || Object.keys(data.component_scores).length === 0) {
                    container.innerHTML = '<div class="alert alert-info">No toxicity data available for this mixture</div>';
                    return;
                }
                
                // Create visualization
                createMixtureToxicityVisualization(data, containerId);
            })
            .catch(error => {
                console.error('Error loading mixture toxicity data:', error);
                container.innerHTML = '<div class="alert alert-danger">Error loading mixture toxicity data</div>';
            });
    }

    /**
     * Initialize toxicity visualization for a page
     */
    function initializePage() {
        // Initialize molecule toxicity visualizations
        document.querySelectorAll('[data-toxicity-molecule-id]').forEach(container => {
            const moleculeId = container.dataset.toxicityMoleculeId;
            if (moleculeId) {
                loadMoleculeToxicityData(moleculeId, container.id);
            }
        });
        
        // Initialize mixture toxicity visualizations
        document.querySelectorAll('[data-toxicity-mixture-id]').forEach(container => {
            const mixtureId = container.dataset.toxicityMixtureId;
            if (mixtureId) {
                loadMixtureToxicityData(mixtureId, container.id);
            }
        });
        
        // Initialize unified score visualizations
        document.querySelectorAll('[data-unified-score-entity]').forEach(container => {
            const entityId = container.dataset.unifiedScoreEntity;
            const entityType = container.dataset.unifiedScoreType || 'molecule';
            const context = container.dataset.unifiedScoreContext || 'general';
            
            if (entityId) {
                loadUnifiedScore(entityId, entityType, context, container.id);
            }
        });
    }

    // Public API
    return {
        createScoreBadge,
        createEndpointVisualization,
        createUnifiedScoreVisualization,
        createMixtureToxicityVisualization,
        createToxicityRadarChart,
        createUnifiedScoreComparisonChart,
        loadMoleculeToxicityData,
        loadUnifiedScore,
        loadMixtureToxicityData,
        initializePage,
        formatToxicityScore,
        getToxicityColor,
        getEndpointColor,
        getScoreColor,
        // Aliases for prompt compatibility
        loadMoleculeToxicity: loadMoleculeToxicityData,
        renderToxicityData: createEndpointVisualization,
        renderToxicityScores: createScoreBadge, // or another suitable function
        renderUnifiedScore: createUnifiedScoreVisualization,
        renderError: function(container, message) {
            container.innerHTML = `<div class="error-message">${message}</div>`;
        }
    };
})();

// Add API methods for toxicity data
if (typeof API !== 'undefined') {
    /**
     * Get toxicity data for a molecule
     * @param {string} moleculeId - Molecule ID
     * @param {object} options - Query options
     * @returns {Promise} - Promise resolving to toxicity data
     */
    API.getToxicityData = function(moleculeId, options = {}) {
        const queryParams = new URLSearchParams();
        
        if (options.assayId) queryParams.append('assay_id', options.assayId);
        if (options.endpoint) queryParams.append('endpoint', options.endpoint);
        if (options.activeOnly) queryParams.append('active_only', options.activeOnly);
        
        const queryString = queryParams.toString() ? `?${queryParams.toString()}` : '';
        
        return this.get(`/api/toxicity/molecule/${moleculeId}${queryString}`);
    };
    
    /**
     * Get toxicity scores for a molecule
     * @param {string} moleculeId - Molecule ID
     * @returns {Promise} - Promise resolving to toxicity scores
     */
    API.getToxicityScores = function(moleculeId) {
        return this.get(`/api/toxicity/scores/molecule/${moleculeId}`);
    };
    
    /**
     * Get toxicity scores for a mixture
     * @param {string} mixtureId - Mixture ID
     * @returns {Promise} - Promise resolving to mixture toxicity scores
     */
    API.getMixtureToxicityScores = function(mixtureId) {
        return this.get(`/api/toxicity/scores/mixture/${mixtureId}`);
    };
    
    /**
     * Get unified score for a molecule
     * @param {string} moleculeId - Molecule ID
     * @param {string} context - Application context
     * @param {string} algorithm - Algorithm to use
     * @param {boolean} recalculate - Whether to recalculate score
     * @returns {Promise} - Promise resolving to unified score
     */
    API.getUnifiedMoleculeScore = function(moleculeId, context = 'general', algorithm = 'random_forest', recalculate = false) {
        const queryParams = new URLSearchParams();
        queryParams.append('context', context);
        queryParams.append('algorithm', algorithm);
        queryParams.append('recalculate', recalculate);
        
        return this.get(`/api/toxicity/unified/molecule/${moleculeId}?${queryParams.toString()}`);
    };
    
    /**
     * Get unified score for a mixture
     * @param {string} mixtureId - Mixture ID
     * @param {string} context - Application context
     * @param {string} algorithm - Algorithm to use
     * @param {boolean} recalculate - Whether to recalculate score
     * @returns {Promise} - Promise resolving to unified score
     */
    API.getUnifiedMixtureScore = function(mixtureId, context = 'general', algorithm = 'random_forest', recalculate = false) {
        const queryParams = new URLSearchParams();
        queryParams.append('context', context);
        queryParams.append('algorithm', algorithm);
        queryParams.append('recalculate', recalculate);
        
        return this.get(`/api/toxicity/unified/mixture/${mixtureId}?${queryParams.toString()}`);
    };
    
    /**
     * Get available application contexts for unified scoring
     * @returns {Promise} - Promise resolving to application contexts
     */
    API.getUnifiedScoreContexts = function() {
        return this.get('/api/toxicity/unified/contexts');
    };
    
    /**
     * Calculate unified scores for multiple entities
     * @param {Array} entityIds - Array of entity IDs
     * @param {string} entityType - Entity type ('molecule' or 'mixture')
     * @param {string} context - Application context
     * @param {string} algorithm - Algorithm to use
     * @param {boolean} recalculate - Whether to recalculate scores
     * @returns {Promise} - Promise resolving to batch results
     */
    API.calculateBatchUnifiedScores = function(entityIds, entityType = 'molecule', context = 'general', algorithm = 'random_forest', recalculate = false) {
        return this.post('/api/toxicity/unified/batch', {
            entity_ids: entityIds,
            entity_type: entityType,
            application_context: context,
            algorithm: algorithm,
            recalculate: recalculate
        });
    };
}

// Initialize on DOM ready
document.addEventListener('DOMContentLoaded', function() {
    ToxicityViz.initializePage();
});
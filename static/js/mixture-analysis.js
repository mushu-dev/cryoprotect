/**
 * CryoProtect Analyzer - Mixture Analysis JavaScript Module
 * 
 * This module provides functions for interacting with the mixture analysis API endpoints.
 */

// Namespace for mixture analysis functions
const MixtureAnalysis = {
    /**
     * Get predicted properties for a mixture.
     * 
     * @param {string} mixtureId - ID of the mixture
     * @returns {Promise} - Promise resolving to the mixture properties
     */
    getProperties: async function(mixtureId) {
        try {
            const response = await fetch(`/api/mixtures/${mixtureId}/properties`);
            if (!response.ok) {
                throw new Error(`Error fetching mixture properties: ${response.statusText}`);
            }
            return await response.json();
        } catch (error) {
            console.error('Error in getProperties:', error);
            throw error;
        }
    },

    /**
     * Analyze compatibility between mixture components.
     * 
     * @param {string} mixtureId - ID of the mixture
     * @returns {Promise} - Promise resolving to the compatibility analysis
     */
    analyzeCompatibility: async function(mixtureId) {
        try {
            const response = await fetch(`/api/mixtures/${mixtureId}/compatibility`);
            if (!response.ok) {
                throw new Error(`Error analyzing mixture compatibility: ${response.statusText}`);
            }
            return await response.json();
        } catch (error) {
            console.error('Error in analyzeCompatibility:', error);
            throw error;
        }
    },

    /**
     * Analyze synergistic or antagonistic effects in a mixture.
     * 
     * @param {string} mixtureId - ID of the mixture
     * @returns {Promise} - Promise resolving to the synergy analysis
     */
    analyzeSynergy: async function(mixtureId) {
        try {
            const response = await fetch(`/api/mixtures/${mixtureId}/synergy`);
            if (!response.ok) {
                throw new Error(`Error analyzing mixture synergy: ${response.statusText}`);
            }
            return await response.json();
        } catch (error) {
            console.error('Error in analyzeSynergy:', error);
            throw error;
        }
    },

    /**
     * Optimize the composition of a mixture.
     * 
     * @param {string} mixtureId - ID of the mixture to optimize
     * @param {Object} options - Optimization options
     * @param {string} [options.targetProperty="Cryoprotection Score"] - Property to optimize
     * @param {number} [options.targetValue] - Target value for the property (if null, maximize the property)
     * @param {Object} [options.constraints] - Constraints for the optimization
     * @returns {Promise} - Promise resolving to the optimized composition
     */
    optimizeComposition: async function(mixtureId, options = {}) {
        try {
            const response = await fetch(`/api/mixtures/${mixtureId}/optimize`, {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json'
                },
                body: JSON.stringify({
                    target_property: options.targetProperty || "Cryoprotection Score",
                    target_value: options.targetValue,
                    constraints: options.constraints || {}
                })
            });
            if (!response.ok) {
                throw new Error(`Error optimizing mixture composition: ${response.statusText}`);
            }
            return await response.json();
        } catch (error) {
            console.error('Error in optimizeComposition:', error);
            throw error;
        }
    },

    /**
     * Optimize the composition of a mixture step by step.
     * 
     * @param {string} mixtureId - ID of the mixture to optimize
     * @param {Object} options - Optimization options
     * @param {string} [options.targetProperty="Cryoprotection Score"] - Property to optimize
     * @param {number} [options.targetValue] - Target value for the property (if null, maximize the property)
     * @param {Object} [options.constraints] - Constraints for the optimization
     * @returns {Promise} - Promise resolving to the optimized composition
     */
    optimizeStepByStep: async function(mixtureId, options = {}) {
        try {
            const response = await fetch(`/api/mixtures/${mixtureId}/optimize-step`, {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json'
                },
                body: JSON.stringify({
                    target_property: options.targetProperty || "Cryoprotection Score",
                    target_value: options.targetValue,
                    constraints: options.constraints || {}
                })
            });
            if (!response.ok) {
                throw new Error(`Error optimizing mixture composition: ${response.statusText}`);
            }
            return await response.json();
        } catch (error) {
            console.error('Error in optimizeStepByStep:', error);
            throw error;
        }
    },

    /**
     * Analyze a mixture and provide recommendations for improvement.
     * 
     * @param {string} mixtureId - ID of the mixture to analyze
     * @returns {Promise} - Promise resolving to the analysis and recommendations
     */
    analyzeMixture: async function(mixtureId) {
        try {
            const response = await fetch(`/api/mixtures/${mixtureId}/analyze`);
            if (!response.ok) {
                throw new Error(`Error analyzing mixture: ${response.statusText}`);
            }
            return await response.json();
        } catch (error) {
            console.error('Error in analyzeMixture:', error);
            throw error;
        }
    },

    /**
     * Recommend new components to add to a mixture.
     * 
     * @param {string} mixtureId - ID of the mixture to improve
     * @param {Object} options - Recommendation options
     * @param {string} [options.targetProperty] - Property to optimize (if null, optimize overall score)
     * @param {number} [options.targetValue] - Target value for the property (if null, maximize the property)
     * @param {number} [options.count=3] - Number of recommendations to provide
     * @returns {Promise} - Promise resolving to the component recommendations
     */
    recommendComponents: async function(mixtureId, options = {}) {
        try {
            const response = await fetch(`/api/mixtures/${mixtureId}/recommend-components`, {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json'
                },
                body: JSON.stringify({
                    target_property: options.targetProperty,
                    target_value: options.targetValue,
                    count: options.count || 3
                })
            });
            if (!response.ok) {
                throw new Error(`Error recommending components: ${response.statusText}`);
            }
            return await response.json();
        } catch (error) {
            console.error('Error in recommendComponents:', error);
            throw error;
        }
    },

    /**
     * Render a compatibility matrix for mixture components.
     * 
     * @param {Object} compatibilityData - Compatibility analysis data
     * @param {string} containerId - ID of the container element
     */
    renderCompatibilityMatrix: function(compatibilityData, containerId) {
        const container = document.getElementById(containerId);
        if (!container) {
            console.error(`Container element with ID ${containerId} not found`);
            return;
        }

        // Clear container
        container.innerHTML = '';

        // Get component pairs
        const pairs = compatibilityData.compatibility.component_pairs;
        if (!pairs || pairs.length === 0) {
            container.innerHTML = '<p>No component pairs to display</p>';
            return;
        }

        // Create responsive container
        const tableResponsive = document.createElement('div');
        tableResponsive.className = 'table-responsive';

        // Create table
        const table = document.createElement('table');
        table.className = 'compatibility-matrix table';

        // Create header row
        const headerRow = document.createElement('tr');
        headerRow.innerHTML = '<th>Component 1</th><th>Component 2</th><th>Compatibility</th><th>Issues</th>';
        table.appendChild(headerRow);

        // Create rows for each pair
        pairs.forEach(pair => {
            const row = document.createElement('tr');
            
            // Component 1
            const cell1 = document.createElement('td');
            cell1.textContent = pair.component1.name;
            row.appendChild(cell1);
            
            // Component 2
            const cell2 = document.createElement('td');
            cell2.textContent = pair.component2.name;
            row.appendChild(cell2);
            
            // Compatibility score
            const cell3 = document.createElement('td');
            const score = pair.compatibility_score || 0;
            cell3.textContent = score.toFixed(2);
            
            // Add color based on score
            if (score >= 0.8) {
                cell3.className = 'high-compatibility';
            } else if (score >= 0.6) {
                cell3.className = 'medium-compatibility';
            } else if (score >= 0.4) {
                cell3.className = 'low-compatibility';
            } else {
                cell3.className = 'incompatible';
            }
            
            row.appendChild(cell3);
            
            // Issues
            const cell4 = document.createElement('td');
            if (pair.issues && pair.issues.length > 0) {
                const ul = document.createElement('ul');
                pair.issues.forEach(issue => {
                    const li = document.createElement('li');
                    li.textContent = issue;
                    ul.appendChild(li);
                });
                cell4.appendChild(ul);
            } else {
                cell4.textContent = 'No issues';
            }
            row.appendChild(cell4);
            
            table.appendChild(row);
        });

        // Add table to responsive container
        tableResponsive.appendChild(table);
        
        // Add responsive container to main container
        container.appendChild(tableResponsive);
    },

    /**
     * Render a synergy analysis for a mixture.
     * 
     * @param {Object} synergyData - Synergy analysis data
     * @param {string} containerId - ID of the container element
     */
    renderSynergyAnalysis: function(synergyData, containerId) {
        const container = document.getElementById(containerId);
        if (!container) {
            console.error(`Container element with ID ${containerId} not found`);
            return;
        }

        // Clear container
        container.innerHTML = '';

        // Create synergy summary
        const summary = document.createElement('div');
        summary.className = 'synergy-summary';
        
        // Synergy type
        const synergyType = document.createElement('h3');
        synergyType.textContent = `Synergy Type: ${synergyData.synergy.synergy_type}`;
        summary.appendChild(synergyType);
        
        // Synergy score
        const synergyScore = document.createElement('p');
        synergyScore.textContent = `Synergy Score: ${synergyData.synergy.synergy_score.toFixed(2)}`;
        summary.appendChild(synergyScore);
        
        // Expected vs Actual
        const scoreComparison = document.createElement('p');
        scoreComparison.textContent = `Expected Score: ${synergyData.synergy.expected_score.toFixed(2)} | Actual Score: ${synergyData.synergy.actual_score.toFixed(2)}`;
        summary.appendChild(scoreComparison);
        
        container.appendChild(summary);

        // Component contributions
        const contributions = synergyData.synergy.component_contributions;
        if (contributions && contributions.length > 0) {
            const contributionsTitle = document.createElement('h3');
            contributionsTitle.textContent = 'Component Contributions';
            container.appendChild(contributionsTitle);
            
            // Create responsive container
            const tableResponsive = document.createElement('div');
            tableResponsive.className = 'table-responsive';
            
            // Create table
            const table = document.createElement('table');
            table.className = 'synergy-contributions table';
            
            // Create header row
            const headerRow = document.createElement('tr');
            headerRow.innerHTML = '<th>Component</th><th>Contribution</th><th>Concentration</th>';
            table.appendChild(headerRow);
            
            // Create rows for each component
            contributions.forEach(component => {
                const row = document.createElement('tr');
                
                // Component name
                const cell1 = document.createElement('td');
                cell1.textContent = component.name;
                row.appendChild(cell1);
                
                // Contribution
                const cell2 = document.createElement('td');
                cell2.textContent = component.contribution.toFixed(3);
                row.appendChild(cell2);
                
                // Concentration
                const cell3 = document.createElement('td');
                cell3.textContent = `${component.concentration.toFixed(2)} ${component.concentration_unit}`;
                row.appendChild(cell3);
                
                table.appendChild(row);
            });
            
            // Add table to responsive container
            tableResponsive.appendChild(table);
            
            // Add responsive container to main container
            container.appendChild(tableResponsive);
        }
    },

    /**
     * Render mixture analysis and recommendations.
     * 
     * @param {Object} analysisData - Mixture analysis data
     * @param {string} containerId - ID of the container element
     */
    renderAnalysis: function(analysisData, containerId) {
        const container = document.getElementById(containerId);
        if (!container) {
            console.error(`Container element with ID ${containerId} not found`);
            return;
        }

        // Clear container
        container.innerHTML = '';

        // Create analysis summary
        const summary = document.createElement('div');
        summary.className = 'analysis-summary';
        
        // Mixture name
        const mixtureName = document.createElement('h2');
        mixtureName.textContent = analysisData.mixture_name;
        summary.appendChild(mixtureName);
        
        // Overall score
        const overallScore = document.createElement('h3');
        overallScore.textContent = `Overall Score: ${analysisData.properties["Cryoprotection Score"].toFixed(2)}`;
        summary.appendChild(overallScore);
        
        // Compatibility status
        const compatibility = document.createElement('p');
        compatibility.textContent = `Compatibility: ${analysisData.compatibility.compatibility_status} (${analysisData.compatibility.overall_compatibility_score.toFixed(2)})`;
        summary.appendChild(compatibility);
        
        // Synergy type
        const synergy = document.createElement('p');
        synergy.textContent = `Synergy: ${analysisData.synergy.synergy_type} (${analysisData.synergy.synergy_score.toFixed(2)})`;
        summary.appendChild(synergy);
        
        container.appendChild(summary);

        // Strengths and weaknesses
        const strengthsWeaknesses = document.createElement('div');
        strengthsWeaknesses.className = 'strengths-weaknesses';
        
        // Strengths
        const strengthsTitle = document.createElement('h3');
        strengthsTitle.textContent = 'Strengths';
        strengthsWeaknesses.appendChild(strengthsTitle);
        
        const strengthsList = document.createElement('ul');
        if (analysisData.strengths && analysisData.strengths.length > 0) {
            analysisData.strengths.forEach(strength => {
                const li = document.createElement('li');
                li.textContent = strength;
                strengthsList.appendChild(li);
            });
        } else {
            const li = document.createElement('li');
            li.textContent = 'No strengths identified';
            strengthsList.appendChild(li);
        }
        strengthsWeaknesses.appendChild(strengthsList);
        
        // Weaknesses
        const weaknessesTitle = document.createElement('h3');
        weaknessesTitle.textContent = 'Weaknesses';
        strengthsWeaknesses.appendChild(weaknessesTitle);
        
        const weaknessesList = document.createElement('ul');
        if (analysisData.weaknesses && analysisData.weaknesses.length > 0) {
            analysisData.weaknesses.forEach(weakness => {
                const li = document.createElement('li');
                li.textContent = weakness;
                weaknessesList.appendChild(li);
            });
        } else {
            const li = document.createElement('li');
            li.textContent = 'No weaknesses identified';
            weaknessesList.appendChild(li);
        }
        strengthsWeaknesses.appendChild(weaknessesList);
        
        container.appendChild(strengthsWeaknesses);

        // Recommendations
        const recommendationsDiv = document.createElement('div');
        recommendationsDiv.className = 'recommendations';
        
        const recommendationsTitle = document.createElement('h3');
        recommendationsTitle.textContent = 'Recommendations';
        recommendationsDiv.appendChild(recommendationsTitle);
        
        if (analysisData.recommendations && analysisData.recommendations.length > 0) {
            const recommendationsList = document.createElement('ul');
            
            analysisData.recommendations.forEach(recommendation => {
                const li = document.createElement('li');
                
                const type = document.createElement('strong');
                type.textContent = recommendation.type + ': ';
                li.appendChild(type);
                
                li.appendChild(document.createTextNode(recommendation.description));
                
                const action = document.createElement('div');
                action.className = 'recommendation-action';
                action.textContent = recommendation.action;
                li.appendChild(action);
                
                recommendationsList.appendChild(li);
            });
            
            recommendationsDiv.appendChild(recommendationsList);
        } else {
            const noRecommendations = document.createElement('p');
            noRecommendations.textContent = 'No recommendations available';
            recommendationsDiv.appendChild(noRecommendations);
        }
        
        container.appendChild(recommendationsDiv);
    }
};

// Export the module
if (typeof module !== 'undefined' && module.exports) {
    module.exports = MixtureAnalysis;
}
/**
 * CryoProtect Analyzer - RDKit Integration
 * 
 * This module provides functions for interacting with the RDKit API endpoints
 * and displaying molecular structures and properties in the web interface.
 */

// API endpoint base URL
const API_BASE_URL = '/api/v1';

/**
 * Calculate molecular properties using the RDKit API.
 * 
 * @param {string} moleculeData - Molecular data (SMILES, MOL, SDF)
 * @param {string} inputFormat - Input format ('smiles', 'mol', 'sdf')
 * @returns {Promise<Object>} - Promise resolving to property data
 */
async function calculateProperties(moleculeData, inputFormat = 'smiles') {
    try {
        const response = await fetch(`${API_BASE_URL}/rdkit/properties`, {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
                'Authorization': `Bearer ${getAuthToken()}`
            },
            body: JSON.stringify({
                molecule_data: moleculeData,
                input_format: inputFormat
            })
        });

        if (!response.ok) {
            const errorData = await response.json();
            throw new Error(errorData.message || 'Failed to calculate properties');
        }

        return await response.json();
    } catch (error) {
        console.error('Error calculating properties:', error);
        throw error;
    }
}

/**
 * Generate a molecular visualization using the RDKit API.
 * 
 * @param {string} moleculeData - Molecular data (SMILES, MOL, SDF)
 * @param {string} inputFormat - Input format ('smiles', 'mol', 'sdf')
 * @param {number} width - Image width in pixels
 * @param {number} height - Image height in pixels
 * @param {Array<number>} highlightAtoms - Optional array of atom indices to highlight
 * @returns {Promise<Object>} - Promise resolving to visualization data
 */
async function generateVisualization(moleculeData, inputFormat = 'smiles', width = 400, height = 300, highlightAtoms = null) {
    try {
        const requestData = {
            molecule_data: moleculeData,
            input_format: inputFormat,
            width: width,
            height: height
        };

        if (highlightAtoms) {
            requestData.highlight_atoms = highlightAtoms;
        }

        const response = await fetch(`${API_BASE_URL}/rdkit/visualization`, {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
                'Authorization': `Bearer ${getAuthToken()}`
            },
            body: JSON.stringify(requestData)
        });

        if (!response.ok) {
            const errorData = await response.json();
            throw new Error(errorData.message || 'Failed to generate visualization');
        }

        return await response.json();
    } catch (error) {
        console.error('Error generating visualization:', error);
        throw error;
    }
}

/**
 * Perform a substructure search using the RDKit API.
 * 
 * @param {string} queryMolData - Query molecule data (SMARTS pattern)
 * @param {string} targetMolData - Target molecule data to search in
 * @param {string} queryFormat - Format of the query data ('smarts', 'smiles')
 * @param {string} targetFormat - Format of the target data ('smiles', 'mol', 'sdf')
 * @returns {Promise<Object>} - Promise resolving to search results
 */
async function performSubstructureSearch(queryMolData, targetMolData, queryFormat = 'smarts', targetFormat = 'smiles') {
    try {
        const response = await fetch(`${API_BASE_URL}/rdkit/substructure`, {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
                'Authorization': `Bearer ${getAuthToken()}`
            },
            body: JSON.stringify({
                query_mol_data: queryMolData,
                target_mol_data: targetMolData,
                query_format: queryFormat,
                target_format: targetFormat
            })
        });

        if (!response.ok) {
            const errorData = await response.json();
            throw new Error(errorData.message || 'Failed to perform substructure search');
        }

        return await response.json();
    } catch (error) {
        console.error('Error performing substructure search:', error);
        throw error;
    }
}

/**
 * Calculate molecular similarity using the RDKit API.
 * 
 * @param {string} mol1Data - First molecule data
 * @param {string} mol2Data - Second molecule data
 * @param {string} mol1Format - Format of the first molecule data
 * @param {string} mol2Format - Format of the second molecule data
 * @param {string} fingerprintType - Type of fingerprint to use ('morgan', 'maccs', 'topological')
 * @returns {Promise<Object>} - Promise resolving to similarity data
 */
async function calculateSimilarity(mol1Data, mol2Data, mol1Format = 'smiles', mol2Format = 'smiles', fingerprintType = 'morgan') {
    try {
        const response = await fetch(`${API_BASE_URL}/rdkit/similarity`, {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
                'Authorization': `Bearer ${getAuthToken()}`
            },
            body: JSON.stringify({
                mol1_data: mol1Data,
                mol2_data: mol2Data,
                mol1_format: mol1Format,
                mol2_format: mol2Format,
                fingerprint_type: fingerprintType
            })
        });

        if (!response.ok) {
            const errorData = await response.json();
            throw new Error(errorData.message || 'Failed to calculate similarity');
        }

        return await response.json();
    } catch (error) {
        console.error('Error calculating similarity:', error);
        throw error;
    }
}

/**
 * Calculate properties for a molecule in the database.
 * 
 * @param {string} moleculeId - ID of the molecule in the database
 * @returns {Promise<Object>} - Promise resolving to result data
 */
async function calculateMoleculeProperties(moleculeId) {
    try {
        const response = await fetch(`${API_BASE_URL}/molecules/${moleculeId}/calculate-properties`, {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
                'Authorization': `Bearer ${getAuthToken()}`
            }
        });

        if (!response.ok) {
            const errorData = await response.json();
            throw new Error(errorData.message || 'Failed to calculate molecule properties');
        }

        return await response.json();
    } catch (error) {
        console.error('Error calculating molecule properties:', error);
        throw error;
    }
}

/**
 * Display molecular properties in a container element.
 * 
 * @param {Object} properties - Property data from the API
 * @param {HTMLElement} container - Container element to display properties in
 */
function displayMolecularProperties(properties, container) {
    // Clear the container
    container.innerHTML = '';
    
    // Create property sections
    const sections = [
        { title: 'Basic Properties', properties: ['molecular_weight', 'exact_mass', 'logp', 'tpsa'] },
        { title: 'Hydrogen Bonding', properties: ['hydrogen_bonding'] },
        { title: 'Molecular Properties', properties: ['molecular_properties'] },
        { title: 'Functional Groups', properties: ['functional_groups'] },
        { title: 'Permeability', properties: ['permeability'] }
    ];
    
    // Create a card for each section
    sections.forEach(section => {
        const sectionDiv = document.createElement('div');
        sectionDiv.className = 'card mb-3';
        
        const cardHeader = document.createElement('div');
        cardHeader.className = 'card-header';
        cardHeader.textContent = section.title;
        sectionDiv.appendChild(cardHeader);
        
        const cardBody = document.createElement('div');
        cardBody.className = 'card-body';
        
        // Add properties for this section
        section.properties.forEach(propKey => {
            if (propKey in properties) {
                const propValue = properties[propKey];
                
                if (typeof propValue === 'object') {
                    // For nested objects (like hydrogen_bonding)
                    const subTable = document.createElement('table');
                    subTable.className = 'table table-sm';
                    
                    const tbody = document.createElement('tbody');
                    
                    Object.entries(propValue).forEach(([key, value]) => {
                        const row = document.createElement('tr');
                        
                        const keyCell = document.createElement('th');
                        keyCell.textContent = key.replace(/_/g, ' ').replace(/\b\w/g, l => l.toUpperCase());
                        keyCell.style.width = '50%';
                        row.appendChild(keyCell);
                        
                        const valueCell = document.createElement('td');
                        valueCell.textContent = value;
                        row.appendChild(valueCell);
                        
                        tbody.appendChild(row);
                    });
                    
                    subTable.appendChild(tbody);
                    cardBody.appendChild(subTable);
                } else {
                    // For simple properties
                    const propRow = document.createElement('div');
                    propRow.className = 'row mb-2';
                    
                    const propNameCol = document.createElement('div');
                    propNameCol.className = 'col-6';
                    propNameCol.textContent = propKey.replace(/_/g, ' ').replace(/\b\w/g, l => l.toUpperCase()) + ':';
                    propRow.appendChild(propNameCol);
                    
                    const propValueCol = document.createElement('div');
                    propValueCol.className = 'col-6';
                    propValueCol.textContent = propValue;
                    propRow.appendChild(propValueCol);
                    
                    cardBody.appendChild(propRow);
                }
            }
        });
        
        sectionDiv.appendChild(cardBody);
        container.appendChild(sectionDiv);
    });
}

/**
 * Display a molecular visualization in a container element.
 * 
 * @param {string} svgData - SVG data from the API
 * @param {HTMLElement} container - Container element to display visualization in
 */
function displayMolecularVisualization(svgData, container) {
    // Clear the container
    container.innerHTML = '';
    
    // Set the SVG data
    container.innerHTML = svgData;
    
    // Make sure SVG is responsive
    const svg = container.querySelector('svg');
    if (svg) {
        svg.setAttribute('width', '100%');
        svg.setAttribute('height', 'auto');
        svg.style.maxWidth = '100%';
    }
}

/**
 * Initialize the molecule search form.
 * 
 * @param {string} formId - ID of the search form
 * @param {string} resultsContainerId - ID of the container to display results in
 */
function initMoleculeSearchForm(formId, resultsContainerId) {
    const form = document.getElementById(formId);
    const resultsContainer = document.getElementById(resultsContainerId);
    
    if (!form || !resultsContainer) {
        console.error('Form or results container not found');
        return;
    }
    
    form.addEventListener('submit', async (event) => {
        event.preventDefault();
        
        // Show loading indicator
        resultsContainer.innerHTML = '<div class="text-center"><div class="spinner-border" role="status"><span class="visually-hidden">Loading...</span></div></div>';
        
        // Get form data
        const formData = new FormData(form);
        const searchType = formData.get('searchType');
        const query = formData.get('query');
        
        try {
            let results;
            
            if (searchType === 'substructure') {
                // Perform substructure search
                results = await performSubstructureSearch(query, formData.get('targetMolecule'));
            } else if (searchType === 'similarity') {
                // Calculate similarity
                results = await calculateSimilarity(query, formData.get('targetMolecule'));
            } else {
                throw new Error('Invalid search type');
            }
            
            // Display results
            displaySearchResults(results, searchType, resultsContainer);
        } catch (error) {
            resultsContainer.innerHTML = `<div class="alert alert-danger">${error.message}</div>`;
        }
    });
}

/**
 * Display search results in a container element.
 * 
 * @param {Object} results - Search results from the API
 * @param {string} searchType - Type of search ('substructure' or 'similarity')
 * @param {HTMLElement} container - Container element to display results in
 */
function displaySearchResults(results, searchType, container) {
    // Clear the container
    container.innerHTML = '';
    
    // Create a card for the results
    const card = document.createElement('div');
    card.className = 'card';
    
    const cardHeader = document.createElement('div');
    cardHeader.className = 'card-header';
    cardHeader.textContent = searchType === 'substructure' ? 'Substructure Search Results' : 'Similarity Search Results';
    card.appendChild(cardHeader);
    
    const cardBody = document.createElement('div');
    cardBody.className = 'card-body';
    
    if (searchType === 'substructure') {
        // Display substructure search results
        const matchText = document.createElement('p');
        matchText.className = results.match ? 'text-success' : 'text-danger';
        matchText.textContent = results.match ? `Match found! (${results.match_count} matches)` : 'No matches found.';
        cardBody.appendChild(matchText);
        
        if (results.visualization) {
            const visualizationDiv = document.createElement('div');
            visualizationDiv.className = 'mt-3';
            visualizationDiv.innerHTML = results.visualization;
            cardBody.appendChild(visualizationDiv);
        }
    } else if (searchType === 'similarity') {
        // Display similarity search results
        const tanimotoRow = document.createElement('div');
        tanimotoRow.className = 'row mb-2';
        
        const tanimotoLabelCol = document.createElement('div');
        tanimotoLabelCol.className = 'col-6';
        tanimotoLabelCol.textContent = 'Tanimoto Similarity:';
        tanimotoRow.appendChild(tanimotoLabelCol);
        
        const tanimotoValueCol = document.createElement('div');
        tanimotoValueCol.className = 'col-6';
        tanimotoValueCol.textContent = results.tanimoto.toFixed(4);
        tanimotoRow.appendChild(tanimotoValueCol);
        
        cardBody.appendChild(tanimotoRow);
        
        const diceRow = document.createElement('div');
        diceRow.className = 'row mb-2';
        
        const diceLabelCol = document.createElement('div');
        diceLabelCol.className = 'col-6';
        diceLabelCol.textContent = 'Dice Similarity:';
        diceRow.appendChild(diceLabelCol);
        
        const diceValueCol = document.createElement('div');
        diceValueCol.className = 'col-6';
        diceValueCol.textContent = results.dice.toFixed(4);
        diceRow.appendChild(diceValueCol);
        
        cardBody.appendChild(diceRow);
        
        const fingerprintRow = document.createElement('div');
        fingerprintRow.className = 'row mb-2';
        
        const fingerprintLabelCol = document.createElement('div');
        fingerprintLabelCol.className = 'col-6';
        fingerprintLabelCol.textContent = 'Fingerprint Type:';
        fingerprintRow.appendChild(fingerprintLabelCol);
        
        const fingerprintValueCol = document.createElement('div');
        fingerprintValueCol.className = 'col-6';
        fingerprintValueCol.textContent = results.fingerprint_type;
        fingerprintRow.appendChild(fingerprintValueCol);
        
        cardBody.appendChild(fingerprintRow);
    }
    
    card.appendChild(cardBody);
    container.appendChild(card);
}

/**
 * Get the authentication token from local storage.
 * 
 * @returns {string} - Authentication token or empty string if not found
 */
function getAuthToken() {
    return localStorage.getItem('authToken') || '';
}

// Export functions for use in other modules
window.RDKitAPI = {
    calculateProperties,
    generateVisualization,
    performSubstructureSearch,
    calculateSimilarity,
    calculateMoleculeProperties,
    displayMolecularProperties,
    displayMolecularVisualization,
    initMoleculeSearchForm
};
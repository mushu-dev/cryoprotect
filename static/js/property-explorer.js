/**
 * Property Explorer JS
 * 
 * This script handles the functionality for the Property Explorer page,
 * including data visualization, filtering, and comparison tools.
 */

// Initialize visualizations and data when the DOM is loaded
document.addEventListener('DOMContentLoaded', function() {
    // Load dashboard stats
    loadDashboardStats();
    
    // Initialize property type filters
    initializePropertyFilters();
    
    // Initialize visualizations
    initializeVisualizations();
    
    // Initialize molecule grid
    initializeMoleculeGrid();
    
    // Initialize comparison tool
    initializeComparisonTool();
    
    // Initialize correlation analysis
    initializeCorrelationAnalysis();
    
    // Initialize export tool
    initializeExportTool();
    
    // Set up event listeners
    setupEventListeners();
});

/**
 * Load dashboard statistics
 */
async function loadDashboardStats() {
    try {
        // In a real implementation, this would fetch data from the API
        const response = await fetch('/api/v1/stats/dashboard');
        const data = await response.json();
        
        // For demo purposes, use placeholder data
        const stats = {
            totalMolecules: 5247,
            totalProperties: 32,
            totalCryoprotectants: 158,
            dataCoverage: '87%'
        };
        
        // Update dashboard cards
        document.getElementById('totalMolecules').textContent = stats.totalMolecules.toLocaleString();
        document.getElementById('totalProperties').textContent = stats.totalProperties.toLocaleString();
        document.getElementById('totalCryoprotectants').textContent = stats.totalCryoprotectants.toLocaleString();
        document.getElementById('dataCoverage').textContent = stats.dataCoverage;
    } catch (error) {
        console.error('Error loading dashboard stats:', error);
    }
}

/**
 * Initialize property filters
 */
function initializePropertyFilters() {
    // In a real implementation, this would fetch property types from the API
    // For demo purposes, use placeholder data
    const propertyTypes = [
        { id: 'physical', name: 'Physical Properties', checked: true },
        { id: 'topological', name: 'Topological Properties', checked: true },
        { id: 'cryoprotectant', name: 'Cryoprotectant Properties', checked: true },
        { id: 'descriptors', name: 'Molecular Descriptors', checked: false },
        { id: 'toxicity', name: 'Toxicity Properties', checked: false }
    ];
    
    // Generate property type checkboxes
    const propertyTypeFilters = document.getElementById('propertyTypeFilters');
    propertyTypeFilters.innerHTML = '';
    
    propertyTypes.forEach(type => {
        const checkboxDiv = document.createElement('div');
        checkboxDiv.className = 'form-check';
        
        const checkbox = document.createElement('input');
        checkbox.className = 'form-check-input';
        checkbox.type = 'checkbox';
        checkbox.id = type.id + 'Properties';
        checkbox.checked = type.checked;
        
        const label = document.createElement('label');
        label.className = 'form-check-label';
        label.htmlFor = type.id + 'Properties';
        label.textContent = type.name;
        
        checkboxDiv.appendChild(checkbox);
        checkboxDiv.appendChild(label);
        propertyTypeFilters.appendChild(checkboxDiv);
    });
    
    // Generate property toggles
    const properties = [
        { id: 'mw', name: 'Molecular Weight', checked: true, type: 'physical' },
        { id: 'logp', name: 'LogP', checked: true, type: 'physical' },
        { id: 'tpsa', name: 'TPSA', checked: false, type: 'topological' },
        { id: 'hba', name: 'H-Bond Acceptors', checked: false, type: 'physical' },
        { id: 'hbd', name: 'H-Bond Donors', checked: false, type: 'physical' },
        { id: 'rotb', name: 'Rotatable Bonds', checked: false, type: 'topological' },
        { id: 'fsp3', name: 'Fraction sp3', checked: false, type: 'topological' },
        { id: 'cypc', name: 'Cryoprotection Score', checked: true, type: 'cryoprotectant' },
        { id: 'aroRings', name: 'Aromatic Rings', checked: false, type: 'topological' },
        { id: 'heavyAtoms', name: 'Heavy Atoms', checked: false, type: 'physical' },
        { id: 'qed', name: 'QED Score', checked: false, type: 'descriptors' },
        { id: 'toxicity', name: 'Toxicity Score', checked: false, type: 'toxicity' }
    ];
    
    const propertyToggleContainer = document.getElementById('propertyToggleContainer');
    propertyToggleContainer.innerHTML = '';
    
    properties.forEach(property => {
        const colDiv = document.createElement('div');
        colDiv.className = 'col-md-3 mb-2';
        colDiv.dataset.propertyType = property.type;
        
        const checkboxDiv = document.createElement('div');
        checkboxDiv.className = 'form-check property-toggle';
        
        const checkbox = document.createElement('input');
        checkbox.className = 'form-check-input';
        checkbox.type = 'checkbox';
        checkbox.id = property.id;
        checkbox.checked = property.checked;
        
        const label = document.createElement('label');
        label.className = 'form-check-label';
        label.htmlFor = property.id;
        label.textContent = property.name;
        
        checkboxDiv.appendChild(checkbox);
        checkboxDiv.appendChild(label);
        colDiv.appendChild(checkboxDiv);
        propertyToggleContainer.appendChild(colDiv);
    });
    
    // Filter properties based on property type selection
    const propertyTypeCheckboxes = document.querySelectorAll('#propertyTypeFilters input[type="checkbox"]');
    propertyTypeCheckboxes.forEach(checkbox => {
        checkbox.addEventListener('change', function() {
            const typeId = this.id.replace('Properties', '');
            const propertyElements = document.querySelectorAll(`#propertyToggleContainer div[data-property-type="${typeId}"]`);
            
            propertyElements.forEach(element => {
                element.style.display = this.checked ? 'block' : 'none';
            });
        });
    });
}

/**
 * Initialize visualizations
 */
function initializeVisualizations() {
    // Get visualization container
    const visualizationContainer = document.getElementById('propertyVisualization');
    
    // Create a canvas for Chart.js
    const canvas = document.createElement('canvas');
    canvas.id = 'visualizationCanvas';
    visualizationContainer.innerHTML = '';
    visualizationContainer.appendChild(canvas);
    
    // Create a histogram with placeholder data
    const ctx = canvas.getContext('2d');
    const chart = new Chart(ctx, {
        type: 'bar',
        data: {
            labels: ['0-100', '100-200', '200-300', '300-400', '400-500', '500-600', '600-700', '700-800', '800-900', '900-1000'],
            datasets: [{
                label: 'Molecular Weight Distribution',
                data: [42, 128, 267, 396, 512, 433, 287, 165, 83, 27],
                backgroundColor: 'rgba(54, 162, 235, 0.5)',
                borderColor: 'rgba(54, 162, 235, 1)',
                borderWidth: 1
            }]
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            scales: {
                y: {
                    beginAtZero: true,
                    title: {
                        display: true,
                        text: 'Count'
                    }
                },
                x: {
                    title: {
                        display: true,
                        text: 'Molecular Weight (g/mol)'
                    }
                }
            },
            plugins: {
                title: {
                    display: true,
                    text: 'Molecular Weight Distribution',
                    font: {
                        size: 16
                    }
                },
                legend: {
                    display: false
                }
            }
        }
    });
    
    // Store chart instance for later updates
    window.currentChart = chart;
    
    // Update visualization when visualization type changes
    document.getElementById('visualizationType').addEventListener('change', function() {
        updateVisualization(this.value);
    });
}

/**
 * Update visualization based on type
 * @param {string} type - Visualization type
 */
function updateVisualization(type) {
    // Destroy existing chart
    if (window.currentChart) {
        window.currentChart.destroy();
    }
    
    const canvas = document.getElementById('visualizationCanvas');
    const ctx = canvas.getContext('2d');
    
    // Create new chart based on selected type
    switch (type) {
        case 'histogram':
            window.currentChart = new Chart(ctx, {
                type: 'bar',
                data: {
                    labels: ['0-100', '100-200', '200-300', '300-400', '400-500', '500-600', '600-700', '700-800', '800-900', '900-1000'],
                    datasets: [{
                        label: 'Molecular Weight Distribution',
                        data: [42, 128, 267, 396, 512, 433, 287, 165, 83, 27],
                        backgroundColor: 'rgba(54, 162, 235, 0.5)',
                        borderColor: 'rgba(54, 162, 235, 1)',
                        borderWidth: 1
                    }]
                },
                options: {
                    responsive: true,
                    maintainAspectRatio: false,
                    scales: {
                        y: {
                            beginAtZero: true,
                            title: {
                                display: true,
                                text: 'Count'
                            }
                        },
                        x: {
                            title: {
                                display: true,
                                text: 'Molecular Weight (g/mol)'
                            }
                        }
                    },
                    plugins: {
                        title: {
                            display: true,
                            text: 'Molecular Weight Distribution',
                            font: {
                                size: 16
                            }
                        },
                        legend: {
                            display: false
                        }
                    }
                }
            });
            break;
            
        case 'scatter':
            window.currentChart = new Chart(ctx, {
                type: 'scatter',
                data: {
                    datasets: [{
                        label: 'LogP vs. Molecular Weight',
                        data: [
                            { x: 100, y: -0.5 },
                            { x: 150, y: 0.2 },
                            { x: 200, y: 1.5 },
                            { x: 250, y: 1.8 },
                            { x: 300, y: 2.3 },
                            { x: 350, y: 2.7 },
                            { x: 400, y: 3.1 },
                            { x: 450, y: 3.5 },
                            { x: 500, y: 3.8 },
                            { x: 550, y: 4.2 },
                            { x: 600, y: 4.5 },
                            { x: 650, y: 4.8 },
                            { x: 700, y: 5.1 }
                        ],
                        backgroundColor: 'rgba(75, 192, 192, 0.7)',
                        borderColor: 'rgba(75, 192, 192, 1)',
                        borderWidth: 1,
                        pointRadius: 5,
                        pointHoverRadius: 8
                    }]
                },
                options: {
                    responsive: true,
                    maintainAspectRatio: false,
                    scales: {
                        y: {
                            title: {
                                display: true,
                                text: 'LogP'
                            }
                        },
                        x: {
                            title: {
                                display: true,
                                text: 'Molecular Weight (g/mol)'
                            }
                        }
                    },
                    plugins: {
                        title: {
                            display: true,
                            text: 'LogP vs. Molecular Weight',
                            font: {
                                size: 16
                            }
                        }
                    }
                }
            });
            break;
            
        case 'boxplot':
            // For boxplot, we'd normally use a plugin like chartjs-chart-box-and-violin-plot
            // For this demo, use a bar chart with error bars to simulate boxplots
            window.currentChart = new Chart(ctx, {
                type: 'bar',
                data: {
                    labels: ['Molecular Weight', 'LogP', 'TPSA', 'Cryoprotection Score'],
                    datasets: [{
                        label: 'Property Distribution',
                        data: [400, 3.2, 75, 0.65],
                        backgroundColor: [
                            'rgba(54, 162, 235, 0.5)',
                            'rgba(75, 192, 192, 0.5)',
                            'rgba(255, 206, 86, 0.5)',
                            'rgba(153, 102, 255, 0.5)'
                        ],
                        borderColor: [
                            'rgba(54, 162, 235, 1)',
                            'rgba(75, 192, 192, 1)',
                            'rgba(255, 206, 86, a)',
                            'rgba(153, 102, 255, 1)'
                        ],
                        borderWidth: 1,
                        // Error bars to simulate min/max
                        errorBars: {
                            'Molecular Weight': { plus: 250, minus: 250 },
                            'LogP': { plus: 2, minus: 2 },
                            'TPSA': { plus: 35, minus: 35 },
                            'Cryoprotection Score': { plus: 0.25, minus: 0.35 }
                        }
                    }]
                },
                options: {
                    responsive: true,
                    maintainAspectRatio: false,
                    scales: {
                        y: {
                            beginAtZero: true,
                            title: {
                                display: true,
                                text: 'Value (normalized)'
                            }
                        }
                    },
                    plugins: {
                        title: {
                            display: true,
                            text: 'Property Distributions',
                            font: {
                                size: 16
                            }
                        }
                    }
                }
            });
            break;
            
        case 'radar':
            window.currentChart = new Chart(ctx, {
                type: 'radar',
                data: {
                    labels: ['MW', 'LogP', 'TPSA', 'HBA', 'HBD', 'RotB', 'CypC'],
                    datasets: [{
                        label: 'Average Molecule',
                        data: [0.5, 0.7, 0.4, 0.6, 0.3, 0.5, 0.6],
                        backgroundColor: 'rgba(54, 162, 235, 0.3)',
                        borderColor: 'rgb(54, 162, 235)',
                        pointBackgroundColor: 'rgb(54, 162, 235)',
                        pointHoverBackgroundColor: '#fff',
                        pointHoverBorderColor: 'rgb(54, 162, 235)'
                    }, {
                        label: 'Ideal Cryoprotectant',
                        data: [0.6, 0.8, 0.5, 0.7, 0.9, 0.3, 0.9],
                        backgroundColor: 'rgba(255, 99, 132, 0.3)',
                        borderColor: 'rgb(255, 99, 132)',
                        pointBackgroundColor: 'rgb(255, 99, 132)',
                        pointHoverBackgroundColor: '#fff',
                        pointHoverBorderColor: 'rgb(255, 99, 132)'
                    }]
                },
                options: {
                    responsive: true,
                    maintainAspectRatio: false,
                    elements: {
                        line: {
                            borderWidth: 2
                        }
                    },
                    plugins: {
                        title: {
                            display: true,
                            text: 'Property Radar Chart',
                            font: {
                                size: 16
                            }
                        }
                    },
                    scales: {
                        r: {
                            angleLines: {
                                display: true
                            },
                            suggestedMin: 0,
                            suggestedMax: 1
                        }
                    }
                }
            });
            break;
    }
}

/**
 * Initialize molecule grid
 */
function initializeMoleculeGrid() {
    // Prepare placeholder data
    const data = [
        { id: 1, name: 'Ethanol', formula: 'C2H6O', mw: 46.07, logp: -0.31, tpsa: 20.23, cypc: 0.78, smiles: 'CCO' },
        { id: 2, name: 'Glycerol', formula: 'C3H8O3', mw: 92.09, logp: -1.76, tpsa: 60.69, cypc: 0.92, smiles: 'C(C(CO)O)O' },
        { id: 3, name: 'Dimethyl Sulfoxide', formula: 'C2H6OS', mw: 78.13, logp: -1.35, tpsa: 17.07, cypc: 0.86, smiles: 'CS(=O)C' },
        { id: 4, name: 'Propylene Glycol', formula: 'C3H8O2', mw: 76.10, logp: -0.92, tpsa: 40.46, cypc: 0.81, smiles: 'CC(CO)O' },
        { id: 5, name: 'Trehalose', formula: 'C12H22O11', mw: 342.30, logp: -4.23, tpsa: 189.53, cypc: 0.89, smiles: 'C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@]2([C@H]([C@@H]([C@H](O2)CO)O)O)CO)O)O)O)O' }
    ];
    
    // Initialize GridJS
    const grid = new gridjs.Grid({
        columns: [
            { id: 'id', name: 'ID', width: '60px' },
            { 
                id: 'name', 
                name: 'Name',
                formatter: (cell, row) => {
                    return gridjs.html(`<a href="#" class="molecule-link" data-id="${row.cells[0].data}">${cell}</a>`);
                }
            },
            { id: 'formula', name: 'Formula' },
            { 
                id: 'mw', 
                name: 'MW (g/mol)',
                formatter: (cell) => {
                    return cell.toFixed(2);
                }
            },
            { 
                id: 'logp', 
                name: 'LogP',
                formatter: (cell) => {
                    return cell.toFixed(2);
                }
            },
            { 
                id: 'cypc', 
                name: 'CypC',
                formatter: (cell) => {
                    // Color-code cryoprotection score
                    let color = 'green';
                    if (cell < 0.5) color = 'red';
                    else if (cell < 0.7) color = 'orange';
                    
                    return gridjs.html(`<span style="color: ${color}; font-weight: bold;">${cell.toFixed(2)}</span>`);
                }
            },
            {
                id: 'actions',
                name: 'Actions',
                formatter: (_, row) => {
                    return gridjs.html(`
                        <div class="btn-group btn-group-sm" role="group">
                            <button type="button" class="btn btn-outline-primary view-btn" data-id="${row.cells[0].data}">
                                <i class="bi bi-eye"></i>
                            </button>
                            <button type="button" class="btn btn-outline-success compare-btn" data-id="${row.cells[0].data}">
                                <i class="bi bi-bar-chart"></i>
                            </button>
                        </div>
                    `);
                }
            }
        ],
        data: data,
        search: true,
        pagination: {
            limit: 10
        },
        sort: true,
        className: {
            table: 'table table-striped table-hover'
        }
    }).render(document.getElementById('moleculeGrid'));
    
    // Update result count
    document.getElementById('resultCount').textContent = data.length;
    
    // Add event listeners to view and compare buttons
    grid.on('rendered', () => {
        // View buttons
        document.querySelectorAll('.view-btn').forEach(button => {
            button.addEventListener('click', function() {
                const id = this.dataset.id;
                // In a real implementation, this would open a molecule details view
                alert(`View molecule ${id}`);
            });
        });
        
        // Compare buttons
        document.querySelectorAll('.compare-btn').forEach(button => {
            button.addEventListener('click', function() {
                const id = this.dataset.id;
                // In a real implementation, this would add the molecule to the comparison
                alert(`Add molecule ${id} to comparison`);
            });
        });
        
        // Molecule links
        document.querySelectorAll('.molecule-link').forEach(link => {
            link.addEventListener('click', function(e) {
                e.preventDefault();
                const id = this.dataset.id;
                // In a real implementation, this would open a molecule details view
                alert(`View molecule ${id}`);
            });
        });
    });
}

/**
 * Initialize comparison tool
 */
function initializeComparisonTool() {
    // Load molecules for comparison selector
    const molecules = [
        { id: 1, name: 'Ethanol' },
        { id: 2, name: 'Glycerol' },
        { id: 3, name: 'Dimethyl Sulfoxide' },
        { id: 4, name: 'Propylene Glycol' },
        { id: 5, name: 'Trehalose' },
        { id: 6, name: 'Sucrose' },
        { id: 7, name: 'Glucose' },
        { id: 8, name: 'Methanol' },
        { id: 9, name: 'Formamide' },
        { id: 10, name: 'Acetamide' }
    ];
    
    // Populate comparison molecule select
    const comparisonMoleculeSelect = document.getElementById('comparisonMoleculeSelect');
    
    // Clear existing options except the first one
    while (comparisonMoleculeSelect.options.length > 1) {
        comparisonMoleculeSelect.remove(1);
    }
    
    // Add molecule options
    molecules.forEach(molecule => {
        const option = document.createElement('option');
        option.value = molecule.id;
        option.textContent = molecule.name;
        comparisonMoleculeSelect.appendChild(option);
    });
    
    // Add to comparison button
    document.getElementById('addToComparisonBtn').addEventListener('click', function() {
        const moleculeId = comparisonMoleculeSelect.value;
        if (moleculeId) {
            addMoleculeToComparison(moleculeId);
        }
    });
    
    // Generate comparison button
    document.getElementById('generateComparisonBtn').addEventListener('click', function() {
        generateComparison();
    });
    
    // Clear comparison button
    document.getElementById('clearComparisonBtn').addEventListener('click', function() {
        clearComparison();
    });
}

/**
 * Add molecule to comparison
 * @param {string} moleculeId - Molecule ID
 */
function addMoleculeToComparison(moleculeId) {
    // Get molecule from selection
    const select = document.getElementById('comparisonMoleculeSelect');
    const moleculeName = select.options[select.selectedIndex].text;
    
    // Find an empty comparison target
    const emptyTargets = document.querySelectorAll('.comparison-target.empty');
    if (emptyTargets.length === 0) {
        alert('All comparison slots are filled. Please clear one first.');
        return;
    }
    
    // Fill the first empty target
    const target = emptyTargets[0];
    target.innerHTML = `
        <div class="text-center">
            <div class="fw-bold">${moleculeName}</div>
            <div class="text-muted">ID: ${moleculeId}</div>
            <button class="btn btn-sm btn-outline-danger mt-2 remove-comparison-btn">
                <i class="bi bi-x-circle"></i> Remove
            </button>
        </div>
    `;
    target.classList.remove('empty');
    target.dataset.moleculeId = moleculeId;
    
    // Add remove button event listener
    target.querySelector('.remove-comparison-btn').addEventListener('click', function() {
        target.innerHTML = `
            <div class="text-muted">
                <i class="bi bi-plus-circle"></i>
                <div>Drop molecule here</div>
            </div>
        `;
        target.classList.add('empty');
        delete target.dataset.moleculeId;
    });
}

/**
 * Generate comparison visualization
 */
function generateComparison() {
    // Get selected molecules
    const comparisonTargets = document.querySelectorAll('.comparison-target:not(.empty)');
    if (comparisonTargets.length < 2) {
        alert('Please select at least two molecules for comparison.');
        return;
    }
    
    // Get selected properties
    const propertySelect = document.getElementById('comparisonPropertiesSelect');
    const selectedProperties = Array.from(propertySelect.selectedOptions).map(option => option.value);
    
    if (selectedProperties.length === 0) {
        alert('Please select at least one property for comparison.');
        return;
    }
    
    // Get comparison type
    const comparisonType = document.getElementById('comparisonTypeSelect').value;
    
    // Generate comparison visualization
    const comparisonResults = document.getElementById('comparisonResults');
    
    // For demo purposes, create a simple radar chart
    comparisonResults.innerHTML = '<canvas id="comparisonChart"></canvas>';
    
    const labels = selectedProperties.map(prop => {
        switch (prop) {
            case 'mw': return 'Molecular Weight';
            case 'logp': return 'LogP';
            case 'tpsa': return 'TPSA';
            case 'hba': return 'H-Bond Acceptors';
            case 'hbd': return 'H-Bond Donors';
            case 'rotb': return 'Rotatable Bonds';
            case 'cypc': return 'Cryoprotection Score';
            default: return prop;
        }
    });
    
    // Generate random data for demo
    const datasets = [];
    const colors = ['rgb(54, 162, 235)', 'rgb(255, 99, 132)', 'rgb(255, 205, 86)', 'rgb(75, 192, 192)'];
    
    comparisonTargets.forEach((target, index) => {
        const moleculeId = target.dataset.moleculeId;
        const moleculeName = target.querySelector('.fw-bold').textContent;
        
        // Generate random values between 0 and 1 for demo
        const data = selectedProperties.map(_ => Math.random().toFixed(2));
        
        datasets.push({
            label: moleculeName,
            data: data,
            backgroundColor: `${colors[index % colors.length]}33`,
            borderColor: colors[index % colors.length],
            pointBackgroundColor: colors[index % colors.length],
            pointHoverBackgroundColor: '#fff',
            pointHoverBorderColor: colors[index % colors.length]
        });
    });
    
    // Create chart
    const ctx = document.getElementById('comparisonChart').getContext('2d');
    new Chart(ctx, {
        type: 'radar',
        data: {
            labels: labels,
            datasets: datasets
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            elements: {
                line: {
                    borderWidth: 2
                }
            },
            scales: {
                r: {
                    angleLines: {
                        display: true
                    },
                    suggestedMin: 0,
                    suggestedMax: 1
                }
            }
        }
    });
}

/**
 * Clear comparison
 */
function clearComparison() {
    const comparisonTargets = document.querySelectorAll('.comparison-target');
    
    comparisonTargets.forEach(target => {
        target.innerHTML = `
            <div class="text-muted">
                <i class="bi bi-plus-circle"></i>
                <div>Drop molecule here</div>
            </div>
        `;
        target.classList.add('empty');
        delete target.dataset.moleculeId;
    });
    
    // Clear comparison results
    document.getElementById('comparisonResults').innerHTML = `
        <div class="text-center py-5">
            <p class="text-muted">Add molecules and generate a comparison to see results</p>
        </div>
    `;
}

/**
 * Initialize correlation analysis
 */
function initializeCorrelationAnalysis() {
    // Calculate correlation button
    document.getElementById('calculateCorrelationBtn').addEventListener('click', function() {
        calculateCorrelation();
    });
    
    // Dataset select change
    document.getElementById('datasetSelect').addEventListener('change', function() {
        if (this.value === 'custom') {
            document.getElementById('customSelectionContainer').classList.remove('d-none');
        } else {
            document.getElementById('customSelectionContainer').classList.add('d-none');
        }
    });
}

/**
 * Calculate correlation
 */
function calculateCorrelation() {
    // Get selected properties
    const propertySelect = document.getElementById('correlationPropertiesSelect');
    const selectedProperties = Array.from(propertySelect.selectedOptions).map(option => option.value);
    
    if (selectedProperties.length < 2) {
        alert('Please select at least two properties for correlation analysis.');
        return;
    }
    
    // Get correlation method
    const correlationMethod = document.getElementById('correlationMethodSelect').value;
    
    // Get dataset
    const dataset = document.getElementById('datasetSelect').value;
    
    // For demo purposes, create a simple correlation matrix
    const correlationMatrix = document.getElementById('correlationMatrix');
    correlationMatrix.innerHTML = '';
    
    // Create a table for the correlation matrix
    const table = document.createElement('table');
    table.className = 'table table-bordered correlation-table';
    
    // Create header row
    const thead = document.createElement('thead');
    const headerRow = document.createElement('tr');
    headerRow.appendChild(document.createElement('th')); // Empty corner cell
    
    selectedProperties.forEach(prop => {
        const th = document.createElement('th');
        th.textContent = getPropertyDisplayName(prop);
        headerRow.appendChild(th);
    });
    
    thead.appendChild(headerRow);
    table.appendChild(thead);
    
    // Create body rows
    const tbody = document.createElement('tbody');
    
    selectedProperties.forEach(rowProp => {
        const row = document.createElement('tr');
        
        // Row header
        const rowHeader = document.createElement('th');
        rowHeader.textContent = getPropertyDisplayName(rowProp);
        row.appendChild(rowHeader);
        
        // Correlation values
        selectedProperties.forEach(colProp => {
            const cell = document.createElement('td');
            
            if (rowProp === colProp) {
                // Diagonal is always 1
                cell.textContent = '1.00';
                cell.style.backgroundColor = 'rgb(0, 128, 0, 0.5)';
            } else {
                // Generate random correlation value for demo
                const value = (Math.random() * 2 - 1).toFixed(2);
                const absValue = Math.abs(parseFloat(value));
                
                // Color based on correlation strength
                let color;
                if (value > 0) {
                    color = `rgba(0, 128, 0, ${absValue})`;
                } else {
                    color = `rgba(255, 0, 0, ${absValue})`;
                }
                
                cell.textContent = value;
                cell.style.backgroundColor = color;
            }
            
            row.appendChild(cell);
        });
        
        tbody.appendChild(row);
    });
    
    table.appendChild(tbody);
    correlationMatrix.appendChild(table);
    
    // Generate correlation analysis
    generateCorrelationAnalysis(selectedProperties, correlationMethod, dataset);
}

/**
 * Generate correlation analysis
 * @param {Array} properties - Selected properties
 * @param {string} method - Correlation method
 * @param {string} dataset - Selected dataset
 */
function generateCorrelationAnalysis(properties, method, dataset) {
    const analysisContainer = document.getElementById('correlationAnalysis');
    
    // Generate a simple analysis report
    analysisContainer.innerHTML = `
        <h6>Correlation Analysis</h6>
        <p>Analysis of ${properties.length} properties using ${method} correlation on the ${dataset} dataset.</p>
        
        <h6>Key Findings</h6>
        <ul>
            <li>Strong positive correlation between Molecular Weight and TPSA (0.78)</li>
            <li>Moderate negative correlation between LogP and Cryoprotection Score (-0.62)</li>
            <li>Weak correlation between H-Bond Donors and Rotatable Bonds (0.23)</li>
        </ul>
        
        <h6>Statistical Significance</h6>
        <p>All reported correlations are statistically significant (p-value < 0.05) except for the correlation between H-Bond Donors and Rotatable Bonds (p-value = 0.09).</p>
        
        <h6>Recommendations</h6>
        <p>Based on the analysis, molecules with lower LogP values tend to have higher cryoprotection scores. This suggests that hydrophilic compounds may be better cryoprotectants.</p>
    `;
}

/**
 * Initialize export tool
 */
function initializeExportTool() {
    // Generate export button
    document.getElementById('generateExportBtn').addEventListener('click', function() {
        generateExport();
    });
    
    // Download export button
    document.getElementById('downloadExportBtn').addEventListener('click', function() {
        downloadExport();
    });
}

/**
 * Generate export preview
 */
function generateExport() {
    // Get export settings
    const dataset = document.getElementById('exportDatasetSelect').value;
    const format = document.getElementById('exportFormatSelect').value;
    const propertiesSelect = document.getElementById('exportPropertiesSelect');
    const selectedProperties = Array.from(propertiesSelect.selectedOptions).map(option => option.value);
    
    if (selectedProperties.length === 0) {
        alert('Please select at least one property to include in the export.');
        return;
    }
    
    // Generate export preview
    const exportPreviewGrid = document.getElementById('exportPreviewGrid');
    
    // Prepare data for preview
    const columns = selectedProperties.map(prop => ({
        id: prop,
        name: getPropertyDisplayName(prop)
    }));
    
    // Generate sample data
    const data = [
        { id: 1, name: 'Ethanol', smiles: 'CCO', mw: 46.07, logp: -0.31, tpsa: 20.23, hba: 1, hbd: 1, rotb: 0, cypc: 0.78 },
        { id: 2, name: 'Glycerol', smiles: 'C(C(CO)O)O', mw: 92.09, logp: -1.76, tpsa: 60.69, hba: 3, hbd: 3, rotb: 2, cypc: 0.92 },
        { id: 3, name: 'Dimethyl Sulfoxide', smiles: 'CS(=O)C', mw: 78.13, logp: -1.35, tpsa: 17.07, hba: 1, hbd: 0, rotb: 0, cypc: 0.86 },
        { id: 4, name: 'Propylene Glycol', smiles: 'CC(CO)O', mw: 76.10, logp: -0.92, tpsa: 40.46, hba: 2, hbd: 2, rotb: 1, cypc: 0.81 },
        { id: 5, name: 'Trehalose', smiles: 'C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@]2([C@H]([C@@H]([C@H](O2)CO)O)O)CO)O)O)O)O', mw: 342.30, logp: -4.23, tpsa: 189.53, hba: 11, hbd: 8, rotb: 2, cypc: 0.89 }
    ];
    
    // Filter data based on selected properties
    const filteredData = data.map(item => {
        const filteredItem = {};
        selectedProperties.forEach(prop => {
            filteredItem[prop] = item[prop];
        });
        return filteredItem;
    });
    
    // Initialize GridJS
    const grid = new gridjs.Grid({
        columns: columns,
        data: filteredData,
        pagination: { limit: 10 },
        sort: true,
        className: {
            table: 'table table-striped table-hover'
        }
    }).render(exportPreviewGrid);
    
    // Update record count
    document.getElementById('exportRecordCount').textContent = filteredData.length;
    
    // Enable download button
    document.getElementById('downloadExportBtn').disabled = false;
}

/**
 * Download export
 */
function downloadExport() {
    const format = document.getElementById('exportFormatSelect').value;
    const filename = `cryoprotect_export_${new Date().toISOString().split('T')[0]}.${format}`;
    
    // In a real implementation, this would generate the actual export file
    alert(`Export file ${filename} would be downloaded here.`);
}

/**
 * Set up event listeners
 */
function setupEventListeners() {
    // Apply filters button
    document.getElementById('applyFilters').addEventListener('click', function() {
        applyFilters();
    });
    
    // Reset filters button
    document.getElementById('resetFilters').addEventListener('click', function() {
        resetFilters();
    });
    
    // Export visualization button
    document.getElementById('exportVisualization').addEventListener('click', function() {
        exportVisualization();
    });
    
    // Export results button
    document.getElementById('exportResults').addEventListener('click', function() {
        exportResults();
    });
    
    // Export comparison button
    document.getElementById('exportComparisonBtn').addEventListener('click', function() {
        exportComparison();
    });
    
    // Export correlation button
    document.getElementById('exportCorrelationBtn').addEventListener('click', function() {
        exportCorrelation();
    });
}

/**
 * Apply filters
 */
function applyFilters() {
    alert('Filters applied. Visualization and results would be updated in a real implementation.');
}

/**
 * Reset filters
 */
function resetFilters() {
    // Reset property type checkboxes
    document.querySelectorAll('#propertyTypeFilters input[type="checkbox"]').forEach(checkbox => {
        checkbox.checked = true;
    });
    
    // Reset value range inputs
    document.getElementById('mwMin').value = '';
    document.getElementById('mwMax').value = '';
    document.getElementById('mwRange').value = 500;
    
    document.getElementById('logpMin').value = '';
    document.getElementById('logpMax').value = '';
    document.getElementById('logpRange').value = 0;
    
    // Reset classification checkboxes
    document.getElementById('cryoprotectantsOnly').checked = false;
    document.getElementById('highPerformanceOnly').checked = false;
    
    // Reset structure search
    document.getElementById('substructureSearch').value = '';
    
    alert('Filters reset. Visualization and results would be updated in a real implementation.');
}

/**
 * Export visualization
 */
function exportVisualization() {
    alert('Visualization would be exported in a real implementation.');
}

/**
 * Export results
 */
function exportResults() {
    alert('Results would be exported in a real implementation.');
}

/**
 * Export comparison
 */
function exportComparison() {
    alert('Comparison would be exported in a real implementation.');
}

/**
 * Export correlation
 */
function exportCorrelation() {
    alert('Correlation analysis would be exported in a real implementation.');
}

/**
 * Get property display name
 * @param {string} propertyId - Property ID
 * @returns {string} Display name
 */
function getPropertyDisplayName(propertyId) {
    const propertyNames = {
        'id': 'ID',
        'name': 'Name',
        'smiles': 'SMILES',
        'mw': 'Molecular Weight',
        'logp': 'LogP',
        'tpsa': 'TPSA',
        'hba': 'H-Bond Acceptors',
        'hbd': 'H-Bond Donors',
        'rotb': 'Rotatable Bonds',
        'cypc': 'Cryoprotection Score',
        'fsp3': 'Fraction sp3',
        'aroRings': 'Aromatic Rings',
        'qed': 'QED Score',
        'heavyAtoms': 'Heavy Atoms',
        'toxicity': 'Toxicity Score'
    };
    
    return propertyNames[propertyId] || propertyId;
}
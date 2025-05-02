/**
 * CryoProtect Analyzer - Dashboard Module
 * 
 * This module provides functions for creating an interactive dashboard with key metrics and visualizations
 */

const Dashboard = (function() {
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
    
    // Cache for dashboard data
    let dashboardData = null;
    
    /**
     * Initialize the dashboard
     * 
     * @param {string} containerId - ID of the container element
     * @param {Object} options - Dashboard options
     */
    function initDashboard(containerId, options = {}) {
        const container = document.getElementById(containerId);
        if (!container) {
            console.error(`Container element with ID ${containerId} not found`);
            return;
        }
        
        // Default options
        const defaultOptions = {
            title: 'CryoProtect Analyzer Dashboard',
            refreshInterval: 0, // 0 means no auto-refresh
            layout: 'grid', // 'grid' or 'fluid'
            showFilters: true
        };
        
        // Merge options
        const dashboardOptions = { ...defaultOptions, ...options };
        
        // Create dashboard structure
        createDashboardStructure(container, dashboardOptions);
        
        // Load dashboard data
        loadDashboardData(dashboardOptions);
        
        // Set up auto-refresh if enabled
        if (dashboardOptions.refreshInterval > 0) {
            setInterval(() => {
                loadDashboardData(dashboardOptions);
            }, dashboardOptions.refreshInterval * 1000);
        }
    }
    
    /**
     * Create the dashboard structure
     * 
     * @param {HTMLElement} container - Container element
     * @param {Object} options - Dashboard options
     */
    function createDashboardStructure(container, options) {
        // Clear container
        container.innerHTML = '';
        
        // Create dashboard header
        const header = document.createElement('div');
        header.className = 'dashboard-header mb-4';
        
        const title = document.createElement('h1');
        title.textContent = options.title;
        header.appendChild(title);
        
        // Create filters if enabled
        if (options.showFilters) {
            const filters = document.createElement('div');
            filters.className = 'dashboard-filters mb-4';
            
            const filtersRow = document.createElement('div');
            filtersRow.className = 'row g-3';
            
            // Date range filter
            const dateFilterCol = document.createElement('div');
            dateFilterCol.className = 'col-md-4';
            
            const dateFilterLabel = document.createElement('label');
            dateFilterLabel.className = 'form-label';
            dateFilterLabel.textContent = 'Date Range';
            dateFilterCol.appendChild(dateFilterLabel);
            
            const dateFilterSelect = document.createElement('select');
            dateFilterSelect.className = 'form-select';
            dateFilterSelect.id = 'dashboard-date-filter';
            
            const dateOptions = [
                { value: 'all', text: 'All Time' },
                { value: 'today', text: 'Today' },
                { value: 'week', text: 'This Week' },
                { value: 'month', text: 'This Month' },
                { value: 'year', text: 'This Year' }
            ];
            
            dateOptions.forEach(option => {
                const optionElement = document.createElement('option');
                optionElement.value = option.value;
                optionElement.textContent = option.text;
                dateFilterSelect.appendChild(optionElement);
            });
            
            dateFilterCol.appendChild(dateFilterSelect);
            filtersRow.appendChild(dateFilterCol);
            
            // Category filter
            const categoryFilterCol = document.createElement('div');
            categoryFilterCol.className = 'col-md-4';
            
            const categoryFilterLabel = document.createElement('label');
            categoryFilterLabel.className = 'form-label';
            categoryFilterLabel.textContent = 'Category';
            categoryFilterCol.appendChild(categoryFilterLabel);
            
            const categoryFilterSelect = document.createElement('select');
            categoryFilterSelect.className = 'form-select';
            categoryFilterSelect.id = 'dashboard-category-filter';
            
            const categoryOptions = [
                { value: 'all', text: 'All Categories' },
                { value: 'molecules', text: 'Molecules' },
                { value: 'mixtures', text: 'Mixtures' },
                { value: 'protocols', text: 'Protocols' },
                { value: 'experiments', text: 'Experiments' }
            ];
            
            categoryOptions.forEach(option => {
                const optionElement = document.createElement('option');
                optionElement.value = option.value;
                optionElement.textContent = option.text;
                categoryFilterSelect.appendChild(optionElement);
            });
            
            categoryFilterCol.appendChild(categoryFilterSelect);
            filtersRow.appendChild(categoryFilterCol);
            
            // Apply button
            const applyButtonCol = document.createElement('div');
            applyButtonCol.className = 'col-md-4 d-flex align-items-end';
            
            const applyButton = document.createElement('button');
            applyButton.className = 'btn btn-primary';
            applyButton.textContent = 'Apply Filters';
            applyButton.id = 'dashboard-apply-filters';
            applyButtonCol.appendChild(applyButton);
            
            filtersRow.appendChild(applyButtonCol);
            filters.appendChild(filtersRow);
            
            // Add event listener to apply button
            applyButton.addEventListener('click', () => {
                loadDashboardData(options);
            });
            
            header.appendChild(filters);
        }
        
        container.appendChild(header);
        
        // Create dashboard content
        const content = document.createElement('div');
        content.className = options.layout === 'grid' ? 'dashboard-content' : 'dashboard-content-fluid';
        content.id = 'dashboard-content';
        
        // Create loading indicator
        const loading = document.createElement('div');
        loading.className = 'text-center py-5';
        loading.id = 'dashboard-loading';
        
        const spinner = document.createElement('div');
        spinner.className = 'spinner-border';
        spinner.setAttribute('role', 'status');
        
        const spinnerText = document.createElement('span');
        spinnerText.className = 'visually-hidden';
        spinnerText.textContent = 'Loading...';
        
        spinner.appendChild(spinnerText);
        loading.appendChild(spinner);
        
        const loadingText = document.createElement('p');
        loadingText.className = 'mt-2';
        loadingText.textContent = 'Loading dashboard data...';
        loading.appendChild(loadingText);
        
        content.appendChild(loading);
        container.appendChild(content);
    }
    
    /**
     * Load dashboard data
     * 
     * @param {Object} options - Dashboard options
     */
    function loadDashboardData(options) {
        const loading = document.getElementById('dashboard-loading');
        if (loading) loading.style.display = 'block';
        
        // Get filter values
        const dateFilter = document.getElementById('dashboard-date-filter');
        const categoryFilter = document.getElementById('dashboard-category-filter');
        
        const filters = {
            dateRange: dateFilter ? dateFilter.value : 'all',
            category: categoryFilter ? categoryFilter.value : 'all'
        };
        
        // Fetch dashboard data from API
        fetch('/api/v1/dashboard', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json'
            },
            body: JSON.stringify(filters)
        })
        .then(response => {
            if (!response.ok) {
                throw new Error(`Error fetching dashboard data: ${response.statusText}`);
            }
            return response.json();
        })
        .then(data => {
            // Cache the data
            dashboardData = data;
            
            // Render the dashboard
            renderDashboard(data, options);
        })
        .catch(error => {
            console.error('Error loading dashboard data:', error);
            
            // Show error message
            const content = document.getElementById('dashboard-content');
            if (content) {
                content.innerHTML = `
                    <div class="alert alert-danger">
                        <h4 class="alert-heading">Error Loading Dashboard</h4>
                        <p>${error.message}</p>
                        <hr>
                        <p class="mb-0">Please try again later or contact support if the problem persists.</p>
                    </div>
                `;
            }
        })
        .finally(() => {
            if (loading) loading.style.display = 'none';
        });
    }
    
    /**
     * Render the dashboard with data
     * 
     * @param {Object} data - Dashboard data
     * @param {Object} options - Dashboard options
     */
    function renderDashboard(data, options) {
        const content = document.getElementById('dashboard-content');
        if (!content) return;
        
        // Clear content (except loading indicator)
        const loading = document.getElementById('dashboard-loading');
        content.innerHTML = '';
        if (loading) content.appendChild(loading);
        
        // Create grid layout
        const grid = document.createElement('div');
        grid.className = 'row g-4';
        
        // Add key metrics
        if (data.keyMetrics) {
            createKeyMetricsSection(grid, data.keyMetrics);
        }
        
        // Add molecule visualization
        if (data.topMolecules && data.topMolecules.length > 0) {
            createMoleculeVisualizationSection(grid, data.topMolecules);
        }
        
        // Add property distribution chart
        if (data.propertyDistribution) {
            createPropertyDistributionSection(grid, data.propertyDistribution);
        }
        
        // Add mixture composition chart
        if (data.topMixtures && data.topMixtures.length > 0) {
            createMixtureCompositionSection(grid, data.topMixtures);
        }
        
        // Add protocol timeline chart
        if (data.recentProtocols && data.recentProtocols.length > 0) {
            createProtocolTimelineSection(grid, data.recentProtocols);
        }
        
        content.appendChild(grid);
    }
    
    /**
     * Create key metrics section
     * 
     * @param {HTMLElement} container - Container element
     * @param {Object} metrics - Key metrics data
     */
    function createKeyMetricsSection(container, metrics) {
        const metricsRow = document.createElement('div');
        metricsRow.className = 'col-12';
        
        const metricsCard = document.createElement('div');
        metricsCard.className = 'card';
        
        const metricsHeader = document.createElement('div');
        metricsHeader.className = 'card-header';
        metricsHeader.textContent = 'Key Metrics';
        metricsCard.appendChild(metricsHeader);
        
        const metricsBody = document.createElement('div');
        metricsBody.className = 'card-body';
        
        const metricsGrid = document.createElement('div');
        metricsGrid.className = 'row g-3';
        
        // Create metric cards
        Object.entries(metrics).forEach(([key, value]) => {
            const metricCol = document.createElement('div');
            metricCol.className = 'col-md-3 col-sm-6';
            
            const metricCard = document.createElement('div');
            metricCard.className = 'card h-100 border-0 shadow-sm';
            
            const metricBody = document.createElement('div');
            metricBody.className = 'card-body text-center';
            
            const metricValue = document.createElement('h2');
            metricValue.className = 'display-4 fw-bold';
            metricValue.textContent = value.value;
            metricBody.appendChild(metricValue);
            
            const metricLabel = document.createElement('p');
            metricLabel.className = 'text-muted mb-0';
            metricLabel.textContent = value.label;
            metricBody.appendChild(metricLabel);
            
            metricCard.appendChild(metricBody);
            metricCol.appendChild(metricCard);
            metricsGrid.appendChild(metricCol);
        });
        
        metricsBody.appendChild(metricsGrid);
        metricsCard.appendChild(metricsBody);
        metricsRow.appendChild(metricsCard);
        container.appendChild(metricsRow);
    }
    
    /**
     * Create molecule visualization section
     * 
     * @param {HTMLElement} container - Container element
     * @param {Array} molecules - Top molecules data
     */
    function createMoleculeVisualizationSection(container, molecules) {
        const moleculeRow = document.createElement('div');
        moleculeRow.className = 'col-md-6';
        
        const moleculeCard = document.createElement('div');
        moleculeCard.className = 'card h-100';
        
        const moleculeHeader = document.createElement('div');
        moleculeHeader.className = 'card-header';
        moleculeHeader.textContent = 'Top Molecules';
        moleculeCard.appendChild(moleculeHeader);
        
        const moleculeBody = document.createElement('div');
        moleculeBody.className = 'card-body';
        
        // Create molecule viewer container
        const moleculeViewer = document.createElement('div');
        moleculeViewer.id = 'top-molecule-viewer';
        moleculeViewer.className = 'molecule-viewer';
        moleculeViewer.style.height = '300px';
        moleculeBody.appendChild(moleculeViewer);
        
        // Create molecule selector
        const moleculeSelector = document.createElement('div');
        moleculeSelector.className = 'mt-3';
        
        const moleculeLabel = document.createElement('label');
        moleculeLabel.className = 'form-label';
        moleculeLabel.textContent = 'Select Molecule';
        moleculeSelector.appendChild(moleculeLabel);
        
        const moleculeSelect = document.createElement('select');
        moleculeSelect.className = 'form-select';
        moleculeSelect.id = 'molecule-selector';
        
        molecules.forEach((molecule, index) => {
            const option = document.createElement('option');
            option.value = index;
            option.textContent = molecule.name || `CID: ${molecule.cid}`;
            moleculeSelect.appendChild(option);
        });
        
        moleculeSelector.appendChild(moleculeSelect);
        moleculeBody.appendChild(moleculeSelector);
        
        // Create molecule toolbar
        const moleculeToolbar = document.createElement('div');
        moleculeToolbar.id = 'molecule-viewer-toolbar';
        moleculeToolbar.className = 'btn-toolbar mt-3';
        moleculeBody.appendChild(moleculeToolbar);
        
        moleculeCard.appendChild(moleculeBody);
        moleculeRow.appendChild(moleculeCard);
        container.appendChild(moleculeRow);
        
        // Initialize molecule viewer after rendering
        setTimeout(() => {
            if (molecules[0].smiles) {
                // Initialize 3D viewer
                MolecularViewer.initViewer('top-molecule-viewer');
                MolecularViewer.loadMoleculeFromSmiles('top-molecule-viewer', molecules[0].smiles);
                MolecularViewer.createViewerToolbar('top-molecule-viewer', 'molecule-viewer-toolbar');
                
                // Add event listener to molecule selector
                moleculeSelect.addEventListener('change', (event) => {
                    const index = parseInt(event.target.value);
                    if (molecules[index] && molecules[index].smiles) {
                        MolecularViewer.loadMoleculeFromSmiles('top-molecule-viewer', molecules[index].smiles);
                    }
                });
            }
        }, 100);
    }
    
    /**
     * Create property distribution section
     * 
     * @param {HTMLElement} container - Container element
     * @param {Object} propertyDistribution - Property distribution data
     */
    function createPropertyDistributionSection(container, propertyDistribution) {
        const propertyRow = document.createElement('div');
        propertyRow.className = 'col-md-6';
        
        const propertyCard = document.createElement('div');
        propertyCard.className = 'card h-100';
        
        const propertyHeader = document.createElement('div');
        propertyHeader.className = 'card-header d-flex justify-content-between align-items-center';
        
        const propertyTitle = document.createElement('span');
        propertyTitle.textContent = 'Property Distribution';
        propertyHeader.appendChild(propertyTitle);
        
        // Property selector
        const propertySelect = document.createElement('select');
        propertySelect.className = 'form-select form-select-sm w-auto';
        propertySelect.id = 'property-selector';
        
        Object.keys(propertyDistribution).forEach(property => {
            const option = document.createElement('option');
            option.value = property;
            option.textContent = property;
            propertySelect.appendChild(option);
        });
        
        propertyHeader.appendChild(propertySelect);
        propertyCard.appendChild(propertyHeader);
        
        const propertyBody = document.createElement('div');
        propertyBody.className = 'card-body';
        
        // Create chart container
        const propertyChart = document.createElement('div');
        propertyChart.id = 'property-distribution-chart';
        propertyChart.style.height = '300px';
        propertyBody.appendChild(propertyChart);
        
        propertyCard.appendChild(propertyBody);
        propertyRow.appendChild(propertyCard);
        container.appendChild(propertyRow);
        
        // Initialize chart after rendering
        setTimeout(() => {
            const firstProperty = Object.keys(propertyDistribution)[0];
            if (firstProperty) {
                const propertyData = propertyDistribution[firstProperty];
                
                // Create histogram
                InteractiveCharts.createBarChart(
                    'property-distribution-chart',
                    propertyData.map((value, index) => ({ bin: propertyData.labels[index], count: value })),
                    'bin',
                    'count',
                    {
                        title: `Distribution of ${firstProperty}`,
                        xAxisTitle: firstProperty,
                        yAxisTitle: 'Count',
                        color: COLORS[1]
                    }
                );
                
                // Add event listener to property selector
                propertySelect.addEventListener('change', (event) => {
                    const property = event.target.value;
                    if (propertyDistribution[property]) {
                        const propertyData = propertyDistribution[property];
                        
                        // Update chart
                        InteractiveCharts.createBarChart(
                            'property-distribution-chart',
                            propertyData.map((value, index) => ({ bin: propertyData.labels[index], count: value })),
                            'bin',
                            'count',
                            {
                                title: `Distribution of ${property}`,
                                xAxisTitle: property,
                                yAxisTitle: 'Count',
                                color: COLORS[1]
                            }
                        );
                    }
                });
            }
        }, 100);
    }
    
    /**
     * Create mixture composition section
     * 
     * @param {HTMLElement} container - Container element
     * @param {Array} mixtures - Top mixtures data
     */
    function createMixtureCompositionSection(container, mixtures) {
        const mixtureRow = document.createElement('div');
        mixtureRow.className = 'col-md-6';
        
        const mixtureCard = document.createElement('div');
        mixtureCard.className = 'card h-100';
        
        const mixtureHeader = document.createElement('div');
        mixtureHeader.className = 'card-header d-flex justify-content-between align-items-center';
        
        const mixtureTitle = document.createElement('span');
        mixtureTitle.textContent = 'Mixture Compositions';
        mixtureHeader.appendChild(mixtureTitle);
        
        // Mixture selector
        const mixtureSelect = document.createElement('select');
        mixtureSelect.className = 'form-select form-select-sm w-auto';
        mixtureSelect.id = 'mixture-selector';
        
        mixtures.forEach((mixture, index) => {
            const option = document.createElement('option');
            option.value = index;
            option.textContent = mixture.name || `Mixture ${index + 1}`;
            mixtureSelect.appendChild(option);
        });
        
        mixtureHeader.appendChild(mixtureSelect);
        mixtureCard.appendChild(mixtureHeader);
        
        const mixtureBody = document.createElement('div');
        mixtureBody.className = 'card-body';
        
        // Create chart container
        const mixtureChart = document.createElement('div');
        mixtureChart.id = 'mixture-composition-chart';
        mixtureChart.style.height = '300px';
        mixtureBody.appendChild(mixtureChart);
        
        mixtureCard.appendChild(mixtureBody);
        mixtureRow.appendChild(mixtureCard);
        container.appendChild(mixtureRow);
        
        // Initialize chart after rendering
        setTimeout(() => {
            if (mixtures[0].components) {
                // Create pie chart
                MixtureVisualizer.createCompositionPieChart(
                    'mixture-composition-chart',
                    mixtures[0]
                );
                
                // Add event listener to mixture selector
                mixtureSelect.addEventListener('change', (event) => {
                    const index = parseInt(event.target.value);
                    if (mixtures[index] && mixtures[index].components) {
                        MixtureVisualizer.createCompositionPieChart(
                            'mixture-composition-chart',
                            mixtures[index]
                        );
                    }
                });
            }
        }, 100);
    }
    
    /**
     * Create protocol timeline section
     * 
     * @param {HTMLElement} container - Container element
     * @param {Array} protocols - Recent protocols data
     */
    function createProtocolTimelineSection(container, protocols) {
        const protocolRow = document.createElement('div');
        protocolRow.className = 'col-md-6';
        
        const protocolCard = document.createElement('div');
        protocolCard.className = 'card h-100';
        
        const protocolHeader = document.createElement('div');
        protocolHeader.className = 'card-header d-flex justify-content-between align-items-center';
        
        const protocolTitle = document.createElement('span');
        protocolTitle.textContent = 'Protocol Timeline';
        protocolHeader.appendChild(protocolTitle);
        
        // Protocol selector
        const protocolSelect = document.createElement('select');
        protocolSelect.className = 'form-select form-select-sm w-auto';
        protocolSelect.id = 'protocol-selector';
        
        protocols.forEach((protocol, index) => {
            const option = document.createElement('option');
            option.value = index;
            option.textContent = protocol.name || `Protocol ${index + 1}`;
            protocolSelect.appendChild(option);
        });
        
        protocolHeader.appendChild(protocolSelect);
        protocolCard.appendChild(protocolHeader);
        
        const protocolBody = document.createElement('div');
        protocolBody.className = 'card-body';
        
        // Create chart container
        const protocolChart = document.createElement('div');
        protocolChart.id = 'protocol-timeline-chart';
        protocolChart.style.height = '300px';
        protocolBody.appendChild(protocolChart);
        
        protocolCard.appendChild(protocolBody);
        protocolRow.appendChild(protocolCard);
        container.appendChild(protocolRow);
        
        // Initialize chart after rendering
        setTimeout(() => {
            if (protocols[0].steps) {
                // Create timeline chart
                ProtocolVisualizer.createProtocolTimeline(
                    'protocol-timeline-chart',
                    protocols[0]
                );
                
                // Add event listener to protocol selector
                protocolSelect.addEventListener('change', (event) => {
                    const index = parseInt(event.target.value);
                    if (protocols[index] && protocols[index].steps) {
                        ProtocolVisualizer.createProtocolTimeline(
                            'protocol-timeline-chart',
                            protocols[index]
                        );
                    }
                });
            }
        }, 100);
    }
    
    // Public API
    return {
        initDashboard
    };
})();

// Export the module
if (typeof module !== 'undefined' && module.exports) {
    module.exports = Dashboard;
} else {
    window.Dashboard = Dashboard;
}

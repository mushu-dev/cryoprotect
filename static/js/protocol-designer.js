/**
 * CryoProtect Analyzer - Protocol Designer JavaScript
 * 
 * This file contains the JavaScript code for the protocol designer functionality.
 */

// Protocol Designer Module
const ProtocolDesigner = (function() {
    // Cache DOM elements
    let $protocolForm;
    let $mixtureSelect;
    let $targetConcentration;
    let $sampleTypeSelect;
    let $startingTemperature;
    let $targetTemperature;
    let $stepCount;
    let $designButton;
    let $saveButton;
    let $protocolResults;
    let $protocolVisualization;
    let $protocolSteps;
    let $protocolList;
    let $compareButton;
    
    // Current protocol data
    let currentProtocol = null;
    let selectedProtocols = [];
    
    // Initialize the module
    function init() {
        // Cache DOM elements
        $protocolForm = $('#protocol-form');
        $mixtureSelect = $('#mixture-select');
        $targetConcentration = $('#target-concentration');
        $sampleTypeSelect = $('#sample-type-select');
        $startingTemperature = $('#starting-temperature');
        $targetTemperature = $('#target-temperature');
        $stepCount = $('#step-count');
        $designButton = $('#design-protocol-btn');
        $saveButton = $('#save-protocol-btn');
        $protocolResults = $('#protocol-results');
        $protocolVisualization = $('#protocol-visualization');
        $protocolSteps = $('#protocol-steps');
        $protocolList = $('#protocol-list');
        $compareButton = $('#compare-protocols-btn');
        
        // Bind events
        $designButton.on('click', designProtocol);
        $saveButton.on('click', saveProtocol);
        $compareButton.on('click', compareProtocols);
        
        // Load mixtures
        loadMixtures();
        
        // Load sample types
        loadSampleTypes();
        
        // Initialize tooltips
        initTooltips();
    }
    
    // Load mixtures from API
    function loadMixtures() {
        $.ajax({
            url: '/api/v1/mixtures',
            method: 'GET',
            success: function(response) {
                $mixtureSelect.empty();
                
                if (response && response.length > 0) {
                    response.forEach(function(mixture) {
                        $mixtureSelect.append(
                            $('<option></option>')
                                .attr('value', mixture.id)
                                .text(mixture.name)
                        );
                    });
                    
                    // Load protocols for the first mixture
                    loadProtocols(response[0].id);
                } else {
                    $mixtureSelect.append(
                        $('<option></option>')
                            .attr('value', '')
                            .text('No mixtures available')
                    );
                }
            },
            error: function(xhr, status, error) {
                console.error('Error loading mixtures:', error);
                showAlert('error', 'Failed to load mixtures. Please try again.');
            }
        });
        
        // Add change event to load protocols when mixture changes
        $mixtureSelect.on('change', function() {
            const mixtureId = $(this).val();
            if (mixtureId) {
                loadProtocols(mixtureId);
            }
        });
    }
    
    // Load sample types from API
    function loadSampleTypes() {
        $.ajax({
            url: '/api/v1/protocols/sensitivity-profiles',
            method: 'GET',
            success: function(response) {
                $sampleTypeSelect.empty();
                
                if (response && response.profiles) {
                    Object.keys(response.profiles).forEach(function(sampleType) {
                        const profile = response.profiles[sampleType];
                        const displayName = profile.name || sampleType
                            .split('_')
                            .map(word => word.charAt(0).toUpperCase() + word.slice(1))
                            .join(' ');
                        
                        // Create tooltip content with detailed profile information
                        let tooltipContent = `
                            <strong>${displayName}</strong><br>
                            Osmotic Tolerance: ${profile.osmotic_tolerance * 100}%<br>
                            Max Step Size: ${profile.max_step_size}%<br>
                            Time per Step: ${profile.time_per_step} min<br>
                            Cooling Rate: ${profile.cooling_rate}°C/min
                        `;
                        
                        if (profile.notes) {
                            tooltipContent += `<br><em>${profile.notes}</em>`;
                        }
                        
                        $sampleTypeSelect.append(
                            $('<option></option>')
                                .attr('value', sampleType)
                                .text(displayName)
                                .attr('data-toggle', 'tooltip')
                                .attr('data-html', 'true')
                                .attr('title', tooltipContent)
                        );
                    });
                    
                    // Add event to show sensitivity details when sample type changes
                    $sampleTypeSelect.on('change', function() {
                        const selectedType = $(this).val();
                        const profile = response.profiles[selectedType];
                        
                        if (profile) {
                            // Update form with recommended values based on profile
                            if (profile.max_step_size) {
                                const targetConc = parseFloat($targetConcentration.val()) || 10.0;
                                const recommendedSteps = Math.ceil(targetConc / profile.max_step_size);
                                $stepCount.val(recommendedSteps);
                            }
                            
                            // Show sensitivity info in a card below the select
                            const $infoCard = $('#sensitivity-info');
                            if ($infoCard.length === 0) {
                                // Create info card if it doesn't exist
                                const $newCard = $(`
                                    <div id="sensitivity-info" class="card mt-2 mb-3">
                                        <div class="card-body">
                                            <h5 class="card-title">${profile.name || displayName} Sensitivity Profile</h5>
                                            <p class="card-text">${profile.description || ''}</p>
                                            <ul class="list-group list-group-flush">
                                                <li class="list-group-item">Osmotic Tolerance: ${profile.osmotic_tolerance * 100}%</li>
                                                <li class="list-group-item">Max Step Size: ${profile.max_step_size}%</li>
                                                <li class="list-group-item">Time per Step: ${profile.time_per_step} min</li>
                                                <li class="list-group-item">Cooling Rate: ${profile.cooling_rate}°C/min</li>
                                                <li class="list-group-item">Warming Rate: ${profile.warming_rate}°C/min</li>
                                            </ul>
                                            ${profile.notes ? `<p class="card-text text-muted mt-2"><small>${profile.notes}</small></p>` : ''}
                                        </div>
                                    </div>
                                `);
                                $sampleTypeSelect.parent().after($newCard);
                            } else {
                                // Update existing info card
                                $infoCard.find('.card-title').text(`${profile.name || displayName} Sensitivity Profile`);
                                $infoCard.find('.card-text').first().text(profile.description || '');
                                $infoCard.find('.list-group-item').eq(0).text(`Osmotic Tolerance: ${profile.osmotic_tolerance * 100}%`);
                                $infoCard.find('.list-group-item').eq(1).text(`Max Step Size: ${profile.max_step_size}%`);
                                $infoCard.find('.list-group-item').eq(2).text(`Time per Step: ${profile.time_per_step} min`);
                                $infoCard.find('.list-group-item').eq(3).text(`Cooling Rate: ${profile.cooling_rate}°C/min`);
                                $infoCard.find('.list-group-item').eq(4).text(`Warming Rate: ${profile.warming_rate}°C/min`);
                                
                                if (profile.notes) {
                                    if ($infoCard.find('.text-muted').length) {
                                        $infoCard.find('.text-muted').text(profile.notes);
                                    } else {
                                        $infoCard.find('.card-body').append(`<p class="card-text text-muted mt-2"><small>${profile.notes}</small></p>`);
                                    }
                                } else {
                                    $infoCard.find('.text-muted').remove();
                                }
                            }
                        }
                    });
                } else {
                    $sampleTypeSelect.append(
                        $('<option></option>')
                            .attr('value', 'cell_line')
                            .text('Cell Line (Default)')
                    );
                }
                
                // Initialize tooltips
                initTooltips();
                
                // Trigger change event to show initial sensitivity info
                $sampleTypeSelect.trigger('change');
            },
            error: function(xhr, status, error) {
                console.error('Error loading sample types:', error);
                
                // Add default sample types
                $sampleTypeSelect.empty();
                $sampleTypeSelect.append($('<option></option>').attr('value', 'cell_line').text('Cell Line'));
                $sampleTypeSelect.append($('<option></option>').attr('value', 'primary_cells').text('Primary Cells'));
                $sampleTypeSelect.append($('<option></option>').attr('value', 'tissue').text('Tissue'));
                $sampleTypeSelect.append($('<option></option>').attr('value', 'organoid').text('Organoid'));
                $sampleTypeSelect.append($('<option></option>').attr('value', 'embryo').text('Embryo'));
            }
        });
    }
    
    // Load protocols for a mixture
    function loadProtocols(mixtureId) {
        $.ajax({
            url: `/api/v1/mixtures/${mixtureId}/protocols`,
            method: 'GET',
            success: function(response) {
                $protocolList.empty();
                selectedProtocols = [];
                
                if (response && response.protocols && response.protocols.length > 0) {
                    response.protocols.forEach(function(protocol) {
                        const protocolDate = new Date(protocol.created_at).toLocaleString();
                        const $protocolItem = $(
                            `<div class="protocol-item card mb-2">
                                <div class="card-body">
                                    <div class="form-check">
                                        <input class="form-check-input protocol-checkbox" type="checkbox" value="${protocol.protocol_id}" id="protocol-${protocol.protocol_id}">
                                        <label class="form-check-label" for="protocol-${protocol.protocol_id}">
                                            <h5 class="card-title">${protocol.sample_type} - ${protocol.target_concentration}%</h5>
                                            <h6 class="card-subtitle mb-2 text-muted">${protocolDate}</h6>
                                            <p class="card-text">Steps: ${protocol.step_count}, Duration: ${protocol.total_duration} min</p>
                                            <button class="btn btn-sm btn-primary view-protocol-btn" data-protocol-id="${protocol.protocol_id}">View</button>
                                        </label>
                                    </div>
                                </div>
                            </div>`
                        );
                        
                        $protocolList.append($protocolItem);
                    });
                    
                    // Add event listeners for protocol checkboxes
                    $('.protocol-checkbox').on('change', function() {
                        const protocolId = $(this).val();
                        if ($(this).is(':checked')) {
                            selectedProtocols.push(protocolId);
                        } else {
                            selectedProtocols = selectedProtocols.filter(id => id !== protocolId);
                        }
                        
                        // Enable/disable compare button
                        $compareButton.prop('disabled', selectedProtocols.length < 2);
                    });
                    
                    // Add event listeners for view buttons
                    $('.view-protocol-btn').on('click', function() {
                        const protocolId = $(this).data('protocol-id');
                        viewProtocol(protocolId);
                    });
                    
                    // Show the protocol list
                    $('#protocol-list-container').show();
                } else {
                    $protocolList.html('<p>No protocols available for this mixture.</p>');
                }
            },
            error: function(xhr, status, error) {
                console.error('Error loading protocols:', error);
                $protocolList.html('<p>Failed to load protocols. Please try again.</p>');
            }
        });
    }
    
    // Design a protocol
    function designProtocol(event) {
        event.preventDefault();
        
        // Get form values
        const mixtureId = $mixtureSelect.val();
        const targetConcentration = parseFloat($targetConcentration.val());
        const sampleType = $sampleTypeSelect.val();
        const startingTemperature = parseFloat($startingTemperature.val());
        const targetTemperature = $targetTemperature.val() ? parseFloat($targetTemperature.val()) : null;
        const stepCount = $stepCount.val() ? parseInt($stepCount.val()) : null;
        
        // Validate inputs
        if (!mixtureId) {
            showAlert('error', 'Please select a mixture.');
            return;
        }
        
        if (isNaN(targetConcentration) || targetConcentration <= 0) {
            showAlert('error', 'Please enter a valid target concentration.');
            return;
        }
        
        // Show loading indicator
        showLoading();
        
        // Send request to API
        $.ajax({
            url: `/api/v1/protocols/design/${mixtureId}`,
            method: 'POST',
            contentType: 'application/json',
            data: JSON.stringify({
                target_concentration: targetConcentration,
                sample_type: sampleType,
                starting_temperature: startingTemperature,
                target_temperature: targetTemperature,
                step_count: stepCount
            }),
            success: function(response) {
                // Hide loading indicator
                hideLoading();
                
                // Store current protocol
                datasets: [{
                    label: 'Concentration (%)',
                    data: data,
                    borderColor: 'rgba(75, 192, 192, 1)',
                    backgroundColor: 'rgba(75, 192, 192, 0.2)',
                    borderWidth: 2,
                    fill: true
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
                            text: 'Concentration (%)'
                        }
                    },
                    x: {
                        title: {
                            display: true,
                            text: 'Protocol Steps'
                        }
                    }
                }
            }
        });
        
        $container.append($chart);
    }
    
    // Create temperature profile chart
    function createTemperatureChart($container, protocol) {
        const $chart = $('<div class="chart-section mb-4"></div>');
        $chart.append('<h6>Temperature Profile</h6>');
        
        const $chartContainer = $('<div class="chart-container" style="height: 200px;"></div>');
        $chart.append($chartContainer);
        
        // Extract temperature data
        const labels = [];
        const data = [];
        
        protocol.steps.forEach(function(step, index) {
            const temperature = step.temperature !== undefined ? step.temperature : null;
            
            labels.push(`Step ${index + 1}`);
            data.push(temperature);
        });
        
        // Create chart using Chart.js
        const ctx = $('<canvas></canvas>');
        $chartContainer.append(ctx);
        
        new Chart(ctx, {
            type: 'line',
            data: {
                labels: labels,
                datasets: [{
                    label: 'Temperature (°C)',
                    data: data,
                    borderColor: 'rgba(255, 99, 132, 1)',
                    backgroundColor: 'rgba(255, 99, 132, 0.2)',
                    borderWidth: 2,
                    fill: true
                }]
            },
            options: {
                responsive: true,
                maintainAspectRatio: false,
                scales: {
                    y: {
                        title: {
                            display: true,
                            text: 'Temperature (°C)'
                        }
                    },
                    x: {
                        title: {
                            display: true,
                            text: 'Protocol Steps'
                        }
                    }
                }
            }
        });
        
        $container.append($chart);
    }
    
    // View a specific protocol
    function viewProtocol(protocolId) {
        if (!protocolId) {
            console.error('No protocol ID provided');
            return;
        }
        
        // Show loading indicator
        showLoading();
        
        // Fetch protocol details
        $.ajax({
            url: `/api/v1/protocols/${protocolId}`,
            method: 'GET',
            success: function(response) {
                // Hide loading indicator
                hideLoading();
                
                // Display protocol
                displayProtocol(response);
                
                // Scroll to results
                $('html, body').animate({
                    scrollTop: $protocolResults.offset().top - 20
                }, 500);
            },
            error: function(xhr, status, error) {
                // Hide loading indicator
                hideLoading();
                
                console.error('Error fetching protocol:', error);
                showAlert('error', 'Failed to load protocol. Please try again.');
            }
        });
    }
    
    // Save a protocol
    function saveProtocol(event) {
        event.preventDefault();
        
        if (!currentProtocol) {
            showAlert('error', 'No protocol to save. Please design a protocol first.');
            return;
        }
        
        // Show loading indicator
        showLoading();
        
        // Send request to API
        $.ajax({
            url: '/api/v1/protocols',
            method: 'POST',
            contentType: 'application/json',
            data: JSON.stringify(currentProtocol),
            success: function(response) {
                // Hide loading indicator
                hideLoading();
                
                // Show success message
                showAlert('success', 'Protocol saved successfully.');
                
                // Reload protocols for the current mixture
                const mixtureId = $mixtureSelect.val();
                if (mixtureId) {
                    loadProtocols(mixtureId);
                }
            },
            error: function(xhr, status, error) {
                // Hide loading indicator
                hideLoading();
                
                console.error('Error saving protocol:', error);
                showAlert('error', 'Failed to save protocol. Please try again.');
            }
        });
    }
    
    // Compare selected protocols
    function compareProtocols(event) {
        event.preventDefault();
        
        if (selectedProtocols.length < 2) {
            showAlert('error', 'Please select at least two protocols to compare.');
            return;
        }
        
        // Show loading indicator
        showLoading();
        
        // Send request to API
        $.ajax({
            url: '/api/v1/protocols/compare',
            method: 'POST',
            contentType: 'application/json',
            data: JSON.stringify({
                protocol_ids: selectedProtocols
            }),
            success: function(response) {
                // Hide loading indicator
                hideLoading();
                
                // Display comparison
                displayComparison(response);
            },
            error: function(xhr, status, error) {
                // Hide loading indicator
                hideLoading();
                
                console.error('Error comparing protocols:', error);
                showAlert('error', 'Failed to compare protocols. Please try again.');
            }
        });
    }
    
    // Display protocol comparison
    function displayComparison(comparison) {
        // Create modal for comparison
        const $modal = $(`
            <div class="modal fade" id="protocol-comparison-modal" tabindex="-1" role="dialog" aria-labelledby="protocol-comparison-modal-label" aria-hidden="true">
                <div class="modal-dialog modal-xl" role="document">
                    <div class="modal-content">
                        <div class="modal-header">
                            <h5 class="modal-title" id="protocol-comparison-modal-label">Protocol Comparison</h5>
                            <button type="button" class="close" data-dismiss="modal" aria-label="Close">
                                <span aria-hidden="true">&times;</span>
                            </button>
                        </div>
                        <div class="modal-body">
                            <ul class="nav nav-tabs" id="comparison-tabs" role="tablist">
                                <li class="nav-item">
                                    <a class="nav-link active" id="summary-tab" data-toggle="tab" href="#summary" role="tab" aria-controls="summary" aria-selected="true">Summary</a>
                                </li>
                                <li class="nav-item">
                                    <a class="nav-link" id="parameters-tab" data-toggle="tab" href="#parameters" role="tab" aria-controls="parameters" aria-selected="false">Parameters</a>
                                </li>
                                <li class="nav-item">
                                    <a class="nav-link" id="steps-tab" data-toggle="tab" href="#steps" role="tab" aria-controls="steps" aria-selected="false">Steps</a>
                                </li>
                                <li class="nav-item">
                                    <a class="nav-link" id="visualization-tab" data-toggle="tab" href="#visualization" role="tab" aria-controls="visualization" aria-selected="false">Visualization</a>
                                </li>
                            </ul>
                            <div class="tab-content mt-3" id="comparison-tab-content">
                                <div class="tab-pane fade show active" id="summary" role="tabpanel" aria-labelledby="summary-tab">
                                    <div id="comparison-summary"></div>
                                    <div id="comparison-recommendations" class="mt-3"></div>
                                    <div id="comparison-details" class="mt-3"></div>
                                </div>
                                <div class="tab-pane fade" id="parameters" role="tabpanel" aria-labelledby="parameters-tab">
                                    <div id="comparison-parameters"></div>
                                </div>
                                <div class="tab-pane fade" id="steps" role="tabpanel" aria-labelledby="steps-tab">
                                    <div id="comparison-steps"></div>
                                </div>
                                <div class="tab-pane fade" id="visualization" role="tabpanel" aria-labelledby="visualization-tab">
                                    <div id="comparison-visualization"></div>
                                </div>
                            </div>
                        </div>
                        <div class="modal-footer">
                            <button type="button" class="btn btn-primary" id="export-comparison-btn">Export Report</button>
                            <button type="button" class="btn btn-secondary" data-dismiss="modal">Close</button>
                        </div>
                    </div>
                </div>
            </div>
        `);
        
        // Add modal to body
        $('body').append($modal);
        
        // Create visualization
        if (comparison.protocols && comparison.protocols.length > 0) {
            const $visualization = $('<div class="comparison-chart"></div>');
            
            // Get protocols
            const protocols = comparison.protocols || [];
            
            // Create chart header
            const $chartHeader = $('<div class="chart-header d-flex"></div>');
            $chartHeader.append('<div class="chart-label" style="width: 150px;"><strong>Step</strong></div>');
            
            protocols.forEach(function(protocol, index) {
                const protocolName = protocol.name || `Protocol ${index + 1}`;
                $chartHeader.append(`<div class="chart-column flex-grow-1 text-center"><strong>${protocolName}</strong></div>`);
            });
            
            $visualization.append($chartHeader);
            
            // Create chart rows for concentration
            const $concentrationSection = $('<div class="chart-section mb-3"></div>');
            $concentrationSection.append('<h6>Concentration Profile</h6>');
            
            // Get step comparison data
            const stepComparison = comparison.step_comparison || {};
            
            // Create concentration chart
            Object.keys(stepComparison).sort().forEach(function(stepKey) {
                const step = stepComparison[stepKey];
                const $row = $('<div class="chart-row d-flex"></div>');
                
                $row.append(`<div class="chart-label" style="width: 150px;">Step ${step.step_number || stepKey.replace('step_', '')}</div>`);
                
                step.data.forEach(function(stepData) {
                    if (stepData) {
                        const concentration = stepData.to_concentration;
                        const height = concentration * 5; // Scale for visualization
                        $row.append(`
                            <div class="chart-column flex-grow-1 text-center">
                                <div class="chart-bar bg-primary" style="height: ${height}px;" title="${concentration}%"></div>
                                <small>${concentration}%</small>
                            </div>
                        `);
                    } else {
                        $row.append('<div class="chart-column flex-grow-1 text-center">-</div>');
                    }
                });
                
                $concentrationSection.append($row);
            });
            
            $visualization.append($concentrationSection);
            
            // Add the visualization to the modal
            $modal.find('#comparison-visualization').html($visualization);
        }
        
        // Display summary
        if (comparison.summary) {
            const summary = comparison.summary;
            const $summary = $(`
                <div class="card">
                    <div class="card-header">
                        <h5 class="mb-0">Summary</h5>
                    </div>
                    <div class="card-body">
                        <div class="row">
                            <div class="col-md-6">
                                <p><strong>Protocols:</strong> ${summary.protocol_count}</p>
                                <p><strong>Sample Types:</strong> ${summary.sample_types.join(', ')}</p>
                            </div>
                            <div class="col-md-6">
                                <p><strong>Concentration Range:</strong> ${summary.concentration_range[0]}% - ${summary.concentration_range[1]}%</p>
                                <p><strong>Step Count Range:</strong> ${summary.step_count_range[0]} - ${summary.step_count_range[1]} steps</p>
                                ${summary.total_duration ? `<p><strong>Duration Range:</strong> ${Math.min(...summary.total_duration)} - ${Math.max(...summary.total_duration)} minutes</p>` : ''}
                            </div>
                        </div>
                    </div>
                </div>
            `);
            
            $modal.find('#comparison-summary').html($summary);
        }
        
        // Display recommendations
        if (comparison.recommendations && comparison.recommendations.length > 0) {
            const $recommendations = $(`
                <div class="card">
                    <div class="card-header">
                        <h5 class="mb-0">Recommendations</h5>
                    </div>
                    <div class="card-body">
                        <ul class="list-group list-group-flush">
                            ${comparison.recommendations.map(rec => `<li class="list-group-item">${rec}</li>`).join('')}
                        </ul>
                    </div>
                </div>
            `);
            
            $modal.find('#comparison-recommendations').html($recommendations);
        }
        
        // Display parameter comparison
        if (comparison.parameter_comparison) {
            const paramComparison = comparison.parameter_comparison;
            const $paramTable = $(`
                <div class="card">
                    <div class="card-header">
                        <h5 class="mb-0">Parameter Comparison</h5>
                    </div>
                    <div class="card-body">
                        <div class="table-responsive">
                            <table class="table table-bordered">
                                <thead>
                                    <tr>
                                        <th>Parameter</th>
                                        ${comparison.protocols.map((p, i) => `<th>${p.name || `Protocol ${i+1}`}</th>`).join('')}
                                        <th>Same?</th>
                                    </tr>
                                </thead>
                                <tbody>
                                </tbody>
                            </table>
                        </div>
                    </div>
                </div>
            `);
            
            const $tbody = $paramTable.find('tbody');
            
            // Add rows for each parameter
            const paramLabels = {
                'target_concentration': 'Target Concentration (%)',
                'sample_type': 'Sample Type',
                'starting_temperature': 'Starting Temperature (°C)',
                'target_temperature': 'Target Temperature (°C)',
                'step_count': 'Step Count'
            };
            
            Object.keys(paramComparison).forEach(function(param) {
                const comparison = paramComparison[param];
                const $row = $('<tr></tr>');
                
                $row.append(`<td>${paramLabels[param] || param}</td>`);
                
                comparison.values.forEach(function(value) {
                    $row.append(`<td>${value !== null ? value : '-'}</td>`);
                });
                
                // Add same/different indicator
                if (comparison.same) {
                    $row.append('<td class="table-success text-center"><i class="fas fa-check"></i> Same</td>');
                } else {
                    $row.append('<td class="table-warning text-center"><i class="fas fa-exclamation-triangle"></i> Different</td>');
                }
                
                $tbody.append($row);
            });
            
            $modal.find('#comparison-parameters').html($paramTable);
        }
        
        // Display step comparison
        if (comparison.step_comparison) {
            const stepComparison = comparison.step_comparison;
            const $stepTable = $(`
                <div class="card">
                    <div class="card-header">
                        <h5 class="mb-0">Step-by-Step Comparison</h5>
                    </div>
                    <div class="card-body">
                        <div class="table-responsive">
                            <table class="table table-bordered">
                                <thead>
                                    <tr>
                                        <th>Step</th>
                                        ${comparison.protocols.map((p, i) => `<th colspan="4">${p.name || `Protocol ${i+1}`}</th>`).join('')}
                                    </tr>
                                    <tr>
                                        <th></th>
                                        ${comparison.protocols.map(() => `
                                            <th>Concentration (%)</th>
                                            <th>Temperature (°C)</th>
                                            <th>Time (min)</th>
                                            <th>Instructions</th>
                                        `).join('')}
                                    </tr>
                                </thead>
                                <tbody>
                                </tbody>
                            </table>
                        </div>
                    </div>
                </div>
            `);
            
            const $tbody = $stepTable.find('tbody');
            
            // Add rows for each step
            Object.keys(stepComparison).sort().forEach(function(stepKey) {
                const step = stepComparison[stepKey];
                const $row = $('<tr></tr>');
                
                $row.append(`<td>Step ${stepKey.replace('step_', '')}</td>`);
                
                // For each protocol, add concentration, temperature, and time
                step.data.forEach(function(stepData) {
                    if (stepData) {
                        // Support both new and legacy field names for protocol steps
                        const concentration = stepData.concentration !== undefined ? stepData.concentration
                            : (stepData.to_concentration !== undefined ? stepData.to_concentration : '-');
                        const temperature = stepData.temperature !== undefined ? stepData.temperature : '-';
                        const duration = stepData.duration !== undefined ? stepData.duration
                            : (stepData.hold_time_min !== undefined ? stepData.hold_time_min : '-');
                        const instructions = stepData.instructions !== undefined ? stepData.instructions : '';

                        $row.append(`<td>${concentration}</td>`);
                        $row.append(`<td>${temperature}</td>`);
                        $row.append(`<td>${duration}</td>`);
                        // Display instructions as a tooltip
                        $row.append(`<td>${instructions ? `<span data-toggle="tooltip" title="${instructions}">ℹ️</span>` : '-'}</td>`);
                    } else {
                        $row.append('<td colspan="4" class="text-center">-</td>');
                    }
                });
                
                $tbody.append($row);
            });
            
            $modal.find('#comparison-steps').html($stepTable);
        }
        
        // Add export functionality
        $modal.find('#export-comparison-btn').on('click', function() {
            // Create a printable version of the comparison
            const $printContent = $(`
                <div class="container">
                    <h1 class="my-4">Protocol Comparison Report</h1>
                    <p class="text-muted">Generated on ${new Date().toLocaleString()}</p>
                    ${$modal.find('#comparison-summary').html() || ''}
                    ${$modal.find('#comparison-recommendations').html() || ''}
                    ${$modal.find('#comparison-parameters').html() || ''}
                    ${$modal.find('#comparison-steps').html() || ''}
                </div>
            `);
            
            // Open a new window and write the content
            const printWindow = window.open('', '_blank');
            printWindow.document.write(`
                <!DOCTYPE html>
                <html>
                <head>
                    <title>Protocol Comparison Report</title>
                    <link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.5.2/css/bootstrap.min.css">
                    <style>
                        body { padding: 20px; }
                        @media print {
                            .card { break-inside: avoid; }
                            .no-print { display: none; }
                        }
                    </style>
                </head>
                <body>
                    ${$printContent.html()}
                    <div class="mt-4 no-print">
                        <button class="btn btn-primary" onclick="window.print()">Print Report</button>
                        <button class="btn btn-secondary" onclick="window.close()">Close</button>
                    </div>
                </body>
                </html>
            `);
        });
        } else {
            $modal.find('#comparison-visualization').html('<p>No visualization available.</p>');
        }
        
        // Display details
        const $details = $modal.find('#comparison-details');
        
        if (comparison.protocols && comparison.protocols.length > 0) {
            const $table = $(
                `<table class="table table-striped mt-3">
                    <thead>
                        <tr>
                            <th>Protocol</th>
                            <th>Sample Type</th>
                            <th>Target Concentration (%)</th>
                            <th>Steps</th>
                            <th>Duration (min)</th>
                        </tr>
                    </thead>
                    <tbody></tbody>
                </table>`
            );
            
            const $tbody = $table.find('tbody');
            
            comparison.protocols.forEach(function(protocol) {
                const $row = $(
                    `<tr>
                        <td>${protocol.mixture_name}</td>
                        <td>${protocol.sample_type}</td>
                        <td>${protocol.target_concentration.toFixed(1)}</td>
                        <td>${protocol.step_count}</td>
                        <td>${protocol.total_duration}</td>
                    </tr>`
                );
                
                $tbody.append($row);
            });
            
            $details.append($table);
        } else {
            $details.html('<p>No comparison details available.</p>');
        }
        
        // Show modal
        $modal.modal('show');
        
        // Remove modal when hidden
        $modal.on('hidden.bs.modal', function() {
            $(this).remove();
        });
    }
    
    // Show loading indicator
    function showLoading() {
        $('#loading-overlay').show();
    }
    
    // Hide loading indicator
    function hideLoading() {
        $('#loading-overlay').hide();
    }
    
    // Show alert
    function showAlert(type, message) {
        const alertClass = type === 'error' ? 'alert-danger' : 'alert-success';
        const $alert = $(
            `<div class="alert ${alertClass} alert-dismissible fade show" role="alert">
                ${message}
                <button type="button" class="close" data-dismiss="alert" aria-label="Close">
                    <span aria-hidden="true">&times;</span>
                </button>
            </div>`
        );
        
        $('#alert-container').empty().append($alert);
        
        // Auto-dismiss after 5 seconds
        setTimeout(function() {
            $alert.alert('close');
        }, 5000);
    }
    
    // Initialize tooltips
    function initTooltips() {
        $('[data-toggle="tooltip"]').tooltip();
    }
    
    // Return public methods
    return {
        init: init
    };
})();

// Initialize when document is ready
$(document).ready(function() {
    ProtocolDesigner.init();
});

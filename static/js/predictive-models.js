/**
 * CryoProtect Analyzer - Predictive Models JavaScript
 * 
 * This file contains the JavaScript code for the predictive models functionality.
 */

// Initialize the page when the DOM is loaded
document.addEventListener('DOMContentLoaded', function() {
    // Initialize UI components
    initializePredictiveModelsUI();
    
    // Load available models
    loadAvailableModels();
    
    // Load mixtures for prediction
    loadMixturesForPrediction();
});

/**
 * Initialize UI components for predictive models
 */
function initializePredictiveModelsUI() {
    // Set up event listeners for the train model form
    const trainModelForm = document.getElementById('train-model-form');
    if (trainModelForm) {
        trainModelForm.addEventListener('submit', function(event) {
            event.preventDefault();
            trainModel();
        });
    }
    
    // Set up event listeners for the predict form
    const predictForm = document.getElementById('predict-form');
    if (predictForm) {
        predictForm.addEventListener('submit', function(event) {
            event.preventDefault();
            makePrediction();
        });
    }
    
    // Set up algorithm selection change handler
    const algorithmSelect = document.getElementById('algorithm');
    if (algorithmSelect) {
        algorithmSelect.addEventListener('change', function() {
            updateHyperparametersForm();
        });
    }
}

/**
 * Load available predictive models
 */
function loadAvailableModels() {
    // Show loading indicator
    const modelsContainer = document.getElementById('models-container');
    if (modelsContainer) {
        modelsContainer.innerHTML = '<div class="text-center"><div class="spinner-border" role="status"><span class="visually-hidden">Loading...</span></div></div>';
    }
    
    // Fetch available models from the API
    fetch('/api/v1/predictive-models')
        .then(response => {
            if (!response.ok) {
                throw new Error('Failed to load models');
            }
            return response.json();
        })
        .then(data => {
            displayAvailableModels(data);
        })
        .catch(error => {
            console.error('Error loading models:', error);
            if (modelsContainer) {
                modelsContainer.innerHTML = `<div class="alert alert-danger">Error loading models: ${error.message}</div>`;
            }
        });
}

/**
 * Display available predictive models
 * 
 * @param {Object} modelsData - Data about available models
 */
function displayAvailableModels(modelsData) {
    const modelsContainer = document.getElementById('models-container');
    if (!modelsContainer) return;
    
    // Clear container
    modelsContainer.innerHTML = '';
    
    // Check if there are any models
    if (Object.keys(modelsData).length === 0) {
        modelsContainer.innerHTML = '<div class="alert alert-info">No models available. Train a new model to get started.</div>';
        return;
    }
    
    // Create a card for each property
    for (const propertyName in modelsData) {
        const models = modelsData[propertyName];
        
        // Create property card
        const propertyCard = document.createElement('div');
        propertyCard.className = 'card mb-4';
        
        // Create card header
        const cardHeader = document.createElement('div');
        cardHeader.className = 'card-header';
        cardHeader.innerHTML = `<h5 class="mb-0">${propertyName}</h5>`;
        propertyCard.appendChild(cardHeader);
        
        // Create card body
        const cardBody = document.createElement('div');
        cardBody.className = 'card-body';
        
        // Create table for models
        const table = document.createElement('table');
        table.className = 'table table-striped';
        
        // Create table header
        const tableHeader = document.createElement('thead');
        tableHeader.innerHTML = `
            <tr>
                <th>Algorithm</th>
                <th>Trained Date</th>
                <th>Validation RMSE</th>
                <th>Validation R²</th>
                <th>Actions</th>
            </tr>
        `;
        table.appendChild(tableHeader);
        
        // Create table body
        const tableBody = document.createElement('tbody');
        
        // Add a row for each model
        models.forEach(model => {
            const row = document.createElement('tr');
            
            // Format metrics
            const validationRMSE = model.metrics?.validation?.rmse?.toFixed(4) || 'N/A';
            const validationR2 = model.metrics?.validation?.r2?.toFixed(4) || 'N/A';
            
            // Format date
            let trainedDate = 'N/A';
            if (model.trained_date) {
                const date = new Date(model.trained_date);
                trainedDate = date.toLocaleString();
            }
            
            row.innerHTML = `
                <td>${model.algorithm_name}</td>
                <td>${trainedDate}</td>
                <td>${validationRMSE}</td>
                <td>${validationR2}</td>
                <td>
                    <button class="btn btn-sm btn-info view-model-btn" data-property="${propertyName}" data-algorithm="${model.algorithm}">
                        <i class="bi bi-info-circle"></i> Details
                    </button>
                    <button class="btn btn-sm btn-danger delete-model-btn" data-property="${propertyName}" data-algorithm="${model.algorithm}">
                        <i class="bi bi-trash"></i> Delete
                    </button>
                </td>
            `;
            
            tableBody.appendChild(row);
        });
        
        table.appendChild(tableBody);
        cardBody.appendChild(table);
        propertyCard.appendChild(cardBody);
        
        modelsContainer.appendChild(propertyCard);
    }
    
    // Add event listeners for view and delete buttons
    const viewButtons = document.querySelectorAll('.view-model-btn');
    viewButtons.forEach(button => {
        button.addEventListener('click', function() {
            const propertyName = this.getAttribute('data-property');
            const algorithm = this.getAttribute('data-algorithm');
            viewModelDetails(propertyName, algorithm);
        });
    });
    
    const deleteButtons = document.querySelectorAll('.delete-model-btn');
    deleteButtons.forEach(button => {
        button.addEventListener('click', function() {
            const propertyName = this.getAttribute('data-property');
            const algorithm = this.getAttribute('data-algorithm');
            deleteModel(propertyName, algorithm);
        });
    });
}

/**
 * View details of a specific model
 * 
 * @param {string} propertyName - Name of the property
 * @param {string} algorithm - Algorithm of the model
 */
function viewModelDetails(propertyName, algorithm) {
    // Fetch model details
    fetch(`/api/v1/predictive-models/${encodeURIComponent(propertyName)}/${encodeURIComponent(algorithm)}`)
        .then(response => {
            if (!response.ok) {
                throw new Error('Failed to load model details');
            }
            return response.json();
        })
        .then(model => {
            // Create and show modal with model details
            const modalTitle = `${propertyName} - ${model.algorithm_name}`;
            
            // Format metrics
            const trainMetrics = model.metrics?.train || {};
            const validationMetrics = model.metrics?.validation || {};
            
            let metricsHtml = '<h5>Metrics</h5>';
            metricsHtml += '<table class="table table-sm">';
            metricsHtml += '<thead><tr><th></th><th>Training</th><th>Validation</th></tr></thead>';
            metricsHtml += '<tbody>';
            metricsHtml += `<tr><td>RMSE</td><td>${trainMetrics.rmse?.toFixed(4) || 'N/A'}</td><td>${validationMetrics.rmse?.toFixed(4) || 'N/A'}</td></tr>`;
            metricsHtml += `<tr><td>MAE</td><td>${trainMetrics.mae?.toFixed(4) || 'N/A'}</td><td>${validationMetrics.mae?.toFixed(4) || 'N/A'}</td></tr>`;
            metricsHtml += `<tr><td>R²</td><td>${trainMetrics.r2?.toFixed(4) || 'N/A'}</td><td>${validationMetrics.r2?.toFixed(4) || 'N/A'}</td></tr>`;
            metricsHtml += '</tbody></table>';
            
            // Format hyperparameters
            let hyperparametersHtml = '<h5>Hyperparameters</h5>';
            hyperparametersHtml += '<table class="table table-sm">';
            hyperparametersHtml += '<thead><tr><th>Parameter</th><th>Value</th></tr></thead>';
            hyperparametersHtml += '<tbody>';
            
            for (const param in model.hyperparameters) {
                let value = model.hyperparameters[param];
                if (typeof value === 'object') {
                    value = JSON.stringify(value);
                }
                hyperparametersHtml += `<tr><td>${param}</td><td>${value}</td></tr>`;
            }
            
            hyperparametersHtml += '</tbody></table>';
            
            // Format feature importance if available
            let featureImportanceHtml = '';
            if (model.feature_importance && Object.keys(model.feature_importance).length > 0) {
                featureImportanceHtml = '<h5>Feature Importance</h5>';
                
                // Sort features by importance
                const sortedFeatures = Object.entries(model.feature_importance)
                    .sort((a, b) => Math.abs(b[1]) - Math.abs(a[1]));
                
                featureImportanceHtml += '<table class="table table-sm">';
                featureImportanceHtml += '<thead><tr><th>Feature</th><th>Importance</th></tr></thead>';
                featureImportanceHtml += '<tbody>';
                
                sortedFeatures.forEach(([feature, importance]) => {
                    featureImportanceHtml += `<tr><td>${feature}</td><td>${importance.toFixed(4)}</td></tr>`;
                });
                
                featureImportanceHtml += '</tbody></table>';
            }
            
            const modalContent = `
                <div>
                    <p><strong>Property:</strong> ${propertyName}</p>
                    <p><strong>Algorithm:</strong> ${model.algorithm_name}</p>
                    <p><strong>Trained Date:</strong> ${new Date(model.trained_date).toLocaleString()}</p>
                    
                    ${metricsHtml}
                    ${hyperparametersHtml}
                    ${featureImportanceHtml}
                </div>
            `;
            
            // Show modal
            showModal(modalTitle, modalContent);
        })
        .catch(error => {
            console.error('Error loading model details:', error);
            showAlert('danger', `Error loading model details: ${error.message}`);
        });
}

/**
 * Delete a specific model
 * 
 * @param {string} propertyName - Name of the property
 * @param {string} algorithm - Algorithm of the model
 */
function deleteModel(propertyName, algorithm) {
    // Confirm deletion
    if (!confirm(`Are you sure you want to delete the ${algorithm} model for ${propertyName}?`)) {
        return;
    }
    
    // Delete the model
    fetch(`/api/v1/predictive-models/${encodeURIComponent(propertyName)}/${encodeURIComponent(algorithm)}`, {
        method: 'DELETE',
        headers: {
            'Content-Type': 'application/json'
        }
    })
        .then(response => {
            if (!response.ok) {
                throw new Error('Failed to delete model');
            }
            return response.json();
        })
        .then(data => {
            showAlert('success', data.message);
            loadAvailableModels();
        })
        .catch(error => {
            console.error('Error deleting model:', error);
            showAlert('danger', `Error deleting model: ${error.message}`);
        });
}

/**
 * Load mixtures for prediction
 */
function loadMixturesForPrediction() {
    // Get the mixture select element
    const mixtureSelect = document.getElementById('mixture-id');
    if (!mixtureSelect) return;
    
    // Show loading indicator
    mixtureSelect.innerHTML = '<option value="">Loading mixtures...</option>';
    
    // Fetch mixtures from the API
    fetch('/api/v1/mixtures')
        .then(response => {
            if (!response.ok) {
                throw new Error('Failed to load mixtures');
            }
            return response.json();
        })
        .then(mixtures => {
            // Clear select
            mixtureSelect.innerHTML = '<option value="">Select a mixture</option>';
            
            // Add options for each mixture
            mixtures.forEach(mixture => {
                const option = document.createElement('option');
                option.value = mixture.id;
                option.textContent = mixture.name;
                mixtureSelect.appendChild(option);
            });
        })
        .catch(error => {
            console.error('Error loading mixtures:', error);
            mixtureSelect.innerHTML = '<option value="">Error loading mixtures</option>';
        });
}

/**
 * Update hyperparameters form based on selected algorithm
 */
function updateHyperparametersForm() {
    const algorithm = document.getElementById('algorithm').value;
    const hyperparametersContainer = document.getElementById('hyperparameters-container');
    
    if (!hyperparametersContainer) return;
    
    // Clear container
    hyperparametersContainer.innerHTML = '';
    
    // Define hyperparameters for each algorithm
    const hyperparameters = {
        'linear_regression': [],
        'ridge_regression': [
            { name: 'alpha', label: 'Alpha (Regularization Strength)', type: 'number', default: 1.0, min: 0.0, step: 0.1 }
        ],
        'lasso_regression': [
            { name: 'alpha', label: 'Alpha (Regularization Strength)', type: 'number', default: 1.0, min: 0.0, step: 0.1 }
        ],
        'random_forest': [
            { name: 'n_estimators', label: 'Number of Trees', type: 'number', default: 100, min: 1, step: 10 },
            { name: 'max_depth', label: 'Maximum Depth', type: 'number', default: 10, min: 1, step: 1 }
        ],
        'gradient_boosting': [
            { name: 'n_estimators', label: 'Number of Boosting Stages', type: 'number', default: 100, min: 1, step: 10 },
            { name: 'learning_rate', label: 'Learning Rate', type: 'number', default: 0.1, min: 0.001, max: 1.0, step: 0.01 },
            { name: 'max_depth', label: 'Maximum Depth', type: 'number', default: 3, min: 1, step: 1 }
        ],
        'neural_network': [
            { name: 'hidden_layer_sizes', label: 'Hidden Layer Sizes (comma-separated)', type: 'text', default: '100' },
            { name: 'max_iter', label: 'Maximum Iterations', type: 'number', default: 1000, min: 100, step: 100 }
        ]
    };
    
    // Get hyperparameters for the selected algorithm
    const params = hyperparameters[algorithm] || [];
    
    // Create form fields for each hyperparameter
    if (params.length > 0) {
        const heading = document.createElement('h5');
        heading.textContent = 'Hyperparameters';
        hyperparametersContainer.appendChild(heading);
        
        params.forEach(param => {
            const formGroup = document.createElement('div');
            formGroup.className = 'mb-3';
            
            const label = document.createElement('label');
            label.htmlFor = `param-${param.name}`;
            label.className = 'form-label';
            label.textContent = param.label;
            
            const input = document.createElement('input');
            input.type = param.type;
            input.className = 'form-control';
            input.id = `param-${param.name}`;
            input.name = `param-${param.name}`;
            input.value = param.default;
            
            if (param.type === 'number') {
                if (param.min !== undefined) input.min = param.min;
                if (param.max !== undefined) input.max = param.max;
                if (param.step !== undefined) input.step = param.step;
            }
            
            formGroup.appendChild(label);
            formGroup.appendChild(input);
            
            hyperparametersContainer.appendChild(formGroup);
        });
    }
}

/**
 * Train a new predictive model
 */
function trainModel() {
    // Get form data
    const propertyName = document.getElementById('property-name').value;
    const algorithm = document.getElementById('algorithm').value;
    
    if (!propertyName) {
        showAlert('danger', 'Please enter a property name');
        return;
    }
    
    // Get hyperparameters
    const hyperparameters = {};
    const hyperparameterInputs = document.querySelectorAll('[id^="param-"]');
    
    hyperparameterInputs.forEach(input => {
        const paramName = input.id.replace('param-', '');
        let value = input.value;
        
        // Convert to appropriate type
        if (input.type === 'number') {
            value = parseFloat(value);
        } else if (paramName === 'hidden_layer_sizes') {
            // Parse hidden layer sizes as an array of integers
            value = value.split(',').map(size => parseInt(size.trim()));
        }
        
        hyperparameters[paramName] = value;
    });
    
    // Show loading indicator
    const trainButton = document.querySelector('#train-model-form button[type="submit"]');
    const originalButtonText = trainButton.innerHTML;
    trainButton.disabled = true;
    trainButton.innerHTML = '<span class="spinner-border spinner-border-sm" role="status" aria-hidden="true"></span> Training...';
    
    // Send request to train model
    fetch('/api/v1/predictive-models', {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json'
        },
        body: JSON.stringify({
            property_name: propertyName,
            algorithm: algorithm,
            hyperparameters: hyperparameters
        })
    })
        .then(response => {
            if (!response.ok) {
                throw new Error('Failed to train model');
            }
            return response.json();
        })
        .then(data => {
            showAlert('success', `Model trained successfully with validation RMSE: ${data.metrics?.validation?.rmse?.toFixed(4) || 'N/A'}`);
            loadAvailableModels();
            
            // Reset form
            document.getElementById('train-model-form').reset();
            updateHyperparametersForm();
        })
        .catch(error => {
            console.error('Error training model:', error);
            showAlert('danger', `Error training model: ${error.message}`);
        })
        .finally(() => {
            // Restore button
            trainButton.disabled = false;
            trainButton.innerHTML = originalButtonText;
        });
}

/**
 * Make a prediction for a mixture
 */
function makePrediction() {
    // Get form data
    const mixtureId = document.getElementById('mixture-id').value;
    const algorithm = document.getElementById('prediction-algorithm').value;
    
    if (!mixtureId) {
        showAlert('danger', 'Please select a mixture');
        return;
    }
    
    // Show loading indicator
    const predictButton = document.querySelector('#predict-form button[type="submit"]');
    const originalButtonText = predictButton.innerHTML;
    predictButton.disabled = true;
    predictButton.innerHTML = '<span class="spinner-border spinner-border-sm" role="status" aria-hidden="true"></span> Predicting...';
    
    // Clear previous results
    const resultsContainer = document.getElementById('prediction-results');
    if (resultsContainer) {
        resultsContainer.innerHTML = '<div class="text-center"><div class="spinner-border" role="status"><span class="visually-hidden">Loading...</span></div></div>';
    }
    
    // Send request to make prediction
    fetch(`/api/v1/mixtures/${encodeURIComponent(mixtureId)}/predict?algorithm=${encodeURIComponent(algorithm)}`)
        .then(response => {
            if (!response.ok) {
                throw new Error('Failed to make prediction');
            }
            return response.json();
        })
        .then(data => {
            displayPredictionResults(data);
        })
        .catch(error => {
            console.error('Error making prediction:', error);
            if (resultsContainer) {
                resultsContainer.innerHTML = `<div class="alert alert-danger">Error making prediction: ${error.message}</div>`;
            }
        })
        .finally(() => {
            // Restore button
            predictButton.disabled = false;
            predictButton.innerHTML = originalButtonText;
        });
}

/**
 * Display prediction results
 * 
 * @param {Object} data - Prediction results
 */
function displayPredictionResults(data) {
    const resultsContainer = document.getElementById('prediction-results');
    if (!resultsContainer) return;
    
    // Format confidence interval
    const confidenceInterval = data.confidence_interval || [0, 0];
    const confidenceIntervalStr = `${confidenceInterval[0].toFixed(2)} - ${confidenceInterval[1].toFixed(2)}`;
    
    // Create results card
    resultsContainer.innerHTML = `
        <div class="card">
            <div class="card-header">
                <h5 class="mb-0">Prediction Results</h5>
            </div>
            <div class="card-body">
                <div class="row">
                    <div class="col-md-6">
                        <p><strong>Property:</strong> ${data.property_name}</p>
                        <p><strong>Algorithm:</strong> ${data.algorithm_name}</p>
                        <p><strong>Prediction:</strong> ${data.prediction.toFixed(2)}</p>
                        <p><strong>Confidence:</strong> ${(data.confidence).toFixed(2)}</p>
                        <p><strong>Confidence Interval:</strong> ${confidenceIntervalStr}</p>
                    </div>
                    <div class="col-md-6">
                        <div class="prediction-gauge" id="prediction-gauge"></div>
                    </div>
                </div>
            </div>
        </div>
    `;
    
    // Create gauge chart
    createGaugeChart('prediction-gauge', data.prediction);
}

/**
 * Create a gauge chart for the prediction
 * 
 * @param {string} elementId - ID of the element to create the chart in
 * @param {number} value - Value to display in the gauge
 */
function createGaugeChart(elementId, value) {
    // Check if Chart.js is available
    if (typeof Chart === 'undefined') {
        console.error('Chart.js is not loaded');
        return;
    }
    
    const ctx = document.getElementById(elementId);
    if (!ctx) return;
    
    // Create gauge chart
    new Chart(ctx, {
        type: 'doughnut',
        data: {
            datasets: [{
                data: [value, 100 - value],
                backgroundColor: [
                    getColorForScore(value),
                    '#f0f0f0'
                ],
                borderWidth: 0
            }]
        },
        options: {
            circumference: 180,
            rotation: -90,
            cutout: '70%',
            plugins: {
                tooltip: {
                    enabled: false
                },
                legend: {
                    display: false
                },
                title: {
                    display: true,
                    text: 'Cryoprotection Score',
                    position: 'bottom',
                    font: {
                        size: 16
                    }
                }
            },
            layout: {
                padding: 20
            }
        },
        plugins: [{
            id: 'centerText',
            afterDraw: function(chart) {
                const width = chart.width;
                const height = chart.height;
                const ctx = chart.ctx;
                
                ctx.restore();
                ctx.font = 'bold 24px Arial';
                ctx.textBaseline = 'middle';
                ctx.textAlign = 'center';
                
                const text = `${Math.round(value)}`;
                const textX = width / 2;
                const textY = height - (height / 4);
                
                ctx.fillText(text, textX, textY);
                ctx.save();
            }
        }]
    });
}

/**
 * Get color for a score value
 * 
 * @param {number} score - Score value
 * @returns {string} - Color in hex format
 */
function getColorForScore(score) {
    if (score >= 80) {
        return '#28a745'; // Green
    } else if (score >= 60) {
        return '#17a2b8'; // Teal
    } else if (score >= 40) {
        return '#ffc107'; // Yellow
    } else if (score >= 20) {
        return '#fd7e14'; // Orange
    } else {
        return '#dc3545'; // Red
    }
}

/**
 * Show a modal with the given title and content
 * 
 * @param {string} title - Modal title
 * @param {string} content - Modal content (HTML)
 */
function showModal(title, content) {
    // Create modal elements
    const modalContainer = document.createElement('div');
    modalContainer.className = 'modal fade';
    modalContainer.id = 'dynamicModal';
    modalContainer.tabIndex = '-1';
    modalContainer.setAttribute('aria-labelledby', 'dynamicModalLabel');
    modalContainer.setAttribute('aria-hidden', 'true');
    
    modalContainer.innerHTML = `
        <div class="modal-dialog modal-lg">
            <div class="modal-content">
                <div class="modal-header">
                    <h5 class="modal-title" id="dynamicModalLabel">${title}</h5>
                    <button type="button" class="btn-close" data-bs-dismiss="modal" aria-label="Close"></button>
                </div>
                <div class="modal-body">
                    ${content}
                </div>
                <div class="modal-footer">
                    <button type="button" class="btn btn-secondary" data-bs-dismiss="modal">Close</button>
                </div>
            </div>
        </div>
    `;
    
    // Add modal to the document
    document.body.appendChild(modalContainer);
    
    // Create and show the modal
    const modal = new bootstrap.Modal(modalContainer);
    modal.show();
    
    // Remove modal from DOM when hidden
    modalContainer.addEventListener('hidden.bs.modal', function() {
        document.body.removeChild(modalContainer);
    });
}

/**
 * Show an alert message
 * 
 * @param {string} type - Alert type (success, danger, warning, info)
 * @param {string} message - Alert message
 */
function showAlert(type, message) {
    const alertsContainer = document.getElementById('alerts-container');
    if (!alertsContainer) return;
    
    // Create alert element
    const alert = document.createElement('div');
    alert.className = `alert alert-${type} alert-dismissible fade show`;
    alert.role = 'alert';
    
    alert.innerHTML = `
        ${message}
        <button type="button" class="btn-close" data-bs-dismiss="alert" aria-label="Close"></button>
    `;
    
    // Add alert to the container
    alertsContainer.appendChild(alert);
    
    // Auto-dismiss after 5 seconds
    setTimeout(() => {
        if (alert.parentNode === alertsContainer) {
            const bsAlert = new bootstrap.Alert(alert);
            bsAlert.close();
        }
    }, 5000);
}
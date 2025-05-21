/**
 * Monitoring Dashboard JS
 * 
 * This script handles the functionality for the Monitoring Dashboard page,
 * including real-time updates of system status, metrics, and logs.
 */

// Initialize dashboard when the DOM is loaded
document.addEventListener('DOMContentLoaded', function() {
    // Initialize charts
    initializeResourcesChart();
    initializeTrafficChart();
    
    // Load initial data
    loadSystemStatus();
    loadDatabaseHealth();
    loadSystemLogs();
    loadEndpointPerformance();
    
    // Set up refresh intervals
    setInterval(loadSystemStatus, 60000); // Refresh status every minute
    setInterval(loadSystemLogs, 300000); // Refresh logs every 5 minutes
    
    // Set up event listeners
    setupEventListeners();
});

/**
 * Set up event listeners for dashboard controls
 */
function setupEventListeners() {
    // Resource time range change
    document.getElementById('resourceTimeRange').addEventListener('change', function() {
        loadResourceMetrics(this.value);
    });
    
    // Traffic time range change
    document.getElementById('trafficTimeRange').addEventListener('change', function() {
        loadTrafficMetrics(this.value);
    });
    
    // Database health refresh button
    document.getElementById('refreshDbHealthBtn').addEventListener('click', function() {
        loadDatabaseHealth();
    });
    
    // Run database update button
    document.getElementById('runDbUpdateBtn').addEventListener('click', function() {
        runDatabaseUpdate();
    });
    
    // Logs refresh button
    document.getElementById('refreshLogsBtn').addEventListener('click', function() {
        loadSystemLogs();
    });
    
    // Log level filter change
    document.getElementById('logLevelFilter').addEventListener('change', function() {
        loadSystemLogs();
    });
    
    // Load more logs button
    document.getElementById('loadMoreLogsBtn').addEventListener('click', function() {
        loadMoreLogs();
    });
    
    // Endpoint performance refresh button
    document.getElementById('refreshEndpointBtn').addEventListener('click', function() {
        loadEndpointPerformance();
    });
}

/**
 * Initialize system resources chart
 */
function initializeResourcesChart() {
    const ctx = document.getElementById('resourcesChart').getContext('2d');
    
    window.resourcesChart = new Chart(ctx, {
        type: 'line',
        data: {
            labels: [], // Will be populated with time labels
            datasets: [
                {
                    label: 'CPU Usage (%)',
                    data: [],
                    borderColor: 'rgba(54, 162, 235, 1)',
                    backgroundColor: 'rgba(54, 162, 235, 0.1)',
                    borderWidth: 2,
                    tension: 0.4,
                    fill: true
                },
                {
                    label: 'Memory Usage (%)',
                    data: [],
                    borderColor: 'rgba(255, 99, 132, 1)',
                    backgroundColor: 'rgba(255, 99, 132, 0.1)',
                    borderWidth: 2,
                    tension: 0.4,
                    fill: true
                }
            ]
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            scales: {
                y: {
                    beginAtZero: true,
                    max: 100,
                    title: {
                        display: true,
                        text: 'Usage (%)'
                    }
                },
                x: {
                    title: {
                        display: true,
                        text: 'Time'
                    }
                }
            },
            interaction: {
                mode: 'index',
                intersect: false
            },
            plugins: {
                tooltip: {
                    enabled: true
                },
                legend: {
                    position: 'top'
                }
            }
        }
    });
    
    // Load initial data with 'hour' timeframe
    loadResourceMetrics('hour');
}

/**
 * Initialize traffic chart
 */
function initializeTrafficChart() {
    const ctx = document.getElementById('trafficChart').getContext('2d');
    
    window.trafficChart = new Chart(ctx, {
        type: 'line',
        data: {
            labels: [], // Will be populated with time labels
            datasets: [
                {
                    label: 'Requests',
                    data: [],
                    borderColor: 'rgba(75, 192, 192, 1)',
                    backgroundColor: 'rgba(75, 192, 192, 0.1)',
                    borderWidth: 2,
                    tension: 0.4,
                    fill: true,
                    yAxisID: 'y'
                },
                {
                    label: 'Response Time (ms)',
                    data: [],
                    borderColor: 'rgba(153, 102, 255, 1)',
                    backgroundColor: 'rgba(153, 102, 255, 0.1)',
                    borderWidth: 2,
                    tension: 0.4,
                    fill: false,
                    yAxisID: 'y1'
                },
                {
                    label: 'Error Rate (%)',
                    data: [],
                    borderColor: 'rgba(255, 159, 64, 1)',
                    backgroundColor: 'rgba(255, 159, 64, 0.1)',
                    borderWidth: 2,
                    tension: 0.4,
                    fill: false,
                    yAxisID: 'y2'
                }
            ]
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            scales: {
                y: {
                    type: 'linear',
                    display: true,
                    position: 'left',
                    beginAtZero: true,
                    title: {
                        display: true,
                        text: 'Requests'
                    }
                },
                y1: {
                    type: 'linear',
                    display: true,
                    position: 'right',
                    beginAtZero: true,
                    grid: {
                        drawOnChartArea: false
                    },
                    title: {
                        display: true,
                        text: 'Response Time (ms)'
                    }
                },
                y2: {
                    type: 'linear',
                    display: true,
                    position: 'right',
                    beginAtZero: true,
                    max: 10, // Assuming max error rate of 10%
                    grid: {
                        drawOnChartArea: false
                    },
                    title: {
                        display: true,
                        text: 'Error Rate (%)'
                    }
                },
                x: {
                    title: {
                        display: true,
                        text: 'Time'
                    }
                }
            },
            interaction: {
                mode: 'index',
                intersect: false
            },
            plugins: {
                tooltip: {
                    enabled: true
                },
                legend: {
                    position: 'top'
                }
            }
        }
    });
    
    // Load initial data with 'hour' timeframe
    loadTrafficMetrics('hour');
}

/**
 * Load system status from API
 */
function loadSystemStatus() {
    fetch('/api/v1/system/status')
        .then(response => response.json())
        .then(data => {
            // Update API status card
            updateStatusCard('apiStatusCard', data.status, data.uptime_formatted);
            
            // Update database status
            const dbStatus = data.database_status === 'connected' ? 'healthy' : 'critical';
            updateStatusCard('dbStatusCard', dbStatus, 'Database: ' + data.database_status);
            
            // Update metrics
            updateResourceMetrics(data.cpu_usage, data.memory_usage.percent);
        })
        .catch(error => {
            console.error('Error loading system status:', error);
            // Mark API as critical if we can't even fetch status
            updateStatusCard('apiStatusCard', 'critical', 'Unable to connect to API');
        });
    
    // In a real implementation, we would fetch more detailed status info
    // For this demo, we'll simulate response time and security status
    simulatePerformanceStatus();
    simulateSecurityStatus();
}

/**
 * Update status card with specified status
 * @param {string} cardId - ID of the card to update
 * @param {string} status - Status (healthy, warning, critical)
 * @param {string} statusText - Text to display
 */
function updateStatusCard(cardId, status, statusText) {
    const card = document.getElementById(cardId);
    if (!card) return;
    
    // Update status icon color
    const icon = card.querySelector('.bi-circle-fill');
    if (icon) {
        // Remove existing colors
        icon.classList.remove('text-success', 'text-warning', 'text-danger', 'text-secondary');
        
        // Add appropriate color
        if (status === 'healthy' || status === 'ok') {
            icon.classList.add('text-success');
        } else if (status === 'warning' || status === 'degraded') {
            icon.classList.add('text-warning');
        } else if (status === 'critical' || status === 'error') {
            icon.classList.add('text-danger');
        } else {
            icon.classList.add('text-secondary');
        }
    }
    
    // Update status text
    const statusTextElement = card.querySelector('.status-text');
    if (statusTextElement) {
        if (status === 'healthy' || status === 'ok') {
            statusTextElement.textContent = 'Healthy';
        } else if (status === 'warning' || status === 'degraded') {
            statusTextElement.textContent = 'Warning';
        } else if (status === 'critical' || status === 'error') {
            statusTextElement.textContent = 'Critical';
        } else {
            statusTextElement.textContent = 'Unknown';
        }
    }
    
    // Update uptime or additional text
    const uptimeText = card.querySelector('.uptime-text, .avg-response-time, .last-audit-text');
    if (uptimeText && statusText) {
        uptimeText.textContent = statusText;
    }
}

/**
 * Update resource metrics display
 * @param {number} cpuUsage - CPU usage percentage
 * @param {number} memoryUsage - Memory usage percentage
 */
function updateResourceMetrics(cpuUsage, memoryUsage) {
    // Update CPU usage
    const cpuValue = document.getElementById('cpuUsageValue');
    const cpuBar = document.getElementById('cpuUsageBar');
    
    if (cpuValue && cpuBar) {
        cpuValue.textContent = `${cpuUsage.toFixed(1)}%`;
        cpuBar.style.width = `${cpuUsage}%`;
        cpuBar.setAttribute('aria-valuenow', cpuUsage);
        
        // Update color based on usage
        cpuBar.className = 'progress-bar';
        if (cpuUsage > 90) {
            cpuBar.classList.add('bg-danger');
        } else if (cpuUsage > 70) {
            cpuBar.classList.add('bg-warning');
        } else {
            cpuBar.classList.add('bg-success');
        }
    }
    
    // Update memory usage
    const memoryValue = document.getElementById('memoryUsageValue');
    const memoryBar = document.getElementById('memoryUsageBar');
    
    if (memoryValue && memoryBar) {
        memoryValue.textContent = `${memoryUsage.toFixed(1)}%`;
        memoryBar.style.width = `${memoryUsage}%`;
        memoryBar.setAttribute('aria-valuenow', memoryUsage);
        
        // Update color based on usage
        memoryBar.className = 'progress-bar';
        if (memoryUsage > 90) {
            memoryBar.classList.add('bg-danger');
        } else if (memoryUsage > 70) {
            memoryBar.classList.add('bg-warning');
        } else {
            memoryBar.classList.add('bg-success');
        }
    }
}

/**
 * Simulate performance status
 */
function simulatePerformanceStatus() {
    // Generate random response time between 150-350ms
    const responseTime = Math.floor(Math.random() * 200) + 150;
    
    // Determine status based on response time
    let status;
    if (responseTime < 200) {
        status = 'healthy';
    } else if (responseTime < 300) {
        status = 'warning';
    } else {
        status = 'critical';
    }
    
    // Update performance status card
    updateStatusCard('performanceStatusCard', status, `Avg. Response: ${responseTime}ms`);
    
    // Update response time value
    const responseTimeValue = document.getElementById('responseTimeValue');
    if (responseTimeValue) {
        responseTimeValue.textContent = `${responseTime} ms`;
    }
}

/**
 * Simulate security status
 */
function simulateSecurityStatus() {
    // Generate a random date within the last 7 days
    const now = new Date();
    const daysAgo = Math.floor(Math.random() * 7);
    const lastAudit = new Date(now.getTime() - (daysAgo * 24 * 60 * 60 * 1000));
    const lastAuditText = `Last audit: ${lastAudit.toLocaleDateString()}`;
    
    // Determine status based on days since last audit
    let status;
    if (daysAgo < 2) {
        status = 'healthy';
    } else if (daysAgo < 5) {
        status = 'warning';
    } else {
        status = 'critical';
    }
    
    // Update security status card
    updateStatusCard('securityStatusCard', status, lastAuditText);
}

/**
 * Load resource metrics for the specified timeframe
 * @param {string} timeframe - Timeframe (hour, day, week)
 */
function loadResourceMetrics(timeframe) {
    fetch(`/api/v1/system/metrics?period=${timeframe}&metrics=cpu,memory`)
        .then(response => response.json())
        .then(data => {
            if (data && data.metrics) {
                // Update resources chart
                updateResourcesChart(data.metrics, timeframe);
            }
        })
        .catch(error => {
            console.error('Error loading resource metrics:', error);
        });
}

/**
 * Update resources chart with new data
 * @param {object} metrics - Metrics data from API
 * @param {string} timeframe - Timeframe (hour, day, week)
 */
function updateResourcesChart(metrics, timeframe) {
    if (!window.resourcesChart) return;
    
    // Get data for CPU and memory
    const cpuData = metrics.cpu || { labels: [], values: [] };
    const memoryData = metrics.memory || { labels: [], values: [] };
    
    // Format labels based on timeframe
    let labels = [];
    if (cpuData.labels && cpuData.labels.length) {
        // In a real implementation, these would be proper timestamps
        // For this demo, we'll format them based on timeframe
        if (timeframe === 'hour') {
            labels = cpuData.labels.map(l => `${l}m`);
        } else if (timeframe === 'day') {
            labels = cpuData.labels.map(l => `${l}h`);
        } else {
            labels = cpuData.labels.map(l => `Day ${l}`);
        }
    }
    
    // Update chart data
    window.resourcesChart.data.labels = labels;
    window.resourcesChart.data.datasets[0].data = cpuData.values || [];
    window.resourcesChart.data.datasets[1].data = memoryData.values || [];
    
    // Update chart
    window.resourcesChart.update();
}

/**
 * Load traffic metrics for the specified timeframe
 * @param {string} timeframe - Timeframe (hour, day, week)
 */
function loadTrafficMetrics(timeframe) {
    fetch(`/api/v1/system/metrics?period=${timeframe}&metrics=requests,response_time,errors`)
        .then(response => response.json())
        .then(data => {
            if (data && data.metrics) {
                // Update traffic chart
                updateTrafficChart(data.metrics, timeframe);
                
                // Update metric values
                updateTrafficMetricValues(data.metrics);
            }
        })
        .catch(error => {
            console.error('Error loading traffic metrics:', error);
        });
}

/**
 * Update traffic chart with new data
 * @param {object} metrics - Metrics data from API
 * @param {string} timeframe - Timeframe (hour, day, week)
 */
function updateTrafficChart(metrics, timeframe) {
    if (!window.trafficChart) return;
    
    // Get data for requests, response time, and errors
    const requestsData = metrics.requests || { labels: [], values: [] };
    const responseTimeData = metrics.response_time || { labels: [], values: [] };
    const errorsData = metrics.errors || { labels: [], values: [] };
    
    // Format labels based on timeframe
    let labels = [];
    if (requestsData.labels && requestsData.labels.length) {
        // In a real implementation, these would be proper timestamps
        // For this demo, we'll format them based on timeframe
        if (timeframe === 'hour') {
            labels = requestsData.labels.map(l => `${l}m`);
        } else if (timeframe === 'day') {
            labels = requestsData.labels.map(l => `${l}h`);
        } else {
            labels = requestsData.labels.map(l => `Day ${l}`);
        }
    }
    
    // Calculate error rates as percentage of requests
    const errorRates = [];
    if (requestsData.values && requestsData.values.length && errorsData.values && errorsData.values.length) {
        for (let i = 0; i < requestsData.values.length; i++) {
            const requests = requestsData.values[i] || 1; // Avoid division by zero
            const errors = errorsData.values[i] || 0;
            const errorRate = (errors / requests) * 100;
            errorRates.push(errorRate);
        }
    }
    
    // Update chart data
    window.trafficChart.data.labels = labels;
    window.trafficChart.data.datasets[0].data = requestsData.values || [];
    window.trafficChart.data.datasets[1].data = responseTimeData.values || [];
    window.trafficChart.data.datasets[2].data = errorRates;
    
    // Update chart
    window.trafficChart.update();
}

/**
 * Update traffic metric values display
 * @param {object} metrics - Metrics data from API
 */
function updateTrafficMetricValues(metrics) {
    // Update request rate
    const requestRateValue = document.getElementById('requestRateValue');
    if (requestRateValue && metrics.requests && metrics.requests.values) {
        // Calculate average request rate
        const sum = metrics.requests.values.reduce((a, b) => a + b, 0);
        const avg = Math.round(sum / metrics.requests.values.length);
        requestRateValue.textContent = avg;
    }
    
    // Update response time
    const responseTimeValue = document.getElementById('responseTimeValue');
    if (responseTimeValue && metrics.response_time && metrics.response_time.values) {
        // Calculate average response time
        const sum = metrics.response_time.values.reduce((a, b) => a + b, 0);
        const avg = Math.round(sum / metrics.response_time.values.length);
        responseTimeValue.textContent = `${avg} ms`;
    }
    
    // Update error rate
    const errorRateValue = document.getElementById('errorRateValue');
    if (errorRateValue && metrics.requests && metrics.requests.values && metrics.errors && metrics.errors.values) {
        // Calculate total requests and errors
        const totalRequests = metrics.requests.values.reduce((a, b) => a + b, 0);
        const totalErrors = metrics.errors.values.reduce((a, b) => a + b, 0);
        const errorRate = totalRequests > 0 ? (totalErrors / totalRequests) * 100 : 0;
        errorRateValue.textContent = `${errorRate.toFixed(2)}%`;
    }
}

/**
 * Load database health information
 */
function loadDatabaseHealth() {
    fetch('/health/database')
        .then(response => response.json())
        .then(data => {
            updateDatabaseHealthDisplay(data);
        })
        .catch(error => {
            console.error('Error loading database health:', error);
            
            // Clear database health display
            const statusElements = [
                'schemaStatusValue', 'integrityStatusValue', 
                'dbPerformanceStatusValue', 'overallHealthStatusValue'
            ];
            
            statusElements.forEach(id => {
                const element = document.getElementById(id);
                if (element) {
                    element.textContent = 'Unknown';
                }
                
                const iconElement = document.getElementById(id.replace('Value', 'Icon'));
                if (iconElement) {
                    iconElement.className = 'bi bi-question-circle fs-2';
                }
            });
        });
    
    // Also fetch database update status
    fetch('/api/v1/system/database-update')
        .then(response => response.json())
        .then(data => {
            updateDatabaseUpdateStatus(data);
        })
        .catch(error => {
            console.error('Error loading database update status:', error);
            
            // Clear database update status
            const dbUpdateStatus = document.getElementById('dbUpdateStatus');
            if (dbUpdateStatus) {
                dbUpdateStatus.textContent = 'Unable to fetch database update status';
            }
        });
}

/**
 * Update database health display
 * @param {object} data - Database health data from API
 */
function updateDatabaseHealthDisplay(data) {
    // Update schema status
    updateHealthStatusElement('schemaStatusValue', 'schemaStatusIcon', data.schema_status);
    
    // Update integrity status
    updateHealthStatusElement('integrityStatusValue', 'integrityStatusIcon', data.integrity_status);
    
    // Update performance status
    updateHealthStatusElement('dbPerformanceStatusValue', 'dbPerformanceStatusIcon', data.performance_status);
    
    // Update overall health status
    updateHealthStatusElement('overallHealthStatusValue', 'overallHealthStatusIcon', data.overall_health);
}

/**
 * Update health status element
 * @param {string} valueId - ID of the value element
 * @param {string} iconId - ID of the icon element
 * @param {string} status - Status value
 */
function updateHealthStatusElement(valueId, iconId, status) {
    const valueElement = document.getElementById(valueId);
    const iconElement = document.getElementById(iconId);
    
    if (!valueElement || !iconElement) return;
    
    // Update value text
    if (status === 'passed' || status === 'healthy') {
        valueElement.textContent = 'Healthy';
        iconElement.className = 'bi bi-check-circle-fill fs-2 text-success';
    } else if (status === 'warning') {
        valueElement.textContent = 'Warning';
        iconElement.className = 'bi bi-exclamation-circle-fill fs-2 text-warning';
    } else if (status === 'error' || status === 'failed') {
        valueElement.textContent = 'Error';
        iconElement.className = 'bi bi-x-circle-fill fs-2 text-danger';
    } else if (status === 'not_checked') {
        valueElement.textContent = 'Not Checked';
        iconElement.className = 'bi bi-dash-circle-fill fs-2 text-secondary';
    } else {
        valueElement.textContent = 'Unknown';
        iconElement.className = 'bi bi-question-circle fs-2 text-secondary';
    }
}

/**
 * Update database update status
 * @param {object} data - Database update status data from API
 */
function updateDatabaseUpdateStatus(data) {
    const dbUpdateStatus = document.getElementById('dbUpdateStatus');
    if (!dbUpdateStatus) return;
    
    if (data.status === 'success') {
        // Format date
        const lastUpdate = new Date(data.last_update);
        const formattedDate = lastUpdate.toLocaleDateString();
        const formattedTime = lastUpdate.toLocaleTimeString();
        
        dbUpdateStatus.innerHTML = `
            <div class="text-success mb-1">
                <i class="bi bi-check-circle-fill"></i> Last update completed successfully
            </div>
            <div class="small text-muted">
                ${formattedDate} at ${formattedTime}
            </div>
        `;
    } else if (data.status === 'unknown') {
        dbUpdateStatus.innerHTML = `
            <div class="text-secondary">
                <i class="bi bi-question-circle"></i> ${data.message}
            </div>
        `;
    } else {
        dbUpdateStatus.innerHTML = `
            <div class="text-danger mb-1">
                <i class="bi bi-x-circle-fill"></i> Last update failed
            </div>
            <div class="small text-muted">
                ${data.message}
            </div>
        `;
    }
}

/**
 * Run database update job
 */
function runDatabaseUpdate() {
    // Disable button during update
    const updateBtn = document.getElementById('runDbUpdateBtn');
    if (updateBtn) {
        updateBtn.disabled = true;
        updateBtn.innerHTML = '<span class="spinner-border spinner-border-sm" role="status" aria-hidden="true"></span> Running...';
    }
    
    // Update status display
    const dbUpdateStatus = document.getElementById('dbUpdateStatus');
    if (dbUpdateStatus) {
        dbUpdateStatus.innerHTML = `
            <div class="text-primary">
                <i class="bi bi-arrow-repeat spinner-border spinner-border-sm"></i> Database update in progress...
            </div>
        `;
    }
    
    // Call API to run update
    fetch('/api/v1/system/database-update', {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json'
        }
    })
        .then(response => response.json())
        .then(data => {
            // Update status display
            if (dbUpdateStatus) {
                if (data.status === 'success') {
                    dbUpdateStatus.innerHTML = `
                        <div class="text-success mb-1">
                            <i class="bi bi-check-circle-fill"></i> Database update completed successfully
                        </div>
                        <div class="small text-muted">
                            ${new Date().toLocaleString()}
                        </div>
                    `;
                } else if (data.status === 'warning') {
                    dbUpdateStatus.innerHTML = `
                        <div class="text-warning mb-1">
                            <i class="bi bi-exclamation-circle-fill"></i> Database update completed with warnings
                        </div>
                        <div class="small text-muted">
                            ${data.message}
                        </div>
                    `;
                } else {
                    dbUpdateStatus.innerHTML = `
                        <div class="text-danger mb-1">
                            <i class="bi bi-x-circle-fill"></i> Database update failed
                        </div>
                        <div class="small text-muted">
                            ${data.message}
                        </div>
                    `;
                }
            }
            
            // Re-enable button
            if (updateBtn) {
                updateBtn.disabled = false;
                updateBtn.innerHTML = 'Run Database Update';
            }
            
            // Refresh database health after update
            loadDatabaseHealth();
        })
        .catch(error => {
            console.error('Error running database update:', error);
            
            // Update status display
            if (dbUpdateStatus) {
                dbUpdateStatus.innerHTML = `
                    <div class="text-danger mb-1">
                        <i class="bi bi-x-circle-fill"></i> Error running database update
                    </div>
                    <div class="small text-muted">
                        ${error.message}
                    </div>
                `;
            }
            
            // Re-enable button
            if (updateBtn) {
                updateBtn.disabled = false;
                updateBtn.innerHTML = 'Run Database Update';
            }
        });
}

// Variables for log pagination
let logOffset = 0;
const logLimit = 50;

/**
 * Load system logs
 */
function loadSystemLogs() {
    // Reset pagination
    logOffset = 0;
    
    // Get log level filter
    const logLevel = document.getElementById('logLevelFilter').value;
    
    // Show loading indicator
    const logContainer = document.getElementById('logContainer');
    if (logContainer) {
        logContainer.innerHTML = `
            <div class="text-center py-4">
                <div class="spinner-border text-primary" role="status">
                    <span class="visually-hidden">Loading...</span>
                </div>
                <p class="mt-2">Loading logs...</p>
            </div>
        `;
    }
    
    // Construct URL with filters
    let url = `/api/v1/system/logs?limit=${logLimit}&offset=${logOffset}`;
    if (logLevel !== 'all') {
        url += `&level=${logLevel}`;
    }
    
    // Fetch logs from API
    fetch(url)
        .then(response => response.json())
        .then(data => {
            // Display logs
            displayLogs(data.logs, data.total_count);
        })
        .catch(error => {
            console.error('Error loading logs:', error);
            
            // Show error message
            if (logContainer) {
                logContainer.innerHTML = `
                    <div class="text-center py-4 text-danger">
                        <i class="bi bi-exclamation-circle fs-1"></i>
                        <p class="mt-2">Error loading logs: ${error.message}</p>
                    </div>
                `;
            }
        });
}

/**
 * Load more logs (pagination)
 */
function loadMoreLogs() {
    // Increment offset
    logOffset += logLimit;
    
    // Get log level filter
    const logLevel = document.getElementById('logLevelFilter').value;
    
    // Disable load more button
    const loadMoreBtn = document.getElementById('loadMoreLogsBtn');
    if (loadMoreBtn) {
        loadMoreBtn.disabled = true;
        loadMoreBtn.innerHTML = '<span class="spinner-border spinner-border-sm" role="status" aria-hidden="true"></span> Loading...';
    }
    
    // Construct URL with filters
    let url = `/api/v1/system/logs?limit=${logLimit}&offset=${logOffset}`;
    if (logLevel !== 'all') {
        url += `&level=${logLevel}`;
    }
    
    // Fetch logs from API
    fetch(url)
        .then(response => response.json())
        .then(data => {
            // Append logs
            appendLogs(data.logs, data.total_count);
            
            // Re-enable load more button
            if (loadMoreBtn) {
                loadMoreBtn.disabled = false;
                loadMoreBtn.innerHTML = 'Load More';
                
                // Hide button if no more logs
                if (logOffset + logLimit >= data.total_count) {
                    loadMoreBtn.style.display = 'none';
                }
            }
        })
        .catch(error => {
            console.error('Error loading more logs:', error);
            
            // Show error message and re-enable button
            if (loadMoreBtn) {
                loadMoreBtn.disabled = false;
                loadMoreBtn.innerHTML = 'Load More';
            }
        });
}

/**
 * Display logs in the log container
 * @param {Array} logs - Array of log entries
 * @param {number} totalCount - Total count of logs
 */
function displayLogs(logs, totalCount) {
    const logContainer = document.getElementById('logContainer');
    const logCount = document.getElementById('logCount');
    
    if (!logContainer) return;
    
    // Clear existing logs
    logContainer.innerHTML = '';
    
    // Display logs
    if (logs && logs.length) {
        logs.forEach(log => {
            // Create log entry
            const logEntry = document.createElement('div');
            logEntry.className = `log-entry py-1 log-${log.level.toLowerCase()}`;
            
            // Format timestamp
            const timestamp = new Date(log.timestamp);
            const formattedTime = timestamp.toLocaleTimeString();
            
            // Create log content
            logEntry.innerHTML = `
                <span class="text-muted me-2">[${formattedTime}]</span>
                <span class="me-2">${log.level}</span>
                <span class="me-2">${log.source}</span>
                <span>${log.message}</span>
            `;
            
            // Add log entry to container
            logContainer.appendChild(logEntry);
        });
        
        // Update log count
        if (logCount) {
            logCount.textContent = totalCount;
        }
        
        // Update load more button visibility
        const loadMoreBtn = document.getElementById('loadMoreLogsBtn');
        if (loadMoreBtn) {
            if (logs.length < logLimit || logOffset + logLimit >= totalCount) {
                loadMoreBtn.style.display = 'none';
            } else {
                loadMoreBtn.style.display = 'block';
            }
        }
    } else {
        // No logs found
        logContainer.innerHTML = `
            <div class="text-center py-4 text-muted">
                <i class="bi bi-journal fs-1"></i>
                <p class="mt-2">No logs found</p>
            </div>
        `;
        
        // Update log count
        if (logCount) {
            logCount.textContent = '0';
        }
        
        // Hide load more button
        const loadMoreBtn = document.getElementById('loadMoreLogsBtn');
        if (loadMoreBtn) {
            loadMoreBtn.style.display = 'none';
        }
    }
}

/**
 * Append logs to the log container
 * @param {Array} logs - Array of log entries
 * @param {number} totalCount - Total count of logs
 */
function appendLogs(logs, totalCount) {
    const logContainer = document.getElementById('logContainer');
    const logCount = document.getElementById('logCount');
    
    if (!logContainer) return;
    
    // Display logs
    if (logs && logs.length) {
        logs.forEach(log => {
            // Create log entry
            const logEntry = document.createElement('div');
            logEntry.className = `log-entry py-1 log-${log.level.toLowerCase()}`;
            
            // Format timestamp
            const timestamp = new Date(log.timestamp);
            const formattedTime = timestamp.toLocaleTimeString();
            
            // Create log content
            logEntry.innerHTML = `
                <span class="text-muted me-2">[${formattedTime}]</span>
                <span class="me-2">${log.level}</span>
                <span class="me-2">${log.source}</span>
                <span>${log.message}</span>
            `;
            
            // Add log entry to container
            logContainer.appendChild(logEntry);
        });
        
        // Update log count
        if (logCount) {
            logCount.textContent = totalCount;
        }
    }
}

/**
 * Load endpoint performance data
 */
function loadEndpointPerformance() {
    // For this demo, we'll use simulated data
    // In a real implementation, this would fetch from an API
    simulateEndpointPerformance();
}

/**
 * Simulate endpoint performance data
 */
function simulateEndpointPerformance() {
    // Define sample endpoints
    const endpoints = [
        { path: '/api/v1/molecules', method: 'GET', requests: 1250, avg_time: 120, error_rate: 0.1 },
        { path: '/api/v1/molecules/{id}', method: 'GET', requests: 3780, avg_time: 85, error_rate: 0.2 },
        { path: '/api/v1/properties/data', method: 'GET', requests: 950, avg_time: 230, error_rate: 1.5 },
        { path: '/api/v1/mixtures', method: 'GET', requests: 420, avg_time: 180, error_rate: 0.3 },
        { path: '/api/v1/system/status', method: 'GET', requests: 320, avg_time: 45, error_rate: 0 },
        { path: '/api/v1/stats/dashboard', method: 'GET', requests: 520, avg_time: 95, error_rate: 0.2 },
        { path: '/api/v1/properties/correlation', method: 'POST', requests: 180, avg_time: 350, error_rate: 2.1 },
        { path: '/api/v1/molecules', method: 'POST', requests: 85, avg_time: 220, error_rate: 1.8 },
        { path: '/api/v1/system/database/health', method: 'GET', requests: 120, avg_time: 310, error_rate: 0.5 },
        { path: '/api/v1/system/database-update', method: 'POST', requests: 12, avg_time: 2500, error_rate: 8.3 }
    ];
    
    // Display endpoint performance
    displayEndpointPerformance(endpoints);
}

/**
 * Display endpoint performance data
 * @param {Array} endpoints - Array of endpoint performance data
 */
function displayEndpointPerformance(endpoints) {
    const tableBody = document.getElementById('endpointTableBody');
    
    if (!tableBody) return;
    
    // Clear existing rows
    tableBody.innerHTML = '';
    
    // Sort endpoints by highest error rate
    endpoints.sort((a, b) => b.error_rate - a.error_rate);
    
    // Display endpoints
    endpoints.forEach(endpoint => {
        // Create row
        const row = document.createElement('tr');
        
        // Determine status based on error rate and response time
        let status;
        let statusClass;
        
        if (endpoint.error_rate > 5 || endpoint.avg_time > 1000) {
            status = 'Critical';
            statusClass = 'text-danger';
        } else if (endpoint.error_rate > 1 || endpoint.avg_time > 300) {
            status = 'Warning';
            statusClass = 'text-warning';
        } else {
            status = 'Healthy';
            statusClass = 'text-success';
        }
        
        // Create row content
        row.innerHTML = `
            <td>${endpoint.path}</td>
            <td>${endpoint.method}</td>
            <td>${endpoint.requests.toLocaleString()}</td>
            <td>${endpoint.avg_time} ms</td>
            <td>${endpoint.error_rate.toFixed(1)}%</td>
            <td class="${statusClass}">
                <i class="bi bi-circle-fill me-1"></i>${status}
            </td>
        `;
        
        // Add row to table
        tableBody.appendChild(row);
    });
}
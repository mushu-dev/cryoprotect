/**
 * CryoProtect Analyzer - Enhanced Export and Sharing Functionality
 *
 * This module provides advanced functionality for exporting data in various formats,
 * batch exporting, visualization exports with templates, and comprehensive sharing
 * options including permissions management and analytics.
 */

// Export functionality
const ExportManager = {
    /**
     * Initialize export functionality
     */
    init: function() {
        this.bindEvents();
        this.loadSavedTemplates();
        this.setupAccessibilityFeatures();
        
        // Initialize offline support
        this.offlineQueue = this.loadOfflineQueue();
        this.setupOfflineSupport();
    },

    /**
     * Bind event listeners
     */
    bindEvents: function() {
        // Export buttons
        document.querySelectorAll('.export-btn').forEach(btn => {
            btn.addEventListener('click', this.handleExportClick.bind(this));
        });

        // Export format selection
        document.querySelectorAll('.export-format-select').forEach(select => {
            select.addEventListener('change', this.updateExportOptions.bind(this));
        });

        // Share buttons
        document.querySelectorAll('.share-btn').forEach(btn => {
            btn.addEventListener('click', this.handleShareClick.bind(this));
        });
        
        // Unified export dialog
        document.querySelectorAll('.open-unified-export-btn').forEach(btn => {
            btn.addEventListener('click', this.openUnifiedExportDialog.bind(this));
        });
        
        // Batch export button
        document.querySelectorAll('.batch-export-btn').forEach(btn => {
            btn.addEventListener('click', this.handleBatchExportClick.bind(this));
        });
        
        // Visualization preview button
        document.querySelectorAll('.preview-visualization-btn').forEach(btn => {
            btn.addEventListener('click', this.previewVisualization.bind(this));
        });
        
        // Save visualization template button
        document.querySelectorAll('.save-visualization-template-btn').forEach(btn => {
            btn.addEventListener('click', this.saveVisualizationTemplate.bind(this));
        });
        
        // Load visualization template selection
        document.querySelectorAll('.visualization-template-select').forEach(select => {
            select.addEventListener('change', this.loadVisualizationTemplate.bind(this));
        });
    },

    /**
     * Load saved templates from localStorage or server
     */
    loadSavedTemplates: function() {
        // Try to load from localStorage first (for offline support)
        const savedTemplates = localStorage.getItem('visualization-templates');
        if (savedTemplates) {
            this.visualizationTemplates = JSON.parse(savedTemplates);
            this.populateTemplateSelects();
        }
        
        // Then try to load from server (will override localStorage if successful)
        fetch('/api/v1/templates/visualization', {
            method: 'GET',
            headers: {
                'Authorization': `Bearer ${this.getAuthToken()}`
            }
        })
        .then(response => {
            if (!response.ok) {
                throw new Error('Failed to load templates');
            }
            return response.json();
        })
        .then(data => {
            this.visualizationTemplates = data.templates;
            localStorage.setItem('visualization-templates', JSON.stringify(this.visualizationTemplates));
            this.populateTemplateSelects();
        })
        .catch(error => {
            console.warn('Could not load templates from server, using cached templates:', error);
        });
    },
    
    /**
     * Populate template select dropdowns
     */
    populateTemplateSelects: function() {
        if (!this.visualizationTemplates) return;
        
        document.querySelectorAll('.visualization-template-select').forEach(select => {
            // Clear existing options except the default
            Array.from(select.options).forEach(option => {
                if (option.value !== 'none') {
                    select.removeChild(option);
                }
            });
            
            // Add template options
            this.visualizationTemplates.forEach(template => {
                const option = document.createElement('option');
                option.value = template.id;
                option.textContent = template.name;
                select.appendChild(option);
            });
        });
    },
    
    /**
     * Setup accessibility features
     */
    setupAccessibilityFeatures: function() {
        // Add ARIA attributes to export and share controls
        document.querySelectorAll('.export-btn').forEach(btn => {
            btn.setAttribute('aria-label', 'Export data');
            btn.setAttribute('role', 'button');
        });
        
        document.querySelectorAll('.share-btn').forEach(btn => {
            btn.setAttribute('aria-label', 'Share data');
            btn.setAttribute('role', 'button');
        });
        
        // Add keyboard navigation
        document.querySelectorAll('.export-format-select, .share-type-select').forEach(select => {
            select.addEventListener('keydown', (e) => {
                if (e.key === 'Enter') {
                    e.target.click();
                }
            });
        });
    },
    
    /**
     * Setup offline support
     */
    setupOfflineSupport: function() {
        // Check for online status changes
        window.addEventListener('online', this.processOfflineQueue.bind(this));
        window.addEventListener('offline', () => {
            this.showNotification('info', 'You are offline. Export and share operations will be queued for when you reconnect.');
        });
        
        // Process any pending operations if we're online
        if (navigator.onLine) {
            this.processOfflineQueue();
        }
    },
    
    /**
     * Load offline queue from localStorage
     */
    loadOfflineQueue: function() {
        const queue = localStorage.getItem('offline-export-queue');
        return queue ? JSON.parse(queue) : [];
    },
    
    /**
     * Save offline queue to localStorage
     */
    saveOfflineQueue: function() {
        localStorage.setItem('offline-export-queue', JSON.stringify(this.offlineQueue));
    },
    
    /**
     * Add operation to offline queue
     */
    addToOfflineQueue: function(operation) {
        this.offlineQueue.push(operation);
        this.saveOfflineQueue();
        this.showNotification('info', 'Operation added to offline queue');
    },
    
    /**
     * Process offline queue when back online
     */
    processOfflineQueue: function() {
        if (this.offlineQueue.length === 0) return;
        
        this.showNotification('info', `Processing ${this.offlineQueue.length} pending operations...`);
        
        const queue = [...this.offlineQueue];
        this.offlineQueue = [];
        this.saveOfflineQueue();
        
        let processed = 0;
        
        queue.forEach(operation => {
            switch (operation.type) {
                case 'export':
                    this.exportData(
                        operation.dataType,
                        operation.itemId,
                        operation.format,
                        operation.includeRelated,
                        operation.filters,
                        operation.metadata
                    );
                    break;
                case 'visualization':
                    this.exportVisualization(
                        operation.chartType,
                        operation.chartData,
                        operation.format,
                        operation.options
                    );
                    break;
                case 'report':
                    this.generateReport(
                        operation.title,
                        operation.sections,
                        operation.options
                    );
                    break;
                case 'batch':
                    this.batchExport(
                        operation.items,
                        operation.format,
                        operation.options
                    );
                    break;
            }
            processed++;
        });
        
        this.showNotification('success', `Processed ${processed} operations from offline queue`);
    },

    /**
     * Handle export button click
     * @param {Event} event - Click event
     */
    handleExportClick: function(event) {
        event.preventDefault();
        const btn = event.currentTarget;
        const dataType = btn.dataset.type;
        const itemId = btn.dataset.id;
        const format = btn.closest('.export-container').querySelector('.export-format-select').value;
        const includeRelated = btn.closest('.export-container').querySelector('.include-related-checkbox')?.checked || false;

        // Get any additional filters or metadata if available
        const filters = {};
        const metadata = {};
        
        const filterElements = btn.closest('.export-container').querySelectorAll('.filter-option select, .filter-option input');
        if (filterElements.length > 0) {
            filterElements.forEach(el => {
                if (el.type === 'checkbox') {
                    filters[el.id] = el.checked;
                } else {
                    filters[el.id] = el.value;
                }
            });
        }
        
        const metadataElements = btn.closest('.export-container').querySelectorAll('.metadata-option input, .metadata-option textarea');
        if (metadataElements.length > 0) {
            metadataElements.forEach(el => {
                if (el.type === 'checkbox') {
                    metadata[el.id] = el.checked;
                } else {
                    metadata[el.id] = el.value;
                }
            });
        }

        this.exportData(dataType, itemId, format, includeRelated, filters, metadata);
    },

    /**
     * Update export options based on selected format
     * @param {Event} event - Change event
     */
    updateExportOptions: function(event) {
        const format = event.target.value;
        const container = event.target.closest('.export-container');
        
        // Show/hide format-specific options
        container.querySelectorAll('.format-option').forEach(option => {
            option.style.display = 'none';
        });
        
        container.querySelectorAll(`.format-option.${format}`).forEach(option => {
            option.style.display = 'block';
        });
        
        // Update preview if available
        const previewBtn = container.querySelector('.preview-export-btn');
        if (previewBtn) {
            previewBtn.click();
        }
    },
    
    /**
     * Handle batch export button click
     * @param {Event} event - Click event
     */
    handleBatchExportClick: function(event) {
        event.preventDefault();
        const btn = event.currentTarget;
        const dataType = btn.dataset.type;
        const format = btn.closest('.batch-export-container').querySelector('.export-format-select').value;
        
        // Get selected items
        const selectedItems = Array.from(
            document.querySelectorAll('.batch-export-item-checkbox:checked')
        ).map(checkbox => checkbox.value);
        
        if (selectedItems.length === 0) {
            this.showNotification('error', 'Please select at least one item to export');
            return;
        }
        
        // Get batch export options
        const options = {
            includeRelated: btn.closest('.batch-export-container').querySelector('.include-related-checkbox')?.checked || false,
            combineFiles: btn.closest('.batch-export-container').querySelector('.combine-files-checkbox')?.checked || false,
            zipFiles: btn.closest('.batch-export-container').querySelector('.zip-files-checkbox')?.checked || false
        };
        
        this.batchExport(dataType, selectedItems, format, options);
    },

    /**
     * Open unified export dialog
     * @param {Event} event - Click event
     */
    openUnifiedExportDialog: function(event) {
        event.preventDefault();
        const btn = event.currentTarget;
        const dataType = btn.dataset.type;
        const itemId = btn.dataset.id;
        
        // Create dialog element
        const dialog = document.createElement('div');
        dialog.className = 'export-dialog';
        dialog.setAttribute('role', 'dialog');
        dialog.setAttribute('aria-labelledby', 'export-dialog-title');
        
        // Create dialog content
        dialog.innerHTML = `
            <div class="export-dialog-content">
                <h3 id="export-dialog-title">Export Data</h3>
                <div class="export-dialog-tabs">
                    <button class="export-tab-btn active" data-tab="format">Format</button>
                    <button class="export-tab-btn" data-tab="filters">Filters</button>
                    <button class="export-tab-btn" data-tab="metadata">Metadata</button>
                    <button class="export-tab-btn" data-tab="preview">Preview</button>
                </div>
                
                <div class="export-tab-content format active">
                    <div class="export-option">
                        <label for="unified-export-format">Export Format:</label>
                        <select id="unified-export-format" class="export-format-select">
                            <option value="csv">CSV</option>
                            <option value="json">JSON</option>
                            <option value="excel">Excel</option>
                            <option value="pdf">PDF</option>
                        </select>
                    </div>
                    
                    <div class="export-option format-option csv">
                        <label for="csv-delimiter">CSV Delimiter:</label>
                        <select id="csv-delimiter">
                            <option value="comma">Comma (,)</option>
                            <option value="tab">Tab</option>
                            <option value="semicolon">Semicolon (;)</option>
                        </select>
                    </div>
                    
                    <div class="export-option format-option excel">
                        <label for="excel-sheet-name">Sheet Name:</label>
                        <input type="text" id="excel-sheet-name" value="Data">
                    </div>
                    
                    <div class="export-option format-option pdf">
                        <label for="pdf-orientation">Orientation:</label>
                        <select id="pdf-orientation">
                            <option value="portrait">Portrait</option>
                            <option value="landscape">Landscape</option>
                        </select>
                    </div>
                    
                    <div class="include-related-container">
                        <input type="checkbox" id="unified-include-related" class="include-related-checkbox">
                        <label for="unified-include-related">Include Related Data</label>
                    </div>
                </div>
                
                <div class="export-tab-content filters">
                    <div class="filter-option">
                        <label for="filter-date-range">Date Range:</label>
                        <select id="filter-date-range">
                            <option value="all">All Time</option>
                            <option value="today">Today</option>
                            <option value="week">This Week</option>
                            <option value="month">This Month</option>
                            <option value="year">This Year</option>
                            <option value="custom">Custom Range</option>
                        </select>
                    </div>
                    
                    <div class="filter-option custom-date-range" style="display: none;">
                        <label for="filter-date-start">Start Date:</label>
                        <input type="date" id="filter-date-start">
                        
                        <label for="filter-date-end">End Date:</label>
                        <input type="date" id="filter-date-end">
                    </div>
                    
                    <div class="filter-option">
                        <label for="filter-properties">Properties to Include:</label>
                        <div class="property-checkboxes" id="filter-properties">
                            <!-- Will be populated dynamically based on data type -->
                        </div>
                    </div>
                </div>
                
                <div class="export-tab-content metadata">
                    <div class="metadata-option">
                        <label for="metadata-title">Title:</label>
                        <input type="text" id="metadata-title" placeholder="Export Title">
                    </div>
                    
                    <div class="metadata-option">
                        <label for="metadata-description">Description:</label>
                        <textarea id="metadata-description" placeholder="Export Description"></textarea>
                    </div>
                    
                    <div class="metadata-option">
                        <label for="metadata-author">Author:</label>
                        <input type="text" id="metadata-author" placeholder="Author Name">
                    </div>
                    
                    <div class="metadata-option">
                        <label for="metadata-keywords">Keywords (comma-separated):</label>
                        <input type="text" id="metadata-keywords" placeholder="keyword1, keyword2, ...">
                    </div>
                    
                    <div class="metadata-option">
                        <input type="checkbox" id="metadata-include-timestamp" checked>
                        <label for="metadata-include-timestamp">Include Timestamp</label>
                    </div>
                </div>
                
                <div class="export-tab-content preview">
                    <div class="preview-container">
                        <p>Preview will be shown here after clicking "Generate Preview"</p>
                    </div>
                    <button class="generate-preview-btn">Generate Preview</button>
                </div>
                
                <div class="export-dialog-actions">
                    <button class="cancel-export-btn">Cancel</button>
                    <button class="confirm-export-btn" data-type="${dataType}" data-id="${itemId}">Export</button>
                </div>
            </div>
        `;
        
        // Add event listeners
        dialog.querySelector('.cancel-export-btn').addEventListener('click', () => {
            dialog.remove();
        });
        
        dialog.querySelector('.confirm-export-btn').addEventListener('click', (e) => {
            const btn = e.currentTarget;
            const dataType = btn.dataset.type;
            const itemId = btn.dataset.id;
            const format = dialog.querySelector('#unified-export-format').value;
            const includeRelated = dialog.querySelector('#unified-include-related').checked;
            
            // Get format-specific options
            const formatOptions = {};
            if (format === 'csv') {
                formatOptions.delimiter = dialog.querySelector('#csv-delimiter').value;
            } else if (format === 'excel') {
                formatOptions.sheetName = dialog.querySelector('#excel-sheet-name').value;
            } else if (format === 'pdf') {
                formatOptions.orientation = dialog.querySelector('#pdf-orientation').value;
            }
            
            // Get filters
            const filters = {
                dateRange: dialog.querySelector('#filter-date-range').value
            };
            
            if (filters.dateRange === 'custom') {
                filters.startDate = dialog.querySelector('#filter-date-start').value;
                filters.endDate = dialog.querySelector('#filter-date-end').value;
            }
            
            // Get selected properties
            filters.properties = Array.from(dialog.querySelectorAll('#filter-properties input:checked')).map(cb => cb.value);
            
            // Get metadata
            const metadata = {
                title: dialog.querySelector('#metadata-title').value,
                description: dialog.querySelector('#metadata-description').value,
                author: dialog.querySelector('#metadata-author').value,
                keywords: dialog.querySelector('#metadata-keywords').value.split(',').map(k => k.trim()),
                includeTimestamp: dialog.querySelector('#metadata-include-timestamp').checked
            };
            
            // Export with all options
            this.exportData(dataType, itemId, format, includeRelated, filters, metadata, formatOptions);
            
            // Close dialog
            dialog.remove();
        });
        
        // Tab switching
        dialog.querySelectorAll('.export-tab-btn').forEach(btn => {
            btn.addEventListener('click', (e) => {
                // Update active tab button
                dialog.querySelectorAll('.export-tab-btn').forEach(b => b.classList.remove('active'));
                e.currentTarget.classList.add('active');
                
                // Show selected tab content
                const tabName = e.currentTarget.dataset.tab;
                dialog.querySelectorAll('.export-tab-content').forEach(content => {
                    content.classList.remove('active');
                });
                dialog.querySelector(`.export-tab-content.${tabName}`).classList.add('active');
            });
        });
        
        // Custom date range toggle
        dialog.querySelector('#filter-date-range').addEventListener('change', (e) => {
            const customDateRange = dialog.querySelector('.custom-date-range');
            customDateRange.style.display = e.target.value === 'custom' ? 'block' : 'none';
        });
        
        // Generate preview button
        dialog.querySelector('.generate-preview-btn').addEventListener('click', () => {
            const previewContainer = dialog.querySelector('.preview-container');
            previewContainer.innerHTML = '<p>Loading preview...</p>';
            
            const dataType = dialog.querySelector('.confirm-export-btn').dataset.type;
            const itemId = dialog.querySelector('.confirm-export-btn').dataset.id;
            const format = dialog.querySelector('#unified-export-format').value;
            
            // Get preview from API
            fetch(`/api/v1/export/preview?data_type=${dataType}&id=${itemId}&format=${format}`, {
                headers: {
                    'Authorization': `Bearer ${this.getAuthToken()}`
                }
            })
            .then(response => response.json())
            .then(data => {
                if (data.preview) {
                    if (format === 'csv' || format === 'json') {
                        previewContainer.innerHTML = `<pre>${data.preview}</pre>`;
                    } else {
                        previewContainer.innerHTML = data.preview;
                    }
                } else {
                    previewContainer.innerHTML = '<p>Preview not available for this format</p>';
                }
            })
            .catch(error => {
                previewContainer.innerHTML = `<p>Error generating preview: ${error.message}</p>`;
            });
        });
        
        // Populate property checkboxes based on data type
        this.populatePropertyCheckboxes(dialog.querySelector('#filter-properties'), dataType);
        
        // Add to document
        document.body.appendChild(dialog);
    },
    
    /**
     * Populate property checkboxes based on data type
     * @param {HTMLElement} container - Container element
     * @param {string} dataType - Data type
     */
    populatePropertyCheckboxes: function(container, dataType) {
        // Clear container
        container.innerHTML = '';
        
        // Define properties based on data type
        let properties = [];
        
        switch (dataType) {
            case 'molecules':
                properties = [
                    { id: 'name', label: 'Name', checked: true },
                    { id: 'formula', label: 'Formula', checked: true },
                    { id: 'smiles', label: 'SMILES', checked: true },
                    { id: 'molecular_weight', label: 'Molecular Weight', checked: true },
                    { id: 'logp', label: 'LogP', checked: false },
                    { id: 'h_bond_donors', label: 'H-Bond Donors', checked: false },
                    { id: 'h_bond_acceptors', label: 'H-Bond Acceptors', checked: false },
                    { id: 'rotatable_bonds', label: 'Rotatable Bonds', checked: false },
                    { id: 'polar_surface_area', label: 'Polar Surface Area', checked: false },
                    { id: 'heavy_atoms', label: 'Heavy Atoms', checked: false }
                ];
                break;
            case 'mixtures':
                properties = [
                    { id: 'name', label: 'Name', checked: true },
                    { id: 'description', label: 'Description', checked: true },
                    { id: 'components', label: 'Components', checked: true },
                    { id: 'concentrations', label: 'Concentrations', checked: true },
                    { id: 'created_at', label: 'Creation Date', checked: true },
                    { id: 'updated_at', label: 'Last Updated', checked: false },
                    { id: 'created_by', label: 'Created By', checked: false },
                    { id: 'tags', label: 'Tags', checked: false }
                ];
                break;
            case 'predictions':
                properties = [
                    { id: 'mixture_id', label: 'Mixture ID', checked: true },
                    { id: 'property', label: 'Property', checked: true },
                    { id: 'value', label: 'Value', checked: true },
                    { id: 'confidence', label: 'Confidence', checked: true },
                    { id: 'method', label: 'Method', checked: true },
                    { id: 'created_at', label: 'Creation Date', checked: false },
                    { id: 'model_version', label: 'Model Version', checked: false }
                ];
                break;
            case 'experiments':
                properties = [
                    { id: 'name', label: 'Name', checked: true },
                    { id: 'description', label: 'Description', checked: true },
                    { id: 'mixture_id', label: 'Mixture ID', checked: true },
                    { id: 'protocol', label: 'Protocol', checked: true },
                    { id: 'results', label: 'Results', checked: true },
                    { id: 'date', label: 'Date', checked: true },
                    { id: 'researcher', label: 'Researcher', checked: false },
                    { id: 'notes', label: 'Notes', checked: false },
                    { id: 'tags', label: 'Tags', checked: false }
                ];
                break;
            default:
                properties = [
                    { id: 'id', label: 'ID', checked: true },
                    { id: 'name', label: 'Name', checked: true },
                    { id: 'description', label: 'Description', checked: true },
                    { id: 'created_at', label: 'Creation Date', checked: true },
                    { id: 'updated_at', label: 'Last Updated', checked: false }
                ];
        }
        
        // Create checkboxes
        properties.forEach(prop => {
            const div = document.createElement('div');
            div.className = 'property-checkbox';
            
            const input = document.createElement('input');
            input.type = 'checkbox';
            input.id = `property-${prop.id}`;
            input.value = prop.id;
            input.checked = prop.checked;
            
            const label = document.createElement('label');
            label.htmlFor = `property-${prop.id}`;
            label.textContent = prop.label;
            
            div.appendChild(input);
            div.appendChild(label);
            container.appendChild(div);
        });
    },
    
    /**
     * Export data in the specified format
     * @param {string} dataType - Type of data to export (molecules, mixtures, etc.)
     * @param {string} itemId - ID of the item to export (optional)
     * @param {string} format - Export format (csv, json, excel, pdf)
     * @param {boolean} includeRelated - Whether to include related data
     * @param {Object} filters - Filters to apply to the data (optional)
     * @param {Object} metadata - Metadata to include with the export (optional)
     * @param {Object} formatOptions - Format-specific options (optional)
     */
    exportData: function(dataType, itemId, format, includeRelated, filters = {}, metadata = {}, formatOptions = {}) {
        // If offline, add to queue and return
        if (!navigator.onLine) {
            this.addToOfflineQueue({
                type: 'export',
                dataType: dataType,
                itemId: itemId,
                format: format,
                includeRelated: includeRelated,
                filters: filters,
                metadata: metadata,
                formatOptions: formatOptions
            });
            return;
        }
        
        // Show loading indicator with progress
        const loadingIndicator = document.getElementById('loading-indicator') || this.createLoadingIndicator();
        loadingIndicator.style.display = 'block';
        
        // Add progress bar if not present
        if (!loadingIndicator.querySelector('.progress-bar')) {
            const progressContainer = document.createElement('div');
            progressContainer.className = 'progress-container';
            
            const progressBar = document.createElement('div');
            progressBar.className = 'progress-bar';
            progressBar.style.width = '0%';
            
            const progressText = document.createElement('div');
            progressText.className = 'progress-text';
            progressText.textContent = 'Preparing export...';
            
            progressContainer.appendChild(progressBar);
            progressContainer.appendChild(progressText);
            loadingIndicator.appendChild(progressContainer);
        }
        
        // Update progress
        const updateProgress = (percent, message) => {
            const progressBar = loadingIndicator.querySelector('.progress-bar');
            const progressText = loadingIndicator.querySelector('.progress-text');
            
            if (progressBar) progressBar.style.width = `${percent}%`;
            if (progressText) progressText.textContent = message || `Exporting... ${percent}%`;
        };

        // Prepare request data
        const requestData = {
            format: format,
            data_type: dataType,
            include_related: includeRelated
        };

        // Add item ID if provided
        if (itemId) {
            requestData.id = itemId;
        }
        
        // Add filters if provided
        if (Object.keys(filters).length > 0) {
            requestData.filters = filters;
        }
        
        // Add metadata if provided
        if (Object.keys(metadata).length > 0) {
            requestData.metadata = metadata;
        }
        
        // Add format-specific options if provided
        if (Object.keys(formatOptions).length > 0) {
            requestData.format_options = formatOptions;
        }
        
        // Update progress
        updateProgress(10, 'Preparing export request...');

        // Make API request
        fetch('/api/v1/export', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
                'Authorization': `Bearer ${this.getAuthToken()}`
            },
            body: JSON.stringify(requestData)
        })
        .then(response => {
            updateProgress(50, 'Processing export data...');
            
            if (!response.ok) {
                return response.json().then(data => {
                    throw new Error(data.message || 'Export failed');
                });
            }
            
            updateProgress(75, 'Preparing download...');
            
            // For file downloads, get the blob and create a download link
            return response.blob().then(blob => {
                updateProgress(90, 'Creating download link...');
                
                const url = window.URL.createObjectURL(blob);
                const a = document.createElement('a');
                a.style.display = 'none';
                a.href = url;
                
                // Generate filename
                let filename = '';
                if (metadata && metadata.title) {
                    // Use metadata title if provided
                    filename = `${metadata.title.replace(/[^a-z0-9]/gi, '_').toLowerCase()}`;
                } else {
                    filename = `${dataType}`;
                }
                
                // Add item ID if available and not in title
                if (itemId && !filename.includes(itemId)) {
                    filename += `_${itemId.substring(0, 8)}`;
                }
                
                // Add timestamp if requested
                if (!metadata || metadata.includeTimestamp !== false) {
                    const timestamp = new Date().toISOString().replace(/[:.]/g, '-');
                    filename += `_${timestamp}`;
                }
                
                // Add format extension
                filename += `.${format}`;
                
                a.download = filename;
                
                document.body.appendChild(a);
                a.click();
                window.URL.revokeObjectURL(url);
                a.remove();
                
                updateProgress(100, 'Export complete!');
                
                // Log export for analytics
                this.logExportActivity(dataType, itemId, format);
                
                // Show success notification
                this.showNotification('success', `Export complete: ${filename}`);
            });
        })
        .catch(error => {
            console.error('Export error:', error);
            this.showNotification('error', `Export failed: ${error.message}`);
        })
        .finally(() => {
            // Hide loading indicator after a short delay to show completion
            setTimeout(() => {
                if (loadingIndicator) loadingIndicator.style.display = 'none';
                // Remove progress elements
                const progressContainer = loadingIndicator.querySelector('.progress-container');
                if (progressContainer) progressContainer.remove();
            }, 1000);
        });
    },
    
    /**
     * Create loading indicator if it doesn't exist
     * @returns {HTMLElement} Loading indicator element
     */
    createLoadingIndicator: function() {
        let loadingIndicator = document.getElementById('loading-indicator');
        
        if (!loadingIndicator) {
            loadingIndicator = document.createElement('div');
            loadingIndicator.id = 'loading-indicator';
            
            const spinner = document.createElement('div');
            spinner.className = 'spinner';
            
            loadingIndicator.appendChild(spinner);
            document.body.appendChild(loadingIndicator);
        }
        
        return loadingIndicator;
    },
    
    /**
     * Log export activity for analytics
     * @param {string} dataType - Type of data exported
     * @param {string} itemId - ID of the exported item
     * @param {string} format - Export format
     */
    logExportActivity: function(dataType, itemId, format) {
        // Send analytics data to server
        fetch('/api/v1/analytics/log', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
                'Authorization': `Bearer ${this.getAuthToken()}`
            },
            body: JSON.stringify({
                activity_type: 'export',
                data_type: dataType,
                item_id: itemId,
                format: format,
                timestamp: new Date().toISOString()
            })
        }).catch(error => {
            console.warn('Failed to log export activity:', error);
        });
    },
    
    /**
     * Batch export multiple items
     * @param {string} dataType - Type of data to export
     * @param {Array} itemIds - Array of item IDs to export
     * @param {string} format - Export format
     * @param {Object} options - Export options
     */
    batchExport: function(dataType, itemIds, format, options = {}) {
        // If offline, add to queue and return
        if (!navigator.onLine) {
            this.addToOfflineQueue({
                type: 'batch',
                dataType: dataType,
                items: itemIds,
                format: format,
                options: options
            });
            return;
        }
        
        // Show loading indicator with progress
        const loadingIndicator = this.createLoadingIndicator();
        loadingIndicator.style.display = 'block';
        
        // Add progress bar if not present
        if (!loadingIndicator.querySelector('.progress-bar')) {
            const progressContainer = document.createElement('div');
            progressContainer.className = 'progress-container';
            
            const progressBar = document.createElement('div');
            progressBar.className = 'progress-bar';
            progressBar.style.width = '0%';
            
            const progressText = document.createElement('div');
            progressText.className = 'progress-text';
            progressText.textContent = 'Preparing batch export...';
            
            progressContainer.appendChild(progressBar);
            progressContainer.appendChild(progressText);
            loadingIndicator.appendChild(progressContainer);
        }
        
        // Update progress
        const updateProgress = (percent, message) => {
            const progressBar = loadingIndicator.querySelector('.progress-bar');
            const progressText = loadingIndicator.querySelector('.progress-text');
            
            if (progressBar) progressBar.style.width = `${percent}%`;
            if (progressText) progressText.textContent = message || `Exporting... ${percent}%`;
        };
        
        updateProgress(10, `Preparing batch export of ${itemIds.length} items...`);
        
        // Prepare request data
        const requestData = {
            data_type: dataType,
            item_ids: itemIds,
            format: format,
            include_related: options.includeRelated || false,
            combine_files: options.combineFiles || false,
            zip_files: options.zipFiles || false
        };
        
        // Add any additional options
        if (options.filters) requestData.filters = options.filters;
        if (options.metadata) requestData.metadata = options.metadata;
        if (options.formatOptions) requestData.format_options = options.formatOptions;
        
        // Make API request
        fetch('/api/v1/export/batch', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
                'Authorization': `Bearer ${this.getAuthToken()}`
            },
            body: JSON.stringify(requestData)
        })
        .then(response => {
            updateProgress(50, 'Processing batch export...');
            
            if (!response.ok) {
                return response.json().then(data => {
                    throw new Error(data.message || 'Batch export failed');
                });
            }
            
            updateProgress(75, 'Preparing download...');
            
            // For file downloads, get the blob and create a download link
            return response.blob().then(blob => {
                updateProgress(90, 'Creating download link...');
                
                const url = window.URL.createObjectURL(blob);
                const a = document.createElement('a');
                a.style.display = 'none';
                a.href = url;
                
                // Generate filename
                const timestamp = new Date().toISOString().replace(/[:.]/g, '-');
                let filename = `${dataType}_batch_${itemIds.length}_items_${timestamp}`;
                
                // Add appropriate extension
                if (options.zipFiles) {
                    filename += '.zip';
                } else if (options.combineFiles) {
                    filename += `.${format}`;
                } else {
                    filename += `.${format}`;
                }
                
                a.download = filename;
                
                document.body.appendChild(a);
                a.click();
                window.URL.revokeObjectURL(url);
                a.remove();
                
                updateProgress(100, 'Batch export complete!');
                
                // Show success notification
                this.showNotification('success', `Batch export complete: ${itemIds.length} items exported`);
            });
        })
        .catch(error => {
            console.error('Batch export error:', error);
            this.showNotification('error', `Batch export failed: ${error.message}`);
        })
        .finally(() => {
            // Hide loading indicator after a short delay to show completion
            setTimeout(() => {
                if (loadingIndicator) loadingIndicator.style.display = 'none';
                // Remove progress elements
                const progressContainer = loadingIndicator.querySelector('.progress-container');
                if (progressContainer) progressContainer.remove();
            }, 1000);
        });
    },

    /**
     * Preview visualization before export
     * @param {Event} event - Click event
     */
    previewVisualization: function(event) {
        event.preventDefault();
        const btn = event.currentTarget;
        const container = btn.closest('.visualization-export-container');
        const chartType = container.dataset.chartType;
        const chartId = container.dataset.chartId;
        
        // Get visualization options
        const width = parseInt(container.querySelector('#visualization-width').value, 10) || 800;
        const height = parseInt(container.querySelector('#visualization-height').value, 10) || 600;
        const style = container.querySelector('#visualization-style').value;
        const title = container.querySelector('#visualization-title')?.value || '';
        const colorScheme = container.querySelector('#visualization-color-scheme')?.value || 'default';
        const fontFamily = container.querySelector('#visualization-font')?.value || 'Arial';
        const showLegend = container.querySelector('#visualization-show-legend')?.checked || false;
        const showGrid = container.querySelector('#visualization-show-grid')?.checked || false;
        const showDataLabels = container.querySelector('#visualization-show-data-labels')?.checked || false;
        
        // Show loading indicator in preview area
        const previewArea = container.querySelector('.visualization-preview-area') ||
            this.createPreviewArea(container);
        
        previewArea.innerHTML = '<div class="preview-loading">Loading preview...</div>';
        previewArea.style.display = 'block';
        
        // Get chart data from the chart element
        let chartData;
        if (chartId) {
            const chartElement = document.getElementById(chartId);
            if (chartElement && chartElement.__chartData) {
                chartData = chartElement.__chartData;
            } else {
                previewArea.innerHTML = '<div class="preview-error">Chart data not available</div>';
                return;
            }
        } else {
            previewArea.innerHTML = '<div class="preview-error">Chart ID not specified</div>';
            return;
        }
        
        // Make API request for preview
        fetch('/api/v1/export/visualization/preview', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
                'Authorization': `Bearer ${this.getAuthToken()}`
            },
            body: JSON.stringify({
                chart_type: chartType,
                data: chartData,
                width: width,
                height: height,
                title: title,
                style: style,
                color_scheme: colorScheme,
                font_family: fontFamily,
                show_legend: showLegend,
                show_grid: showGrid,
                show_data_labels: showDataLabels
            })
        })
        .then(response => {
            if (!response.ok) {
                return response.json().then(data => {
                    throw new Error(data.message || 'Preview generation failed');
                });
            }
            return response.json();
        })
        .then(data => {
            if (data.preview_html) {
                previewArea.innerHTML = data.preview_html;
            } else {
                previewArea.innerHTML = '<div class="preview-error">Preview not available</div>';
            }
        })
        .catch(error => {
            console.error('Visualization preview error:', error);
            previewArea.innerHTML = `<div class="preview-error">Error: ${error.message}</div>`;
        });
    },
    
    /**
     * Create preview area for visualization
     * @param {HTMLElement} container - Container element
     * @returns {HTMLElement} Preview area element
     */
    createPreviewArea: function(container) {
        const previewArea = document.createElement('div');
        previewArea.className = 'visualization-preview-area';
        previewArea.style.marginTop = '15px';
        previewArea.style.border = '1px solid #ccc';
        previewArea.style.padding = '10px';
        previewArea.style.backgroundColor = '#f9f9f9';
        previewArea.style.borderRadius = '4px';
        previewArea.style.minHeight = '200px';
        previewArea.style.display = 'none';
        
        container.appendChild(previewArea);
        return previewArea;
    },
    
    /**
     * Save visualization template
     * @param {Event} event - Click event
     */
    saveVisualizationTemplate: function(event) {
        event.preventDefault();
        const btn = event.currentTarget;
        const container = btn.closest('.visualization-export-container');
        
        // Create dialog for template name and description
        const dialog = document.createElement('div');
        dialog.className = 'template-dialog';
        
        dialog.innerHTML = `
            <div class="template-dialog-content">
                <h3>Save Visualization Template</h3>
                <div class="template-form">
                    <div class="template-form-group">
                        <label for="template-name">Template Name:</label>
                        <input type="text" id="template-name" placeholder="Enter template name">
                    </div>
                    <div class="template-form-group">
                        <label for="template-description">Description:</label>
                        <textarea id="template-description" placeholder="Enter template description"></textarea>
                    </div>
                </div>
                <div class="template-dialog-actions">
                    <button class="cancel-template-btn">Cancel</button>
                    <button class="save-template-btn">Save Template</button>
                </div>
            </div>
        `;
        
        // Add event listeners
        dialog.querySelector('.cancel-template-btn').addEventListener('click', () => {
            dialog.remove();
        });
        
        dialog.querySelector('.save-template-btn').addEventListener('click', () => {
            const name = dialog.querySelector('#template-name').value.trim();
            const description = dialog.querySelector('#template-description').value.trim();
            
            if (!name) {
                this.showNotification('error', 'Please enter a template name');
                return;
            }
            
            // Get template settings
            const templateSettings = {
                name: name,
                description: description,
                chart_type: container.dataset.chartType,
                width: parseInt(container.querySelector('#visualization-width').value, 10) || 800,
                height: parseInt(container.querySelector('#visualization-height').value, 10) || 600,
                style: container.querySelector('#visualization-style').value,
                title_format: container.querySelector('#visualization-title')?.value || '',
                color_scheme: container.querySelector('#visualization-color-scheme')?.value || 'default',
                font_family: container.querySelector('#visualization-font')?.value || 'Arial',
                show_legend: container.querySelector('#visualization-show-legend')?.checked || false,
                show_grid: container.querySelector('#visualization-show-grid')?.checked || false,
                show_data_labels: container.querySelector('#visualization-show-data-labels')?.checked || false
            };
            
            // Save template to server
            fetch('/api/v1/templates/visualization', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                    'Authorization': `Bearer ${this.getAuthToken()}`
                },
                body: JSON.stringify(templateSettings)
            })
            .then(response => {
                if (!response.ok) {
                    return response.json().then(data => {
                        throw new Error(data.message || 'Failed to save template');
                    });
                }
                return response.json();
            })
            .then(data => {
                // Add to local templates
                if (!this.visualizationTemplates) {
                    this.visualizationTemplates = [];
                }
                
                templateSettings.id = data.template_id;
                this.visualizationTemplates.push(templateSettings);
                
                // Update localStorage
                localStorage.setItem('visualization-templates', JSON.stringify(this.visualizationTemplates));
                
                // Update template selects
                this.populateTemplateSelects();
                
                // Show success notification
                this.showNotification('success', 'Template saved successfully');
                
                // Close dialog
                dialog.remove();
            })
            .catch(error => {
                console.error('Template save error:', error);
                this.showNotification('error', `Failed to save template: ${error.message}`);
            });
        });
        
        // Add to document
        document.body.appendChild(dialog);
    },
    
    /**
     * Load visualization template
     * @param {Event} event - Change event
     */
    loadVisualizationTemplate: function(event) {
        const select = event.target;
        const templateId = select.value;
        
        if (templateId === 'none') return;
        
        const container = select.closest('.visualization-export-container');
        
        // Find template
        const template = this.visualizationTemplates.find(t => t.id === templateId);
        if (!template) {
            this.showNotification('error', 'Template not found');
            return;
        }
        
        // Apply template settings
        if (container.querySelector('#visualization-width')) {
            container.querySelector('#visualization-width').value = template.width;
        }
        
        if (container.querySelector('#visualization-height')) {
            container.querySelector('#visualization-height').value = template.height;
        }
        
        if (container.querySelector('#visualization-style')) {
            container.querySelector('#visualization-style').value = template.style;
        }
        
        if (container.querySelector('#visualization-title')) {
            container.querySelector('#visualization-title').value = template.title_format;
        }
        
        if (container.querySelector('#visualization-color-scheme')) {
            container.querySelector('#visualization-color-scheme').value = template.color_scheme;
        }
        
        if (container.querySelector('#visualization-font')) {
            container.querySelector('#visualization-font').value = template.font_family;
        }
        
        if (container.querySelector('#visualization-show-legend')) {
            container.querySelector('#visualization-show-legend').checked = template.show_legend;
        }
        
        if (container.querySelector('#visualization-show-grid')) {
            container.querySelector('#visualization-show-grid').checked = template.show_grid;
        }
        
        if (container.querySelector('#visualization-show-data-labels')) {
            container.querySelector('#visualization-show-data-labels').checked = template.show_data_labels;
        }
        
        // Update preview if available
        const previewBtn = container.querySelector('.preview-visualization-btn');
        if (previewBtn) {
            previewBtn.click();
        }
        
        this.showNotification('info', `Template "${template.name}" applied`);
    },
    
    /**
     * Export visualization
     * @param {string} chartType - Type of chart (property_comparison, mixture_composition, etc.)
     * @param {Object} chartData - Data for the visualization
     * @param {string} format - Export format (png, svg, pdf)
     * @param {Object} options - Additional options (width, height, title, style)
     */
    exportVisualization: function(chartType, chartData, format, options = {}) {
        // If offline, add to queue and return
        if (!navigator.onLine) {
            this.addToOfflineQueue({
                type: 'visualization',
                chartType: chartType,
                chartData: chartData,
                format: format,
                options: options
            });
            return;
        }
        
        // Show loading indicator with progress
        const loadingIndicator = this.createLoadingIndicator();
        loadingIndicator.style.display = 'block';
        
        // Add progress bar if not present
        if (!loadingIndicator.querySelector('.progress-bar')) {
            const progressContainer = document.createElement('div');
            progressContainer.className = 'progress-container';
            
            const progressBar = document.createElement('div');
            progressBar.className = 'progress-bar';
            progressBar.style.width = '0%';
            
            const progressText = document.createElement('div');
            progressText.className = 'progress-text';
            progressText.textContent = 'Preparing visualization export...';
            
            progressContainer.appendChild(progressBar);
            progressContainer.appendChild(progressText);
            loadingIndicator.appendChild(progressContainer);
        }
        
        // Update progress
        const updateProgress = (percent, message) => {
            const progressBar = loadingIndicator.querySelector('.progress-bar');
            const progressText = loadingIndicator.querySelector('.progress-text');
            
            if (progressBar) progressBar.style.width = `${percent}%`;
            if (progressText) progressText.textContent = message || `Exporting... ${percent}%`;
        };
        
        updateProgress(10, 'Preparing visualization...');

        // Prepare request data
        const requestData = {
            chart_type: chartType,
            data: chartData,
            format: format,
            width: options.width || 800,
            height: options.height || 600,
            title: options.title || '',
            style: options.style || 'default'
        };
        
        // Add additional visualization options if provided
        if (options.colorScheme) requestData.color_scheme = options.colorScheme;
        if (options.fontFamily) requestData.font_family = options.fontFamily;
        if (options.showLegend !== undefined) requestData.show_legend = options.showLegend;
        if (options.showGrid !== undefined) requestData.show_grid = options.showGrid;
        if (options.showDataLabels !== undefined) requestData.show_data_labels = options.showDataLabels;
        if (options.templateId) requestData.template_id = options.templateId;
        
        // Add advanced chart-specific options
        if (options.advanced) {
            requestData.advanced_options = options.advanced;
        }

        // Make API request
        fetch('/api/v1/export/visualization', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
                'Authorization': `Bearer ${this.getAuthToken()}`
            },
            body: JSON.stringify(requestData)
        })
        .then(response => {
            updateProgress(50, 'Processing visualization...');
            
            if (!response.ok) {
                return response.json().then(data => {
                    throw new Error(data.message || 'Visualization export failed');
                });
            }
            
            updateProgress(75, 'Preparing download...');
            
            // For file downloads, get the blob and create a download link
            return response.blob().then(blob => {
                updateProgress(90, 'Creating download link...');
                
                const url = window.URL.createObjectURL(blob);
                const a = document.createElement('a');
                a.style.display = 'none';
                a.href = url;
                
                // Generate filename
                let filename = '';
                if (options.title) {
                    // Use title if provided
                    filename = `${options.title.replace(/[^a-z0-9]/gi, '_').toLowerCase()}`;
                } else {
                    filename = `visualization_${chartType}`;
                }
                
                // Add timestamp
                const timestamp = new Date().toISOString().replace(/[:.]/g, '-');
                filename += `_${timestamp}.${format}`;
                
                a.download = filename;
                
                document.body.appendChild(a);
                a.click();
                window.URL.revokeObjectURL(url);
                a.remove();
                
                updateProgress(100, 'Visualization export complete!');
                
                // Log export for analytics
                this.logExportActivity('visualization', null, format);
                
                // Show success notification
                this.showNotification('success', `Visualization exported as ${format.toUpperCase()}`);
            });
        })
        .catch(error => {
            console.error('Visualization export error:', error);
            this.showNotification('error', `Visualization export failed: ${error.message}`);
        })
        .finally(() => {
            // Hide loading indicator after a short delay to show completion
            setTimeout(() => {
                if (loadingIndicator) loadingIndicator.style.display = 'none';
                // Remove progress elements
                const progressContainer = loadingIndicator.querySelector('.progress-container');
                if (progressContainer) progressContainer.remove();
            }, 1000);
        });
    },

    /**
     * Generate a report
     * @param {string} title - Report title
     * @param {Array} sections - Report sections
     * @param {Object} options - Additional options
     */
    generateReport: function(title, sections, options = {}) {
        // Show loading indicator
        const loadingIndicator = document.getElementById('loading-indicator');
        if (loadingIndicator) loadingIndicator.style.display = 'block';

        // Prepare request data
        const requestData = {
            title: title,
            sections: sections,
            include_visualizations: options.includeVisualizations !== false,
            include_data_tables: options.includeDataTables !== false
        };

        // Add template ID if provided
        if (options.templateId) {
            requestData.template_id = options.templateId;
        }

        // Make API request
        fetch('/api/v1/export/report', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
                'Authorization': `Bearer ${this.getAuthToken()}`
            },
            body: JSON.stringify(requestData)
        })
        .then(response => {
            if (!response.ok) {
                return response.json().then(data => {
                    throw new Error(data.message || 'Report generation failed');
                });
            }
            
            // For file downloads, get the blob and create a download link
            return response.blob().then(blob => {
                const url = window.URL.createObjectURL(blob);
                const a = document.createElement('a');
                a.style.display = 'none';
                a.href = url;
                
                // Generate filename
                const timestamp = new Date().toISOString().replace(/[:.]/g, '-');
                const filename = `report_${timestamp}.pdf`;
                a.download = filename;
                
                document.body.appendChild(a);
                a.click();
                window.URL.revokeObjectURL(url);
                a.remove();
            });
        })
        .catch(error => {
            console.error('Report generation error:', error);
            this.showNotification('error', `Report generation failed: ${error.message}`);
        })
        .finally(() => {
            // Hide loading indicator
            if (loadingIndicator) loadingIndicator.style.display = 'none';
        });
    },

    /**
     * Get authentication token
     * @returns {string} Authentication token
     */
    getAuthToken: function() {
        // Get token from Supabase auth
        const supabase = window.supabase;
        if (supabase && supabase.auth) {
            const session = supabase.auth.session();
            return session ? session.access_token : '';
        }
        return '';
    },

    /**
     * Show notification
     * @param {string} type - Notification type (success, error, info)
     * @param {string} message - Notification message
     */
    showNotification: function(type, message) {
        // Create notification element
        const notification = document.createElement('div');
        notification.className = `notification ${type}`;
        notification.textContent = message;
        
        // Add close button
        const closeBtn = document.createElement('button');
        closeBtn.className = 'notification-close';
        closeBtn.innerHTML = '&times;';
        closeBtn.addEventListener('click', () => {
            notification.remove();
        });
        notification.appendChild(closeBtn);
        
        // Add to document
        document.body.appendChild(notification);
        
        // Auto-remove after 5 seconds
        setTimeout(() => {
            notification.remove();
        }, 5000);
    }
};

// Sharing functionality
const ShareManager = {
    /**
     * Initialize sharing functionality
     */
    init: function() {
        this.bindEvents();
        this.setupAccessibilityFeatures();
        this.loadSharedItems();
    },

    /**
     * Bind event listeners
     */
    bindEvents: function() {
        // Share buttons
        document.querySelectorAll('.share-btn').forEach(btn => {
            btn.addEventListener('click', this.handleShareClick.bind(this));
        });
        
        // Open comprehensive share dialog button
        document.querySelectorAll('.open-share-dialog-btn').forEach(btn => {
            btn.addEventListener('click', this.openShareDialog.bind(this));
        });
        
        // Manage shared items button
        document.querySelectorAll('.manage-shared-items-btn').forEach(btn => {
            btn.addEventListener('click', this.openSharedItemsManager.bind(this));
        });
        
        // View analytics button
        document.querySelectorAll('.view-sharing-analytics-btn').forEach(btn => {
            btn.addEventListener('click', this.openSharingAnalytics.bind(this));
        });

        // Share type selection
        document.querySelectorAll('.share-type-select').forEach(select => {
            select.addEventListener('change', this.updateShareOptions.bind(this));
        });

        // Password protection checkbox
        document.querySelectorAll('.password-protected-checkbox').forEach(checkbox => {
            checkbox.addEventListener('change', this.togglePasswordField.bind(this));
        });
    },
    
    /**
     * Setup accessibility features
     */
    setupAccessibilityFeatures: function() {
        // Add ARIA attributes to share controls
        document.querySelectorAll('.share-btn').forEach(btn => {
            btn.setAttribute('aria-label', 'Share item');
            btn.setAttribute('role', 'button');
        });
        
        // Add keyboard navigation
        document.querySelectorAll('.share-type-select').forEach(select => {
            select.addEventListener('keydown', (e) => {
                if (e.key === 'Enter') {
                    e.target.click();
                }
            });
        });
    },
    
    /**
     * Load shared items from server
     */
    loadSharedItems: function() {
        fetch('/api/v1/share/my-shares', {
            method: 'GET',
            headers: {
                'Authorization': `Bearer ${this.getAuthToken()}`
            }
        })
        .then(response => {
            if (!response.ok) {
                throw new Error('Failed to load shared items');
            }
            return response.json();
        })
        .then(data => {
            this.sharedItems = data.items;
            
            // Update shared items count if element exists
            const sharedItemsCount = document.getElementById('shared-items-count');
            if (sharedItemsCount) {
                sharedItemsCount.textContent = this.sharedItems.length;
            }
        })
        .catch(error => {
            console.warn('Could not load shared items:', error);
        });
    },

    /**
     * Handle share button click
     * @param {Event} event - Click event
     */
    handleShareClick: function(event) {
        event.preventDefault();
        const btn = event.currentTarget;
        const dataType = btn.dataset.type;
        const itemId = btn.dataset.id;
        const shareContainer = btn.closest('.share-container');
        const shareType = shareContainer.querySelector('.share-type-select').value;
        
        // Get share options
        const options = {
            passwordProtected: shareContainer.querySelector('.password-protected-checkbox')?.checked || false,
            password: shareContainer.querySelector('.share-password-input')?.value || '',
            expiration: parseInt(shareContainer.querySelector('.share-expiration-select')?.value || '86400', 10)
        };
        
        // Get recipients for email sharing
        if (shareType === 'email') {
            options.recipients = Array.from(
                shareContainer.querySelector('.share-recipients-input').value.split(',')
            ).map(email => email.trim()).filter(email => email);
            
            options.message = shareContainer.querySelector('.share-message-textarea')?.value || '';
        }
        
        this.shareItem(dataType, itemId, shareType, options);
    },

    /**
     * Update share options based on selected share type
     * @param {Event} event - Change event
     */
    updateShareOptions: function(event) {
        const shareType = event.target.value;
        const container = event.target.closest('.share-container');
        
        // Show/hide share type-specific options
        container.querySelectorAll('.share-option').forEach(option => {
            option.style.display = 'none';
        });
        
        container.querySelectorAll(`.share-option.${shareType}`).forEach(option => {
            option.style.display = 'block';
        });
    },

    /**
     * Toggle password field visibility
     * @param {Event} event - Change event
     */
    togglePasswordField: function(event) {
        const isChecked = event.target.checked;
        const container = event.target.closest('.share-container');
        const passwordField = container.querySelector('.share-password-field');
        
        if (passwordField) {
            passwordField.style.display = isChecked ? 'block' : 'none';
        }
    },

    /**
     * Share an item
     * @param {string} dataType - Type of data to share (molecules, mixtures, etc.)
     * @param {string} itemId - ID of the item to share
     * @param {string} shareType - Share type (link, email, embed)
     * @param {Object} options - Share options
     */
    shareItem: function(dataType, itemId, shareType, options = {}) {
        // Show loading indicator
        const loadingIndicator = document.getElementById('loading-indicator');
        if (loadingIndicator) loadingIndicator.style.display = 'block';

        // Prepare request data
        const requestData = {
            share_type: shareType,
            data_type: dataType,
            id: itemId,
            password_protected: options.passwordProtected || false,
            expiration: options.expiration || 86400 // Default to 24 hours
        };

        // Add password if protected
        if (options.passwordProtected && options.password) {
            requestData.password = options.password;
        }

        // Add recipients and message for email sharing
        if (shareType === 'email') {
            if (!options.recipients || options.recipients.length === 0) {
                this.showNotification('error', 'Please enter at least one recipient email address');
                if (loadingIndicator) loadingIndicator.style.display = 'none';
                return;
            }
            
            requestData.recipients = options.recipients;
            
            if (options.message) {
                requestData.message = options.message;
            }
        }

        // Make API request
        fetch('/api/v1/share', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
                'Authorization': `Bearer ${this.getAuthToken()}`
            },
            body: JSON.stringify(requestData)
        })
        .then(response => {
            if (!response.ok) {
                return response.json().then(data => {
                    throw new Error(data.message || 'Sharing failed');
                });
            }
            
            return response.json();
        })
        .then(data => {
            // Handle successful share
            if (shareType === 'link') {
                this.showShareLinkDialog(data.share_url, data.password_protected);
            } else if (shareType === 'email') {
                this.showNotification('success', `Item shared via email to ${options.recipients.join(', ')}`);
            } else if (shareType === 'embed') {
                this.showEmbedCodeDialog(data.embed_code);
            }
        })
        .catch(error => {
            console.error('Sharing error:', error);
            this.showNotification('error', `Sharing failed: ${error.message}`);
        })
        .finally(() => {
            // Hide loading indicator
            if (loadingIndicator) loadingIndicator.style.display = 'none';
        });
    },

    /**
     * Show share link dialog
     * @param {string} shareUrl - Share URL
     * @param {boolean} isPasswordProtected - Whether the share is password protected
     */
    showShareLinkDialog: function(shareUrl, isPasswordProtected) {
        // Create dialog element
        const dialog = document.createElement('div');
        dialog.className = 'share-dialog';
        
        // Create dialog content
        dialog.innerHTML = `
            <div class="share-dialog-content">
                <h3>Share Link</h3>
                <p>Use the link below to share this item:</p>
                <div class="share-url-container">
                    <input type="text" class="share-url-input" value="${shareUrl}" readonly>
                    <button class="copy-url-btn">Copy</button>
                </div>
                ${isPasswordProtected ? '<p class="password-notice">This link is password protected.</p>' : ''}
                <div class="share-dialog-actions">
                    <button class="share-dialog-close-btn">Close</button>
                </div>
            </div>
        `;
        
        // Add event listeners
        dialog.querySelector('.copy-url-btn').addEventListener('click', () => {
            const input = dialog.querySelector('.share-url-input');
            input.select();
            document.execCommand('copy');
            this.showNotification('success', 'Link copied to clipboard');
        });
        
        dialog.querySelector('.share-dialog-close-btn').addEventListener('click', () => {
            dialog.remove();
        });
        
        // Add to document
        document.body.appendChild(dialog);
    },

    /**
     * Show embed code dialog
     * @param {string} embedCode - Embed code
     */
    showEmbedCodeDialog: function(embedCode) {
        // Create dialog element
        const dialog = document.createElement('div');
        dialog.className = 'share-dialog';
        
        // Create dialog content
        dialog.innerHTML = `
            <div class="share-dialog-content">
                <h3>Embed Code</h3>
                <p>Use the code below to embed this item in your website:</p>
                <div class="embed-code-container">
                    <textarea class="embed-code-textarea" readonly>${embedCode}</textarea>
                    <button class="copy-embed-btn">Copy</button>
                </div>
                <div class="share-dialog-actions">
                    <button class="share-dialog-close-btn">Close</button>
                </div>
            </div>
        `;
        
        // Add event listeners
        dialog.querySelector('.copy-embed-btn').addEventListener('click', () => {
            const textarea = dialog.querySelector('.embed-code-textarea');
            textarea.select();
            document.execCommand('copy');
            this.showNotification('success', 'Embed code copied to clipboard');
        });
        
        dialog.querySelector('.share-dialog-close-btn').addEventListener('click', () => {
            dialog.remove();
        });
        
        // Add to document
        document.body.appendChild(dialog);
    },

    /**
     * Get authentication token
     * @returns {string} Authentication token
     */
    getAuthToken: function() {
        // Get token from Supabase auth
        const supabase = window.supabase;
        if (supabase && supabase.auth) {
            const session = supabase.auth.session();
            return session ? session.access_token : '';
        }
        return '';
    },

    /**
     * Show notification
     * @param {string} type - Notification type (success, error, info)
     * @param {string} message - Notification message
     */
    showNotification: function(type, message) {
        // Create notification element
        const notification = document.createElement('div');
        notification.className = `notification ${type}`;
        notification.textContent = message;
        
        // Add close button
        const closeBtn = document.createElement('button');
        closeBtn.className = 'notification-close';
        closeBtn.innerHTML = '&times;';
        closeBtn.addEventListener('click', () => {
            notification.remove();
        });
        notification.appendChild(closeBtn);
        
        // Add to document
        document.body.appendChild(notification);
        
        // Auto-remove after 5 seconds
        setTimeout(() => {
            notification.remove();
        }, 5000);
    }
};

// Initialize when DOM is ready
document.addEventListener('DOMContentLoaded', function() {
    ExportManager.init();
    ShareManager.init();
});
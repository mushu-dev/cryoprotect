/**
 * CryoProtect Analyzer - Enhanced 3D Molecular Viewer
 *
 * This module provides functions for interactive 3D visualization of molecules using 3Dmol.js
 * with enhanced controls for rotation, zoom, and representation styles.
 *
 * Version: 2.0.0
 */

const MolecularViewer = (function() {
    // Cache for viewers and their states
    const viewers = {};
    const viewerStates = {};
    
    // Constants for control settings
    const ROTATION_SPEED = 0.5;
    const ZOOM_FACTOR = 0.1;
    const PAN_STEP = 1.0;
    
    /**
     * Initialize a 3D molecular viewer in a container
     * 
     * @param {string} containerId - ID of the container element
     * @param {Object} options - Viewer options
     * @returns {Object} - 3Dmol viewer instance
     */
    function initViewer(containerId, options = {}) {
        const container = document.getElementById(containerId);
        if (!container) {
            console.error(`Container element with ID ${containerId} not found`);
            return null;
        }
        
        // Default options
        const defaultOptions = {
            backgroundColor: 'white',
            antialias: true,
            width: container.clientWidth,
            height: 400,
            disableFog: true,
            defaultcolors: $3Dmol.rasmolElementColors
        };
        
        // Merge options
        const viewerOptions = { ...defaultOptions, ...options };
        
        // Create viewer
        const viewer = $3Dmol.createViewer(containerId, viewerOptions);
        
        // Store viewer in cache
        viewers[containerId] = viewer;
        
        // Initialize viewer state
        viewerStates[containerId] = {
            isSpinning: false,
            currentStyle: 'stick',
            showSurface: false,
            showLabels: false,
            showHydrogens: true,
            colorScheme: 'cyanCarbon',
            zoomLevel: 100, // percentage
            measurementMode: null,
            measurements: []
        };
        
        // Set up mouse wheel zoom
        setupMouseWheelZoom(containerId, container);
        
        // Set up keyboard navigation
        setupKeyboardNavigation(containerId);
        
        return viewer;
    }
    
    /**
     * Set up mouse wheel zoom for the viewer
     *
     * @param {string} containerId - ID of the container element
     * @param {HTMLElement} container - Container element
     */
    function setupMouseWheelZoom(containerId, container) {
        container.addEventListener('wheel', function(event) {
            event.preventDefault();
            
            const viewer = viewers[containerId];
            if (!viewer) return;
            
            // Determine zoom direction
            const delta = event.deltaY || event.detail || event.wheelDelta;
            
            if (delta > 0) {
                zoomOut(containerId);
            } else {
                zoomIn(containerId);
            }
        });
    }
    
    /**
     * Set up keyboard navigation for the viewer
     *
     * @param {string} containerId - ID of the container element
     */
    function setupKeyboardNavigation(containerId) {
        const container = document.getElementById(containerId);
        if (!container) return;
        
        // Make container focusable
        if (!container.hasAttribute('tabindex')) {
            container.setAttribute('tabindex', '0');
        }
        
        // Add ARIA label for accessibility
        if (!container.hasAttribute('aria-label')) {
            container.setAttribute('aria-label', 'Molecular 3D viewer. Use arrow keys to rotate, +/- to zoom, WASD to pan');
        }
        
        container.addEventListener('keydown', function(event) {
            const viewer = viewers[containerId];
            if (!viewer) return;
            
            switch (event.key) {
                case 'ArrowLeft':
                    rotateY(containerId, -ROTATION_SPEED);
                    event.preventDefault();
                    break;
                case 'ArrowRight':
                    rotateY(containerId, ROTATION_SPEED);
                    event.preventDefault();
                    break;
                case 'ArrowUp':
                    rotateX(containerId, -ROTATION_SPEED);
                    event.preventDefault();
                    break;
                case 'ArrowDown':
                    rotateX(containerId, ROTATION_SPEED);
                    event.preventDefault();
                    break;
                case '+':
                case '=':
                    zoomIn(containerId);
                    event.preventDefault();
                    break;
                case '-':
                case '_':
                    zoomOut(containerId);
                    event.preventDefault();
                    break;
                case 'w':
                case 'W':
                    pan(containerId, 0, PAN_STEP);
                    event.preventDefault();
                    break;
                case 's':
                case 'S':
                    pan(containerId, 0, -PAN_STEP);
                    event.preventDefault();
                    break;
                case 'a':
                case 'A':
                    pan(containerId, PAN_STEP, 0);
                    event.preventDefault();
                    break;
                case 'd':
                case 'D':
                    pan(containerId, -PAN_STEP, 0);
                    event.preventDefault();
                    break;
                case 'r':
                case 'R':
                    resetView(containerId);
                    event.preventDefault();
                    break;
            }
        });
    }
    
    /**
     * Load a molecule from SMILES string with enhanced visualization options
     *
     * @param {string} containerId - ID of the container element
     * @param {string} smiles - SMILES string
     * @param {Object} options - Visualization options
     */
    function loadMoleculeFromSmiles(containerId, smiles, options = {}) {
        // Get or create viewer
        let viewer = viewers[containerId];
        if (!viewer) {
            viewer = initViewer(containerId);
            if (!viewer) return;
        }
        
        // Clear viewer and measurements
        viewer.clear();
        viewerStates[containerId].measurements = [];
        
        // Default options
        const defaultOptions = {
            style: 'stick',
            colorScheme: 'cyanCarbon',
            showLabels: false,
            showHydrogens: true
        };
        
        // Merge options
        const visualOptions = { ...defaultOptions, ...options };
        
        // Update viewer state
        viewerStates[containerId].currentStyle = visualOptions.style;
        viewerStates[containerId].showLabels = visualOptions.showLabels;
        viewerStates[containerId].showHydrogens = visualOptions.showHydrogens;
        viewerStates[containerId].colorScheme = visualOptions.colorScheme;
        viewerStates[containerId].showSurface = visualOptions.showSurface || false;
        
        try {
            // Show loading indicator
            const container = document.getElementById(containerId);
            if (container) {
                container.innerHTML = `
                    <div class="text-center py-5">
                        <div class="spinner-border" role="status">
                            <span class="visually-hidden">Loading...</span>
                        </div>
                        <p class="mt-2">Generating 3D model...</p>
                    </div>
                `;
            }
            
            // Convert SMILES to 3D model
            const model = $3Dmol.createModelFrom(smiles, "smi", {}, function(model) {
                // Add model to viewer
                viewer.addModel(model);
                
                // Apply style based on the selected representation
                applyStyle(containerId, visualOptions.style);
                
                // Add labels if requested
                if (visualOptions.showLabels) {
                    viewer.addLabels(model.getAtoms(), {
                        font: 'Arial',
                        fontSize: 12,
                        showBackground: false,
                        alignment: 'center'
                    });
                }
                
                // Add surface if requested
                if (visualOptions.showSurface) {
                    viewer.addSurface($3Dmol.SurfaceType.VDW, {
                        opacity: 0.7,
                        color: 'white'
                    });
                }
                
                // Zoom to fit
                viewer.zoomTo();
                
                // Render
                viewer.render();
                
                // Announce to screen readers that the model is loaded
                const ariaLive = document.createElement('div');
                ariaLive.setAttribute('aria-live', 'polite');
                ariaLive.classList.add('sr-only');
                ariaLive.textContent = 'Molecule 3D model loaded successfully';
                document.body.appendChild(ariaLive);
                setTimeout(() => document.body.removeChild(ariaLive), 3000);
            });
        } catch (error) {
            console.error('Error loading molecule from SMILES:', error);
            
            // Show error message in container
            const container = document.getElementById(containerId);
            if (container) {
                container.innerHTML = `<div class="alert alert-danger">Error loading 3D model: ${error.message}</div>`;
            }
        }
    }
    
    /**
     * Load a molecule from PDB ID
     * 
     * @param {string} containerId - ID of the container element
     * @param {string} pdbId - PDB ID
     * @param {Object} options - Visualization options
     */
    function loadMoleculeFromPDB(containerId, pdbId, options = {}) {
        // Get or create viewer
        let viewer = viewers[containerId];
        if (!viewer) {
            viewer = initViewer(containerId);
            if (!viewer) return;
        }
        
        // Clear viewer
        viewer.clear();
        
        // Default options
        const defaultOptions = {
            style: 'cartoon',
            colorScheme: 'spectrum',
            showLabels: false,
            showSurface: false
        };
        
        // Merge options
        const visualOptions = { ...defaultOptions, ...options };
        
        try {
            // Load PDB
            $3Dmol.download(`pdb:${pdbId}`, viewer, {}, function() {
                // Apply style
                if (visualOptions.style === 'cartoon') {
                    viewer.setStyle({}, {cartoon: {color: visualOptions.colorScheme}});
                } else if (visualOptions.style === 'stick') {
                    viewer.setStyle({}, {stick: {radius: 0.15, colorscheme: visualOptions.colorScheme}});
                } else if (visualOptions.style === 'sphere') {
                    viewer.setStyle({}, {sphere: {scale: 0.25, colorscheme: visualOptions.colorScheme}});
                } else if (visualOptions.style === 'line') {
                    viewer.setStyle({}, {line: {colorscheme: visualOptions.colorScheme}});
                } else {
                    viewer.setStyle({}, {cartoon: {color: visualOptions.colorScheme}});
                }
                
                // Add surface if requested
                if (visualOptions.showSurface) {
                    viewer.addSurface($3Dmol.SurfaceType.VDW, {
                        opacity: 0.7,
                        color: 'white'
                    });
                }
                
                // Zoom to fit
                viewer.zoomTo();
                
                // Render
                viewer.render();
            });
        } catch (error) {
            console.error('Error loading molecule from PDB:', error);
            
            // Show error message in container
            const container = document.getElementById(containerId);
            if (container) {
                container.innerHTML = `<div class="alert alert-danger">Error loading 3D model: ${error.message}</div>`;
            }
        }
    }
    
    /**
     * Load a molecule from MOL data
     * 
     * @param {string} containerId - ID of the container element
     * @param {string} molData - MOL data
     * @param {Object} options - Visualization options
     */
    function loadMoleculeFromMOL(containerId, molData, options = {}) {
        // Get or create viewer
        let viewer = viewers[containerId];
        if (!viewer) {
            viewer = initViewer(containerId);
            if (!viewer) return;
        }
        
        // Clear viewer
        viewer.clear();
        
        // Default options
        const defaultOptions = {
            style: 'stick',
            colorScheme: 'cyanCarbon',
            showLabels: false,
            showHydrogens: true
        };
        
        // Merge options
        const visualOptions = { ...defaultOptions, ...options };
        
        try {
            // Add model to viewer
            const model = viewer.addModel(molData, "mol");
            
            // Apply style
            if (visualOptions.style === 'stick') {
                viewer.setStyle({}, {stick: {radius: 0.15, colorscheme: visualOptions.colorScheme}});
            } else if (visualOptions.style === 'sphere') {
                viewer.setStyle({}, {sphere: {scale: 0.25, colorscheme: visualOptions.colorScheme}});
            } else if (visualOptions.style === 'cartoon') {
                viewer.setStyle({}, {cartoon: {colorscheme: visualOptions.colorScheme}});
            } else if (visualOptions.style === 'line') {
                viewer.setStyle({}, {line: {colorscheme: visualOptions.colorScheme}});
            } else {
                viewer.setStyle({}, {stick: {radius: 0.15, colorscheme: visualOptions.colorScheme}});
            }
            
            // Add labels if requested
            if (visualOptions.showLabels) {
                viewer.addLabels(model.getAtoms(), {font: 'Arial', fontSize: 12, showBackground: false});
            }
            
            // Add surface if requested
            if (visualOptions.showSurface) {
                viewer.addSurface($3Dmol.SurfaceType.VDW, {
                    opacity: 0.7,
                    color: 'white'
                });
            }
            
            // Zoom to fit
            viewer.zoomTo();
            
            // Render
            viewer.render();
        } catch (error) {
            console.error('Error loading molecule from MOL:', error);
            
            // Show error message in container
            const container = document.getElementById(containerId);
            if (container) {
                container.innerHTML = `<div class="alert alert-danger">Error loading 3D model: ${error.message}</div>`;
            }
        }
    }
    
    /**
     * Apply a specific style to the molecule
     *
     * @param {string} containerId - ID of the container element
     * @param {string} style - Style to apply
     */
    function applyStyle(containerId, style) {
        const viewer = viewers[containerId];
        if (!viewer) return;
        
        const state = viewerStates[containerId];
        const colorScheme = state.colorScheme || 'cyanCarbon';
        
        // Update current style in state
        state.currentStyle = style;
        
        // Apply the selected style
        switch (style) {
            case 'stick':
                viewer.setStyle({}, {stick: {radius: 0.15, colorscheme: colorScheme}});
                break;
            case 'sphere':
                viewer.setStyle({}, {sphere: {scale: 0.25, colorscheme: colorScheme}});
                break;
            case 'line':
                viewer.setStyle({}, {line: {colorscheme: colorScheme}});
                break;
            case 'cartoon':
                viewer.setStyle({}, {cartoon: {colorscheme: colorScheme}});
                break;
            case 'ball-and-stick':
                viewer.setStyle({}, {
                    stick: {radius: 0.1, colorscheme: colorScheme},
                    sphere: {scale: 0.25, colorscheme: colorScheme}
                });
                break;
            case 'spacefill':
                viewer.setStyle({}, {sphere: {scale: 0.8, colorscheme: colorScheme}});
                break;
            default:
                viewer.setStyle({}, {stick: {radius: 0.15, colorscheme: colorScheme}});
        }
        
        viewer.render();
    }
    
    /**
     * Rotate the molecule around the X axis
     *
     * @param {string} containerId - ID of the container element
     * @param {number} angle - Rotation angle in degrees
     */
    function rotateX(containerId, angle) {
        const viewer = viewers[containerId];
        if (!viewer) return;
        
        viewer.rotate(angle, 'x');
        viewer.render();
    }
    
    /**
     * Rotate the molecule around the Y axis
     *
     * @param {string} containerId - ID of the container element
     * @param {number} angle - Rotation angle in degrees
     */
    function rotateY(containerId, angle) {
        const viewer = viewers[containerId];
        if (!viewer) return;
        
        viewer.rotate(angle, 'y');
        viewer.render();
    }
    
    /**
     * Rotate the molecule around the Z axis
     *
     * @param {string} containerId - ID of the container element
     * @param {number} angle - Rotation angle in degrees
     */
    function rotateZ(containerId, angle) {
        const viewer = viewers[containerId];
        if (!viewer) return;
        
        viewer.rotate(angle, 'z');
        viewer.render();
    }
    
    /**
     * Zoom in the view
     *
     * @param {string} containerId - ID of the container element
     */
    function zoomIn(containerId) {
        const viewer = viewers[containerId];
        if (!viewer) return;
        
        const state = viewerStates[containerId];
        state.zoomLevel = Math.min(state.zoomLevel + 10, 200); // Max 200%
        
        viewer.zoom(1 + ZOOM_FACTOR);
        viewer.render();
    }
    
    /**
     * Zoom out the view
     *
     * @param {string} containerId - ID of the container element
     */
    function zoomOut(containerId) {
        const viewer = viewers[containerId];
        if (!viewer) return;
        
        const state = viewerStates[containerId];
        state.zoomLevel = Math.max(state.zoomLevel - 10, 20); // Min 20%
        
        viewer.zoom(1 - ZOOM_FACTOR);
        viewer.render();
    }
    
    /**
     * Pan the view
     *
     * @param {string} containerId - ID of the container element
     * @param {number} x - X translation
     * @param {number} y - Y translation
     */
    function pan(containerId, x, y) {
        const viewer = viewers[containerId];
        if (!viewer) return;
        
        viewer.translate(x, y);
        viewer.render();
    }
    
    /**
     * Reset the view to default
     *
     * @param {string} containerId - ID of the container element
     */
    function resetView(containerId) {
        const viewer = viewers[containerId];
        if (!viewer) return;
        
        const state = viewerStates[containerId];
        state.zoomLevel = 100;
        
        viewer.zoomTo();
        viewer.render();
    }
    
    /**
     * Toggle spin animation
     *
     * @param {string} containerId - ID of the container element
     * @param {boolean} enabled - Whether to enable or disable spinning
     */
    function toggleSpin(containerId, enabled) {
        const viewer = viewers[containerId];
        if (!viewer) return;
        
        const state = viewerStates[containerId];
        state.isSpinning = enabled;
        
        viewer.spin(enabled);
    }
    
    /**
     * Toggle surface display
     *
     * @param {string} containerId - ID of the container element
     * @param {boolean} show - Whether to show or hide the surface
     */
    function toggleSurface(containerId, show) {
        const viewer = viewers[containerId];
        if (!viewer) return;
        
        const state = viewerStates[containerId];
        state.showSurface = show;
        
        if (show) {
            viewer.addSurface($3Dmol.SurfaceType.VDW, {
                opacity: 0.7,
                color: 'white'
            });
        } else {
            viewer.removeAllSurfaces();
        }
        
        viewer.render();
    }
    
    /**
     * Create a toolbar for the molecular viewer
     * 
     * @param {string} containerId - ID of the container element
     * @param {string} toolbarId - ID of the toolbar element
     * @param {Object} options - Toolbar options
     */
    function createViewerToolbar(containerId, toolbarId, options = {}) {
        const toolbar = document.getElementById(toolbarId);
        if (!toolbar) {
            console.error(`Toolbar element with ID ${toolbarId} not found`);
            return;
        }
        
        // Clear toolbar
        toolbar.innerHTML = '';
        
        // Create style buttons
        const styleGroup = document.createElement('div');
        styleGroup.className = 'btn-group me-2 mb-2';
        styleGroup.setAttribute('role', 'group');
        styleGroup.setAttribute('aria-label', 'Visualization style');
        
        // Style buttons configuration
        const styleButtons = [
            { id: 'stick', icon: 'bi-diagram-3', label: 'Stick' },
            { id: 'ball-and-stick', icon: 'bi-circle-square', label: 'Ball & Stick' },
            { id: 'spacefill', icon: 'bi-circle-fill', label: 'Spacefill' },
            { id: 'line', icon: 'bi-slash-lg', label: 'Line' }
        ];
        
        // Create style buttons
        styleButtons.forEach(btn => {
            const button = document.createElement('button');
            button.type = 'button';
            button.id = `style-${btn.id}-${containerId}`;
            button.className = 'btn btn-outline-primary';
            button.innerHTML = `<i class="bi ${btn.icon}"></i> ${btn.label}`;
            button.setAttribute('aria-label', `Set ${btn.label} visualization style`);
            button.addEventListener('click', () => {
                applyStyle(containerId, btn.id);
                updateActiveButton(button, styleGroup);
            });
            styleGroup.appendChild(button);
        });
        
        toolbar.appendChild(styleGroup);
        
        // Create rotation controls
        const rotationGroup = document.createElement('div');
        rotationGroup.className = 'btn-group me-2 mb-2';
        rotationGroup.setAttribute('role', 'group');
        rotationGroup.setAttribute('aria-label', 'Rotation controls');
        
        // Rotation buttons configuration
        const rotationButtons = [
            { axis: 'x', dir: -1, icon: 'bi-arrow-up', label: 'Rotate Up' },
            { axis: 'x', dir: 1, icon: 'bi-arrow-down', label: 'Rotate Down' },
            { axis: 'y', dir: -1, icon: 'bi-arrow-left', label: 'Rotate Left' },
            { axis: 'y', dir: 1, icon: 'bi-arrow-right', label: 'Rotate Right' }
        ];
        
        // Create rotation buttons
        rotationButtons.forEach(btn => {
            const button = document.createElement('button');
            button.type = 'button';
            button.className = 'btn btn-outline-secondary';
            button.innerHTML = `<i class="${btn.icon}"></i>`;
            button.setAttribute('aria-label', btn.label);
            button.addEventListener('click', () => {
                if (btn.axis === 'x') {
                    rotateX(containerId, btn.dir * ROTATION_SPEED * 10);
                } else if (btn.axis === 'y') {
                    rotateY(containerId, btn.dir * ROTATION_SPEED * 10);
                }
            });
            rotationGroup.appendChild(button);
        });
        
        toolbar.appendChild(rotationGroup);
        
        // Create zoom controls
        const zoomGroup = document.createElement('div');
        zoomGroup.className = 'btn-group me-2 mb-2';
        zoomGroup.setAttribute('role', 'group');
        zoomGroup.setAttribute('aria-label', 'Zoom controls');
        
        // Zoom in button
        const zoomInBtn = document.createElement('button');
        zoomInBtn.type = 'button';
        zoomInBtn.className = 'btn btn-outline-secondary';
        zoomInBtn.innerHTML = '<i class="bi bi-zoom-in"></i>';
        zoomInBtn.setAttribute('aria-label', 'Zoom In');
        zoomInBtn.addEventListener('click', () => zoomIn(containerId));
        zoomGroup.appendChild(zoomInBtn);
        
        // Zoom out button
        const zoomOutBtn = document.createElement('button');
        zoomOutBtn.type = 'button';
        zoomOutBtn.className = 'btn btn-outline-secondary';
        zoomOutBtn.innerHTML = '<i class="bi bi-zoom-out"></i>';
        zoomOutBtn.setAttribute('aria-label', 'Zoom Out');
        zoomOutBtn.addEventListener('click', () => zoomOut(containerId));
        zoomGroup.appendChild(zoomOutBtn);
        
        toolbar.appendChild(zoomGroup);
        
        // Create surface toggle
        const surfaceGroup = document.createElement('div');
        surfaceGroup.className = 'btn-group me-2 mb-2';
        surfaceGroup.setAttribute('role', 'group');
        surfaceGroup.setAttribute('aria-label', 'Surface options');
        
        // Surface toggle button
        const surfaceBtn = document.createElement('button');
        surfaceBtn.type = 'button';
        surfaceBtn.className = 'btn btn-outline-secondary';
        surfaceBtn.innerHTML = '<i class="bi bi-egg"></i> Surface';
        surfaceBtn.setAttribute('aria-label', 'Toggle Surface');
        surfaceBtn.dataset.state = 'off';
        surfaceBtn.addEventListener('click', () => {
            if (surfaceBtn.dataset.state === 'off') {
                toggleSurface(containerId, true);
                surfaceBtn.dataset.state = 'on';
                surfaceBtn.classList.add('active');
            } else {
                toggleSurface(containerId, false);
                surfaceBtn.dataset.state = 'off';
                surfaceBtn.classList.remove('active');
            }
        });
        surfaceGroup.appendChild(surfaceBtn);
        
        toolbar.appendChild(surfaceGroup);
        
        // Create spin toggle
        const spinGroup = document.createElement('div');
        spinGroup.className = 'btn-group me-2 mb-2';
        spinGroup.setAttribute('role', 'group');
        spinGroup.setAttribute('aria-label', 'Animation options');
        
        // Spin toggle button
        const spinBtn = document.createElement('button');
        spinBtn.type = 'button';
        spinBtn.className = 'btn btn-outline-secondary';
        spinBtn.innerHTML = '<i class="bi bi-arrow-repeat"></i> Spin';
        spinBtn.setAttribute('aria-label', 'Toggle Spin Animation');
        spinBtn.dataset.state = 'off';
        spinBtn.addEventListener('click', () => {
            if (spinBtn.dataset.state === 'off') {
                toggleSpin(containerId, true);
                spinBtn.dataset.state = 'on';
                spinBtn.classList.add('active');
            } else {
                toggleSpin(containerId, false);
                spinBtn.dataset.state = 'off';
                spinBtn.classList.remove('active');
            }
        });
        spinGroup.appendChild(spinBtn);
        
        toolbar.appendChild(spinGroup);
        
        // Create reset view button
        const resetGroup = document.createElement('div');
        resetGroup.className = 'btn-group mb-2';
        resetGroup.setAttribute('role', 'group');
        resetGroup.setAttribute('aria-label', 'Reset view');
        
        // Reset button
        const resetBtn = document.createElement('button');
        resetBtn.type = 'button';
        resetBtn.className = 'btn btn-outline-secondary';
        resetBtn.innerHTML = '<i class="bi bi-arrow-counterclockwise"></i> Reset';
        resetBtn.setAttribute('aria-label', 'Reset View');
        resetBtn.addEventListener('click', () => resetView(containerId));
        resetGroup.appendChild(resetBtn);
        
        toolbar.appendChild(resetGroup);
        
        // Set initial active style button
        const initialStyleBtn = document.getElementById(`style-stick-${containerId}`);
        if (initialStyleBtn) {
            updateActiveButton(initialStyleBtn, styleGroup);
        }
    }
    
    /**
     * Update active button in a button group
     *
     * @param {HTMLElement} activeButton - Button to set as active
     * @param {HTMLElement} group - Button group container
     */
    function updateActiveButton(activeButton, group) {
        // Remove active class from all buttons in the group
        const buttons = group.querySelectorAll('button');
        buttons.forEach(button => {
            button.classList.remove('active');
        });
        
        // Add active class to the selected button
        activeButton.classList.add('active');
    }
    
    // Public API
    return {
        initViewer,
        loadMoleculeFromSmiles,
        loadMoleculeFromPDB,
        loadMoleculeFromMOL,
        createViewerToolbar,
        // New rotation and view controls
        rotateX,
        rotateY,
        rotateZ,
        zoomIn,
        zoomOut,
        pan,
        resetView,
        // Style and visualization controls
        applyStyle,
        toggleSpin,
        toggleSurface
    };
})();

// Export the module
if (typeof module !== 'undefined' && module.exports) {
    module.exports = MolecularViewer;
} else {
    window.MolecularViewer = MolecularViewer;
}
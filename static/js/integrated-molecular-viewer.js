/**
 * CryoProtect Analyzer - Integrated Molecular Viewer
 * 
 * This module provides an integrated interface for 2D and 3D visualization of molecules,
 * combining the functionality of RDKit.js (2D) and 3Dmol.js (3D).
 */

const IntegratedMolecularViewer = (function() {
    // Cache for viewers
    const viewers = {};
    
    // Default options
    const defaultOptions = {
        container2D: null,
        container3D: null,
        toolbar: null,
        width: 400,
        height: 400,
        backgroundColor: 'white',
        style: 'stick',
        colorScheme: 'cyanCarbon',
        showLabels: false,
        showHydrogens: true,
        showSurface: false,
        spin: false,
        measurementMode: false
    };
    
    // Current state
    const state = {
        currentMolecule: null,
        currentFormat: null,
        selectedAtoms: [],
        measurements: [],
        viewType: '3D', // '2D' or '3D'
        isInitialized: false
    };
    
    /**
     * Initialize the integrated molecular viewer
     * 
     * @param {Object} options - Viewer options
     * @returns {Object} - Viewer instance
     */
    function init(options = {}) {
        // Merge options with defaults
        const viewerOptions = { ...defaultOptions, ...options };
        
        // Validate required options
        if (!viewerOptions.container2D || !viewerOptions.container3D) {
            console.error('Both container2D and container3D are required');
            return null;
        }
        
        // Get container elements
        const container2D = typeof viewerOptions.container2D === 'string' 
            ? document.getElementById(viewerOptions.container2D) 
            : viewerOptions.container2D;
            
        const container3D = typeof viewerOptions.container3D === 'string' 
            ? document.getElementById(viewerOptions.container3D) 
            : viewerOptions.container3D;
        
        if (!container2D || !container3D) {
            console.error('Container elements not found');
            return null;
        }
        
        // Initialize 3D viewer
        const viewer3D = $3Dmol.createViewer(container3D, {
            backgroundColor: viewerOptions.backgroundColor,
            antialias: true,
            width: viewerOptions.width,
            height: viewerOptions.height
        });
        
        // Store viewer in cache
        const viewerId = Math.random().toString(36).substring(2, 15);
        viewers[viewerId] = {
            options: viewerOptions,
            viewer3D: viewer3D,
            container2D: container2D,
            container3D: container3D
        };
        
        // Create toolbar if specified
        if (viewerOptions.toolbar) {
            createToolbar(viewerId, viewerOptions.toolbar);
        }
        
        // Set state
        state.isInitialized = true;
        
        return {
            id: viewerId,
            loadMolecule: (data, format) => loadMolecule(viewerId, data, format),
            setViewType: (type) => setViewType(viewerId, type),
            setStyle: (style) => setStyle(viewerId, style),
            toggleSurface: () => toggleSurface(viewerId),
            toggleSpin: () => toggleSpin(viewerId),
            toggleMeasurementMode: () => toggleMeasurementMode(viewerId),
            resetView: () => resetView(viewerId),
            exportImage: (type = 'png') => exportImage(viewerId, type),
            dispose: () => dispose(viewerId)
        };
    }
    
    /**
     * Load a molecule into the viewer
     * 
     * @param {string} viewerId - ID of the viewer
     * @param {string} data - Molecular data (SMILES, MOL, PDB)
     * @param {string} format - Data format ('smiles', 'mol', 'pdb')
     */
    async function loadMolecule(viewerId, data, format) {
        const viewer = viewers[viewerId];
        if (!viewer) {
            console.error(`Viewer with ID ${viewerId} not found`);
            return;
        }
        
        // Update state
        state.currentMolecule = data;
        state.currentFormat = format;
        state.selectedAtoms = [];
        state.measurements = [];
        
        try {
            // Load 2D visualization
            await load2D(viewerId, data, format);
            
            // Load 3D visualization
            await load3D(viewerId, data, format);
            
            // Set view type based on current state
            setViewType(viewerId, state.viewType);
        } catch (error) {
            console.error('Error loading molecule:', error);
            
            // Show error message in containers
            viewer.container2D.innerHTML = `<div class="alert alert-danger">Error loading 2D model: ${error.message}</div>`;
            viewer.container3D.innerHTML = `<div class="alert alert-danger">Error loading 3D model: ${error.message}</div>`;
        }
    }
    
    /**
     * Load a 2D visualization of the molecule
     * 
     * @param {string} viewerId - ID of the viewer
     * @param {string} data - Molecular data
     * @param {string} format - Data format
     */
    async function load2D(viewerId, data, format) {
        const viewer = viewers[viewerId];
        if (!viewer) return;
        
        // Convert format if needed
        let inputFormat = format;
        if (format === 'pdb') {
            inputFormat = 'smiles'; // We'll need to convert PDB to SMILES for 2D visualization
            // This would require a server-side conversion in a real implementation
            // For now, we'll just show an error
            viewer.container2D.innerHTML = `<div class="alert alert-warning">2D visualization of PDB files is not supported yet</div>`;
            return;
        }
        
        try {
            // Use RDKit API to generate 2D visualization
            const visualization = await window.RDKitAPI.generateVisualization(
                data, 
                inputFormat, 
                viewer.options.width, 
                viewer.options.height,
                state.selectedAtoms.length > 0 ? state.selectedAtoms : null
            );
            
            // Display the visualization
            window.RDKitAPI.displayMolecularVisualization(visualization.svg, viewer.container2D);
            
            // Add click handler for atom selection
            const svgElement = viewer.container2D.querySelector('svg');
            if (svgElement) {
                // Add event listeners for atom selection
                const atoms = svgElement.querySelectorAll('.atom');
                atoms.forEach(atom => {
                    atom.addEventListener('click', (event) => {
                        const atomId = parseInt(atom.getAttribute('data-atom-id'));
                        if (!isNaN(atomId)) {
                            toggleAtomSelection(viewerId, atomId);
                        }
                    });
                });
            }
        } catch (error) {
            console.error('Error loading 2D visualization:', error);
            viewer.container2D.innerHTML = `<div class="alert alert-danger">Error loading 2D model: ${error.message}</div>`;
        }
    }
    
    /**
     * Load a 3D visualization of the molecule
     * 
     * @param {string} viewerId - ID of the viewer
     * @param {string} data - Molecular data
     * @param {string} format - Data format
     */
    async function load3D(viewerId, data, format) {
        const viewer = viewers[viewerId];
        if (!viewer) return;
        
        // Clear viewer
        viewer.viewer3D.clear();
        
        try {
            if (format === 'smiles') {
                // Convert SMILES to 3D model
                const model = $3Dmol.createModelFrom(data, "smi", {}, function(model) {
                    // Add model to viewer
                    viewer.viewer3D.addModel(model);
                    
                    // Apply style
                    applyStyle(viewerId);
                    
                    // Add surface if requested
                    if (viewer.options.showSurface) {
                        addSurface(viewerId);
                    }
                    
                    // Highlight selected atoms
                    highlightSelectedAtoms(viewerId);
                    
                    // Show measurements
                    showMeasurements(viewerId);
                    
                    // Zoom to fit
                    viewer.viewer3D.zoomTo();
                    
                    // Start spinning if enabled
                    if (viewer.options.spin) {
                        viewer.viewer3D.spin(true);
                    }
                    
                    // Render
                    viewer.viewer3D.render();
                });
            } else if (format === 'mol') {
                // Add model to viewer
                const model = viewer.viewer3D.addModel(data, "mol");
                
                // Apply style
                applyStyle(viewerId);
                
                // Add surface if requested
                if (viewer.options.showSurface) {
                    addSurface(viewerId);
                }
                
                // Highlight selected atoms
                highlightSelectedAtoms(viewerId);
                
                // Show measurements
                showMeasurements(viewerId);
                
                // Zoom to fit
                viewer.viewer3D.zoomTo();
                
                // Start spinning if enabled
                if (viewer.options.spin) {
                    viewer.viewer3D.spin(true);
                }
                
                // Render
                viewer.viewer3D.render();
            } else if (format === 'pdb') {
                // Load PDB
                $3Dmol.download(`pdb:${data}`, viewer.viewer3D, {}, function() {
                    // Apply style
                    applyStyle(viewerId);
                    
                    // Add surface if requested
                    if (viewer.options.showSurface) {
                        addSurface(viewerId);
                    }
                    
                    // Highlight selected atoms
                    highlightSelectedAtoms(viewerId);
                    
                    // Show measurements
                    showMeasurements(viewerId);
                    
                    // Zoom to fit
                    viewer.viewer3D.zoomTo();
                    
                    // Start spinning if enabled
                    if (viewer.options.spin) {
                        viewer.viewer3D.spin(true);
                    }
                    
                    // Render
                    viewer.viewer3D.render();
                });
            } else {
                throw new Error(`Unsupported format: ${format}`);
            }
            
            // Add click handler for atom selection
            viewer.viewer3D.setClickable({}, true, function(atom) {
                toggleAtomSelection(viewerId, atom.serial);
            });
            
            // Add hover handler for atom highlighting
            viewer.viewer3D.setHoverable({}, true, function(atom, viewer, event, isHovered) {
                if (isHovered) {
                    viewer.addLabel(atom.resn + " " + atom.atom, {position: atom, backgroundColor: 'mintcream', fontColor: 'black'});
                } else {
                    viewer.removeAllLabels();
                }
                viewer.render();
            });
        } catch (error) {
            console.error('Error loading 3D visualization:', error);
            viewer.container3D.innerHTML = `<div class="alert alert-danger">Error loading 3D model: ${error.message}</div>`;
        }
    }
    
    /**
     * Apply the current style to the 3D model
     * 
     * @param {string} viewerId - ID of the viewer
     */
    function applyStyle(viewerId) {
        const viewer = viewers[viewerId];
        if (!viewer) return;
        
        const style = viewer.options.style;
        const colorScheme = viewer.options.colorScheme;
        
        if (style === 'stick') {
            viewer.viewer3D.setStyle({}, {stick: {radius: 0.15, colorscheme: colorScheme}});
        } else if (style === 'sphere') {
            viewer.viewer3D.setStyle({}, {sphere: {scale: 0.25, colorscheme: colorScheme}});
        } else if (style === 'cartoon') {
            viewer.viewer3D.setStyle({}, {cartoon: {colorscheme: colorScheme}});
        } else if (style === 'line') {
            viewer.viewer3D.setStyle({}, {line: {colorscheme: colorScheme}});
        } else {
            viewer.viewer3D.setStyle({}, {stick: {radius: 0.15, colorscheme: colorScheme}});
        }
        
        // Add labels if requested
        if (viewer.options.showLabels) {
            const atoms = viewer.viewer3D.getModel().selectedAtoms({});
            viewer.viewer3D.addLabels(atoms, {font: 'Arial', fontSize: 12, showBackground: false});
        }
        
        viewer.viewer3D.render();
    }
    
    /**
     * Add a surface to the 3D model
     * 
     * @param {string} viewerId - ID of the viewer
     */
    function addSurface(viewerId) {
        const viewer = viewers[viewerId];
        if (!viewer) return;
        
        viewer.viewer3D.addSurface($3Dmol.SurfaceType.VDW, {
            opacity: 0.7,
            color: 'white'
        });
    }
    
    /**
     * Toggle the surface display
     * 
     * @param {string} viewerId - ID of the viewer
     */
    function toggleSurface(viewerId) {
        const viewer = viewers[viewerId];
        if (!viewer) return;
        
        viewer.options.showSurface = !viewer.options.showSurface;
        
        if (viewer.options.showSurface) {
            addSurface(viewerId);
        } else {
            viewer.viewer3D.removeAllSurfaces();
        }
        
        viewer.viewer3D.render();
        
        // Update toolbar button if exists
        const surfaceBtn = document.querySelector(`#${viewer.options.toolbar} .surface-btn`);
        if (surfaceBtn) {
            if (viewer.options.showSurface) {
                surfaceBtn.classList.add('active');
            } else {
                surfaceBtn.classList.remove('active');
            }
        }
    }
    
    /**
     * Toggle molecule spinning
     * 
     * @param {string} viewerId - ID of the viewer
     */
    function toggleSpin(viewerId) {
        const viewer = viewers[viewerId];
        if (!viewer) return;
        
        viewer.options.spin = !viewer.options.spin;
        viewer.viewer3D.spin(viewer.options.spin);
        
        // Update toolbar button if exists
        const spinBtn = document.querySelector(`#${viewer.options.toolbar} .spin-btn`);
        if (spinBtn) {
            if (viewer.options.spin) {
                spinBtn.classList.add('active');
            } else {
                spinBtn.classList.remove('active');
            }
        }
    }
    
    /**
     * Toggle measurement mode
     * 
     * @param {string} viewerId - ID of the viewer
     */
    function toggleMeasurementMode(viewerId) {
        const viewer = viewers[viewerId];
        if (!viewer) return;
        
        viewer.options.measurementMode = !viewer.options.measurementMode;
        
        // Update toolbar button if exists
        const measureBtn = document.querySelector(`#${viewer.options.toolbar} .measure-btn`);
        if (measureBtn) {
            if (viewer.options.measurementMode) {
                measureBtn.classList.add('active');
            } else {
                measureBtn.classList.remove('active');
            }
        }
    }
    
    /**
     * Reset the view
     * 
     * @param {string} viewerId - ID of the viewer
     */
    function resetView(viewerId) {
        const viewer = viewers[viewerId];
        if (!viewer) return;
        
        viewer.viewer3D.zoomTo();
        viewer.viewer3D.render();
    }
    
    /**
     * Toggle atom selection
     * 
     * @param {string} viewerId - ID of the viewer
     * @param {number} atomId - Atom ID
     */
    function toggleAtomSelection(viewerId, atomId) {
        const viewer = viewers[viewerId];
        if (!viewer) return;
        
        // Check if atom is already selected
        const index = state.selectedAtoms.indexOf(atomId);
        
        if (index === -1) {
            // Add to selection
            state.selectedAtoms.push(atomId);
        } else {
            // Remove from selection
            state.selectedAtoms.splice(index, 1);
        }
        
        // If in measurement mode and we have 2 or more atoms selected
        if (viewer.options.measurementMode && state.selectedAtoms.length >= 2) {
            // Get the last two selected atoms
            const atom1 = state.selectedAtoms[state.selectedAtoms.length - 2];
            const atom2 = state.selectedAtoms[state.selectedAtoms.length - 1];
            
            // Add measurement
            addMeasurement(viewerId, atom1, atom2);
        }
        
        // Update visualizations
        highlightSelectedAtoms(viewerId);
        
        // Reload 2D view to show selection
        load2D(viewerId, state.currentMolecule, state.currentFormat);
    }
    
    /**
     * Highlight selected atoms in the 3D view
     * 
     * @param {string} viewerId - ID of the viewer
     */
    function highlightSelectedAtoms(viewerId) {
        const viewer = viewers[viewerId];
        if (!viewer || state.selectedAtoms.length === 0) return;
        
        // Reset styles
        applyStyle(viewerId);
        
        // Highlight selected atoms
        viewer.viewer3D.setStyle({serial: state.selectedAtoms}, {
            stick: {radius: 0.2, colorscheme: 'greenCarbon'},
            sphere: {scale: 0.5, colorscheme: 'greenCarbon'}
        });
        
        viewer.viewer3D.render();
    }
    
    /**
     * Add a measurement between two atoms
     * 
     * @param {string} viewerId - ID of the viewer
     * @param {number} atom1 - First atom ID
     * @param {number} atom2 - Second atom ID
     */
    function addMeasurement(viewerId, atom1, atom2) {
        const viewer = viewers[viewerId];
        if (!viewer) return;
        
        // Add to measurements list
        state.measurements.push({atom1, atom2});
        
        // Show measurements
        showMeasurements(viewerId);
    }
    
    /**
     * Show all measurements in the 3D view
     * 
     * @param {string} viewerId - ID of the viewer
     */
    function showMeasurements(viewerId) {
        const viewer = viewers[viewerId];
        if (!viewer || state.measurements.length === 0) return;
        
        // Remove existing measurements
        viewer.viewer3D.removeAllShapes();
        
        // Add each measurement
        state.measurements.forEach(measurement => {
            viewer.viewer3D.addLabel("Distance", {
                position: {x: 0, y: 0, z: 0},
                useScreen: false,
                backgroundColor: "black",
                fontColor: "white",
                backgroundOpacity: 0.5,
                fontSize: 12,
                showBackground: true
            });
            
            viewer.viewer3D.addLine({
                start: {serial: measurement.atom1},
                end: {serial: measurement.atom2},
                color: "yellow",
                dashed: true
            });
        });
        
        viewer.viewer3D.render();
    }
    
    /**
     * Set the view type (2D or 3D)
     * 
     * @param {string} viewerId - ID of the viewer
     * @param {string} type - View type ('2D' or '3D')
     */
    function setViewType(viewerId, type) {
        const viewer = viewers[viewerId];
        if (!viewer) return;
        
        // Store previous view type for comparison
        const previousViewType = state.viewType;
        
        // Update global state
        state.viewType = type;
        
        // Update viewer-specific state if needed
        if (!viewer.viewType) {
            viewer.viewType = type;
        } else {
            viewer.viewType = type;
        }
        
        // Update container visibility
        if (type === '2D') {
            viewer.container2D.style.display = 'block';
            viewer.container3D.style.display = 'none';
        } else {
            viewer.container2D.style.display = 'none';
            viewer.container3D.style.display = 'block';
            
            // If switching from 2D to 3D, ensure 3D view is properly rendered
            if (previousViewType === '2D' && viewer.viewer3D) {
                viewer.viewer3D.resize();
                viewer.viewer3D.render();
            }
        }
        
        // Update toolbar buttons if exists
        const view2DBtn = document.querySelector(`#${viewer.options.toolbar} .view-2d-btn`);
        const view3DBtn = document.querySelector(`#${viewer.options.toolbar} .view-3d-btn`);
        
        if (view2DBtn && view3DBtn) {
            if (type === '2D') {
                view2DBtn.classList.add('active');
                view3DBtn.classList.remove('active');
            } else {
                view2DBtn.classList.remove('active');
                view3DBtn.classList.add('active');
            }
        }
        
        // Dispatch a custom event to notify other components about the view change
        const event = new CustomEvent('viewTypeChanged', {
            detail: { viewerId, type, previousType: previousViewType }
        });
        document.dispatchEvent(event);
    }
    
    /**
     * Set the visualization style
     * 
     * @param {string} viewerId - ID of the viewer
     * @param {string} style - Style name ('stick', 'sphere', 'cartoon', 'line')
     */
    function setStyle(viewerId, style) {
        const viewer = viewers[viewerId];
        if (!viewer) return;
        
        viewer.options.style = style;
        applyStyle(viewerId);
        
        // Update toolbar buttons if exists
        const styleButtons = document.querySelectorAll(`#${viewer.options.toolbar} .style-btn`);
        styleButtons.forEach(btn => {
            if (btn.dataset.style === style) {
                btn.classList.add('active');
            } else {
                btn.classList.remove('active');
            }
        });
    }
    
    /**
     * Create a toolbar for the molecular viewer
     * 
     * @param {string} viewerId - ID of the viewer
     * @param {string|HTMLElement} toolbarId - ID of the toolbar element or the element itself
     */
    function createToolbar(viewerId, toolbarId) {
        const viewer = viewers[viewerId];
        if (!viewer) return;
        
        const toolbar = typeof toolbarId === 'string' 
            ? document.getElementById(toolbarId) 
            : toolbarId;
            
        if (!toolbar) {
            console.error(`Toolbar element with ID ${toolbarId} not found`);
            return;
        }
        
        // Clear toolbar
        toolbar.innerHTML = '';
        
        // Create view toggle buttons
        const viewGroup = document.createElement('div');
        viewGroup.className = 'btn-group me-2';
        viewGroup.setAttribute('role', 'group');
        viewGroup.setAttribute('aria-label', 'View type');
        
        // 2D view button
        const view2DBtn = document.createElement('button');
        view2DBtn.type = 'button';
        view2DBtn.className = 'btn btn-outline-primary view-2d-btn';
        view2DBtn.innerHTML = '<i class="bi bi-square"></i> 2D';
        view2DBtn.addEventListener('click', () => setViewType(viewerId, '2D'));
        viewGroup.appendChild(view2DBtn);
        
        // 3D view button
        const view3DBtn = document.createElement('button');
        view3DBtn.type = 'button';
        view3DBtn.className = 'btn btn-outline-primary view-3d-btn active';
        view3DBtn.innerHTML = '<i class="bi bi-cube"></i> 3D';
        view3DBtn.addEventListener('click', () => setViewType(viewerId, '3D'));
        viewGroup.appendChild(view3DBtn);
        
        toolbar.appendChild(viewGroup);
        
        // Create style buttons
        const styleGroup = document.createElement('div');
        styleGroup.className = 'btn-group me-2';
        styleGroup.setAttribute('role', 'group');
        styleGroup.setAttribute('aria-label', 'Visualization style');
        
        // Stick style button
        const stickBtn = document.createElement('button');
        stickBtn.type = 'button';
        stickBtn.className = 'btn btn-outline-primary style-btn active';
        stickBtn.dataset.style = 'stick';
        stickBtn.innerHTML = '<i class="bi bi-diagram-3"></i> Stick';
        stickBtn.addEventListener('click', () => setStyle(viewerId, 'stick'));
        styleGroup.appendChild(stickBtn);
        
        // Sphere style button
        const sphereBtn = document.createElement('button');
        sphereBtn.type = 'button';
        sphereBtn.className = 'btn btn-outline-primary style-btn';
        sphereBtn.dataset.style = 'sphere';
        sphereBtn.innerHTML = '<i class="bi bi-circle"></i> Sphere';
        sphereBtn.addEventListener('click', () => setStyle(viewerId, 'sphere'));
        styleGroup.appendChild(sphereBtn);
        
        // Line style button
        const lineBtn = document.createElement('button');
        lineBtn.type = 'button';
        lineBtn.className = 'btn btn-outline-primary style-btn';
        lineBtn.dataset.style = 'line';
        lineBtn.innerHTML = '<i class="bi bi-slash-lg"></i> Line';
        lineBtn.addEventListener('click', () => setStyle(viewerId, 'line'));
        styleGroup.appendChild(lineBtn);
        
        toolbar.appendChild(styleGroup);
        
        // Create surface toggle
        const surfaceBtn = document.createElement('button');
        surfaceBtn.type = 'button';
        surfaceBtn.className = 'btn btn-outline-secondary surface-btn me-2';
        surfaceBtn.innerHTML = '<i class="bi bi-egg"></i> Surface';
        surfaceBtn.addEventListener('click', () => toggleSurface(viewerId));
        toolbar.appendChild(surfaceBtn);
        
        // Create measurement toggle
        const measureBtn = document.createElement('button');
        measureBtn.type = 'button';
        measureBtn.className = 'btn btn-outline-secondary measure-btn me-2';
        measureBtn.innerHTML = '<i class="bi bi-rulers"></i> Measure';
        measureBtn.addEventListener('click', () => toggleMeasurementMode(viewerId));
        toolbar.appendChild(measureBtn);
        
        // Create spin toggle
        const spinBtn = document.createElement('button');
        spinBtn.type = 'button';
        spinBtn.className = 'btn btn-outline-secondary spin-btn me-2';
        spinBtn.innerHTML = '<i class="bi bi-arrow-repeat"></i> Spin';
        spinBtn.addEventListener('click', () => toggleSpin(viewerId));
        toolbar.appendChild(spinBtn);
        
        // Create reset view button
        const resetBtn = document.createElement('button');
        resetBtn.type = 'button';
        resetBtn.className = 'btn btn-outline-secondary';
        resetBtn.innerHTML = '<i class="bi bi-arrow-counterclockwise"></i> Reset';
        resetBtn.addEventListener('click', () => resetView(viewerId));
        toolbar.appendChild(resetBtn);
    }
    
    /**
     * Dispose of a viewer instance
     * 
     * @param {string} viewerId - ID of the viewer
     */
    function dispose(viewerId) {
        const viewer = viewers[viewerId];
        if (!viewer) return;
        
        // Clean up 3D viewer
        if (viewer.viewer3D) {
            viewer.viewer3D.removeAllModels();
            viewer.viewer3D.removeAllLabels();
            viewer.viewer3D.removeAllShapes();
            viewer.viewer3D.removeAllSurfaces();
        }
        
        // Clean up containers
        if (viewer.container2D) {
            viewer.container2D.innerHTML = '';
        }
        
        if (viewer.container3D) {
            viewer.container3D.innerHTML = '';
        }
        
        // Remove from cache
        delete viewers[viewerId];
    }
    
    // Public API
    return {
        init,
        dispose
    };
})();

// Export the module
if (typeof module !== 'undefined' && module.exports) {
    module.exports = IntegratedMolecularViewer;
} else {
    window.IntegratedMolecularViewer = IntegratedMolecularViewer;
}

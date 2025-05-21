/**
 * MoleculeViewer3D Component
 * 3D molecular structure viewer for SMILES strings
 */
import React, { useEffect, useRef } from 'react';
import { Card } from '../ui/card';

export function MoleculeViewer3D({
  smiles,
  height = 300,
  name,
  style = 'stick',
  spin = false,
  backgroundColor = 'transparent',
}) {
  const containerRef = useRef(null);
  const viewerRef = useRef(null);

  useEffect(() => {
    let isMounted = true;
    
    const loadMolecule = async () => {
      if (!containerRef.current || typeof window === 'undefined') return;
      
      try {
        // Dynamically import 3DMol.js to ensure it's only loaded in the browser
        const $3Dmol = await import('3dmol');
        
        // Clear the container
        containerRef.current.innerHTML = '';
        
        // Create viewer
        const viewer = $3Dmol.createViewer(
          containerRef.current,
          { 
            backgroundColor: backgroundColor,
            antialias: true,
            height: height,
            width: '100%',
          }
        );
        viewerRef.current = viewer;
        
        if (!smiles) {
          // If no SMILES is provided, show a message
          renderPlaceholder();
          return;
        }
        
        // Add the model from SMILES data
        const model = viewer.addModel();
        model.addMolData(smiles, 'smi');
        
        // Generate 3D coordinates
        $3Dmol.generate3DStructure(model);
        
        // Set style based on prop
        if (style === 'stick') {
          model.setStyle({}, { stick: { radius: 0.15 } });
        } else if (style === 'sphere') {
          model.setStyle({}, { sphere: { scale: 0.25 } });
        } else if (style === 'cartoon') {
          model.setStyle({}, { cartoon: { color: 'spectrum' } });
        } else if (style === 'surface') {
          model.setStyle({}, { surface: { opacity: 0.8 } });
        }
        
        // Set up spin if enabled
        if (spin) {
          viewer.spin(true);
        }
        
        // Zoom to fit the molecule
        viewer.zoomTo();
        
        // Render the scene
        viewer.render();
        
        // Add labels if name is provided
        if (name) {
          // Show the molecule name in the top left corner
          const labelSpec = {
            position: { x: -10, y: 10, z: 0 },
            backgroundColor: "rgba(0, 0, 0, 0.5)",
            fontColor: "white",
            font: "sans-serif",
            fontSize: 14,
            showBackground: true,
            fixedPosition: true,
            alignment: $3Dmol.SpriteAlignment.topLeft,
          };
          
          viewer.addLabel(name, labelSpec);
          viewer.render();
        }
      } catch (error) {
        console.error('Failed to render molecule:', error);
        renderPlaceholder(error instanceof Error ? error.message : 'Failed to render molecule');
      }
    };
    
    const renderPlaceholder = (errorMessage) => {
      if (!containerRef.current) return;
      
      const container = containerRef.current;
      container.innerHTML = '';
      
      // Create placeholder
      const placeholder = document.createElement('div');
      placeholder.className = 'flex items-center justify-center w-full h-full bg-muted/30 rounded-md';
      placeholder.style.height = `${height}px`;
      
      const placeholderText = document.createElement('div');
      placeholderText.className = 'text-center p-4';
      
      const titleEl = document.createElement('p');
      titleEl.className = 'font-medium';
      titleEl.textContent = name || 'Molecule Visualization';
      
      const smilesEl = document.createElement('p');
      smilesEl.className = 'text-xs text-muted-foreground mt-1';
      smilesEl.textContent = `SMILES: ${smiles || 'Not provided'}`;
      
      const infoEl = document.createElement('p');
      infoEl.className = 'text-xs text-muted-foreground mt-3';
      
      if (errorMessage) {
        infoEl.textContent = `Error: ${errorMessage}`;
        infoEl.className += ' text-destructive';
      } else if (!smiles) {
        infoEl.textContent = 'No molecule data provided';
      } else {
        infoEl.textContent = 'Loading 3D molecule rendering...';
      }
      
      placeholderText.appendChild(titleEl);
      placeholderText.appendChild(smilesEl);
      placeholderText.appendChild(infoEl);
      placeholder.appendChild(placeholderText);
      
      container.appendChild(placeholder);
    };
    
    // Show a placeholder while loading
    renderPlaceholder();
    
    // Load the molecule
    loadMolecule();
    
    // Cleanup function
    return () => {
      isMounted = false;
      if (viewerRef.current) {
        try {
          // Stop spinning if enabled
          viewerRef.current.spin(false);
          // Remove any event listeners
          viewerRef.current.removeAllLabels();
          viewerRef.current = null;
        } catch (e) {
          console.error('Error cleaning up 3DMol viewer:', e);
        }
      }
      if (containerRef.current) {
        containerRef.current.innerHTML = '';
      }
    };
  }, [smiles, height, name, style, spin, backgroundColor]);

  return (
    <Card className="overflow-hidden">
      <div className="molecule-viewer-container">
        <div 
          ref={containerRef}
          className="w-full h-full" 
          style={{ height: `${height}px` }}
          role="region" 
          aria-label={`3D visualization of ${name || 'molecule'}`}
        />
      </div>
    </Card>
  );
}
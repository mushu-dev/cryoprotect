import React, { useEffect, useRef, useState } from 'react';

/**
 * A 3D molecule viewer component using Mol* library
 * @param {Object} props - Component props
 * @param {string} props.smiles - SMILES notation of the molecule
 * @param {string} props.name - Name of the molecule
 * @param {string} props.pubchem_cid - PubChem CID for direct 3D structure loading
 */
export default function MoleculeViewer3D({ smiles, name, pubchem_cid }) {
  const containerRef = useRef(null);
  const viewerRef = useRef(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState(null);
  const [molData, setMolData] = useState(null);
  const [viewerMode, setViewerMode] = useState('ball-and-stick');
  const [retryCount, setRetryCount] = useState(0);
  const [loadingStatus, setLoadingStatus] = useState('');

  // Load Mol* library dynamically
  useEffect(() => {
    // Check if container ref is available
    if (!containerRef.current) return;

    const loadMolstarLibrary = async () => {
      try {
        setLoadingStatus('Initializing 3D viewer...');
        
        // In a production app, we would dynamically import Mol* here
        // import('molstar').then(({ Viewer }) => { ... })
        // For this minimal implementation, we'll use a placeholder while
        // preparing for the real integration
        console.log('Preparing Mol* integration for', name || 'molecule');
        
        // Simulate loading the Mol* library
        await new Promise(resolve => setTimeout(resolve, 800));
        
        // Set up mock data for now
        setMolData({
          smiles,
          name,
          pubchem_cid,
          loaded: true
        });
        
        setLoadingStatus('');
        setLoading(false);
      } catch (err) {
        console.error('Failed to load Mol* library:', err);
        setError('Failed to initialize 3D viewer. Please try again later.');
        setLoading(false);
        setLoadingStatus('');
      }
    };

    loadMolstarLibrary();

    // Cleanup function to remove any Mol* instances
    return () => {
      if (viewerRef.current) {
        // In real implementation, we would call viewer.dispose() here
        viewerRef.current = null;
      }
    };
  }, []);

  // Load molecule data with retry logic
  useEffect(() => {
    if (!smiles && !pubchem_cid) {
      setError('No molecule data provided');
      setLoading(false);
      return;
    }

    const loadMoleculeStructure = async () => {
      try {
        setLoading(true);
        setError(null);
        
        // Set informative loading message
        if (pubchem_cid) {
          setLoadingStatus(`Loading 3D structure for ${name || 'molecule'} (PubChem CID: ${pubchem_cid})...`);
        } else {
          setLoadingStatus(`Generating 3D structure from SMILES...`);
        }
        
        // Prioritize loading from PubChem CID if available
        if (pubchem_cid) {
          console.log(`Loading 3D structure from PubChem CID: ${pubchem_cid}`);
          // In real implementation, we would fetch from PubChem API or cached data
          // const response = await fetch(`https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${pubchem_cid}/record/JSON`);
          
          // Simulate the API call
          await new Promise(resolve => setTimeout(resolve, 600));
        } else if (smiles) {
          console.log(`Converting SMILES to 3D structure: ${smiles}`);
          // In real implementation, we would use RDKit or similar service
          // For example:
          // const response = await fetch(`${API_ENDPOINTS.RDKIT_API}/convert3d`, {
          //   method: 'POST',
          //   headers: { 'Content-Type': 'application/json' },
          //   body: JSON.stringify({ smiles })
          // });
          
          // Simulate the conversion
          await new Promise(resolve => setTimeout(resolve, 800));
        }
        
        // Update molecule data
        setMolData({
          smiles,
          name,
          pubchem_cid,
          loaded: true
        });
        
        setLoadingStatus('');
        setLoading(false);
      } catch (err) {
        console.error('Failed to load 3D structure:', err);
        
        // Implement exponential backoff for retries
        if (retryCount < 3) {
          const backoffTime = Math.pow(2, retryCount) * 1000; // 1s, 2s, 4s
          setLoadingStatus(`Connection issue. Retrying in ${backoffTime/1000}s... (attempt ${retryCount + 1}/3)`);
          
          setTimeout(() => {
            setRetryCount(prev => prev + 1);
          }, backoffTime);
        } else {
          setError('Failed to load 3D structure after multiple attempts. Please try again later.');
          setLoading(false);
          setLoadingStatus('');
        }
      }
    };

    if (!molData?.loaded || retryCount > 0) {
      loadMoleculeStructure();
    }
  }, [smiles, name, pubchem_cid, retryCount]);

  // Handle viewer controls
  const handleRotate = () => {
    // In a real implementation with Mol* library, this would trigger rotation
    // viewer.canvas3d.requestCameraReset({ durationMs: 1000, snapshot: viewer.canvas3d.camera.getSnapshot() });
    console.log('Rotate molecule');
  };

  const handleZoom = (factor = 1.2) => {
    // In a real implementation with Mol* library, this would trigger zoom
    // viewer.canvas3d.requestCameraZoom({ delta: factor > 1 ? 1 : -1, durationMs: 100 });
    console.log('Zoom molecule by factor', factor);
  };

  const handleReset = () => {
    // In a real implementation with Mol* library, this would reset the view
    // viewer.canvas3d.resetCamera();
    console.log('Reset view');
  };

  const handleViewModeChange = (e) => {
    const newMode = e.target.value;
    setViewerMode(newMode);
    
    // In a real implementation with Mol* library, this would change the representation
    // Example:
    // const visual = viewer.plugin.build();
    // visual.update(component => {
    //   // Remove existing representations
    //   component.representations = [];
    //   
    //   // Add new representation based on mode
    //   switch(newMode) {
    //     case 'ball-and-stick':
    //       component.addRepresentation('ball-and-stick');
    //       break;
    //     case 'space-filling':
    //       component.addRepresentation('spacefill');
    //       break;
    //     case 'wireframe':
    //       component.addRepresentation('wire');
    //       break;
    //     case 'cartoon':
    //       component.addRepresentation('cartoon');
    //       break;
    //     case 'surface':
    //       component.addRepresentation('surface');
    //       break;
    //   }
    // });
    
    console.log('Change view mode to', newMode);
  };

  // Handle retry button click
  const handleRetry = () => {
    setRetryCount(0); // Reset retry count to trigger a new attempt
    setError(null);
  };

  return (
    <div className="molecule-viewer-3d">
      <div 
        className="molecule-viewer-container" 
        ref={containerRef}
        role="application"
        aria-label={`3D visualization of ${name || 'molecule'}`}
      >
        {loading ? (
          <div className="molecule-viewer-loading">
            <div className="viewer-spinner" aria-hidden="true"></div>
            <p>{loadingStatus || 'Loading 3D structure...'}</p>
          </div>
        ) : error ? (
          <div className="molecule-viewer-error">
            <p>{error}</p>
            <button 
              className="retry-button"
              onClick={handleRetry}
              aria-label="Retry loading 3D structure"
            >
              Try Again
            </button>
          </div>
        ) : (
          <div className="molecule-viewer-canvas">
            {/* This will be replaced with Mol* viewer in the full implementation */}
            <div className="molecule-viewer-placeholder">
              <div className="molecule-3d-model" aria-hidden="true">
                <div className="atom atom-1"></div>
                <div className="atom atom-2"></div>
                <div className="atom atom-3"></div>
                <div className="atom atom-4"></div>
                <div className="bond bond-1"></div>
                <div className="bond bond-2"></div>
                <div className="bond bond-3"></div>
              </div>
              <div className="molecule-info">
                <p className="molecule-name">{name || 'Molecule'}</p>
                {pubchem_cid && (
                  <p className="molecule-pubchem">
                    <a 
                      href={`https://pubchem.ncbi.nlm.nih.gov/compound/${pubchem_cid}`}
                      target="_blank"
                      rel="noopener noreferrer"
                      className="pubchem-link"
                      aria-label={`View ${name || 'molecule'} on PubChem (opens in new tab)`}
                    >
                      View on PubChem
                    </a>
                  </p>
                )}
                {smiles && (
                  <p className="molecule-smiles" title={smiles}>
                    {smiles.length > 30 ? `${smiles.substring(0, 30)}...` : smiles}
                  </p>
                )}
                <p className="viewer-mode-info">
                  Viewing mode: <strong>{viewerMode}</strong>
                </p>
              </div>
            </div>
          </div>
        )}
      </div>
      
      <div className="molecule-viewer-controls" role="toolbar" aria-label="Molecule viewer controls">
        <button 
          className="control-button"
          onClick={handleRotate}
          aria-label="Rotate molecule"
          disabled={loading || !!error}
        >
          Rotate
        </button>
        <div className="zoom-controls">
          <button 
            className="control-button zoom-button"
            onClick={() => handleZoom(1.2)}
            aria-label="Zoom in"
            disabled={loading || !!error}
          >
            Zoom +
          </button>
          <button 
            className="control-button zoom-button"
            onClick={() => handleZoom(0.8)}
            aria-label="Zoom out"
            disabled={loading || !!error}
          >
            Zoom -
          </button>
        </div>
        <button 
          className="control-button"
          onClick={handleReset}
          aria-label="Reset view"
          disabled={loading || !!error}
        >
          Reset
        </button>
        <div className="select-wrapper">
          <select 
            className="control-select"
            value={viewerMode}
            onChange={handleViewModeChange}
            aria-label="Change molecule representation style"
            disabled={loading || !!error}
          >
            <option value="ball-and-stick">Ball and Stick</option>
            <option value="space-filling">Space Filling</option>
            <option value="wireframe">Wireframe</option>
            <option value="cartoon">Cartoon</option>
            <option value="surface">Surface</option>
          </select>
        </div>
      </div>
      
      <style jsx>{`
        .molecule-viewer-3d {
          width: 100%;
          margin-top: 15px;
        }
        
        .molecule-viewer-container {
          height: 400px;
          background-color: #f8f9fa;
          border-radius: 8px;
          border: 1px solid #e2e8f0;
          overflow: hidden;
          position: relative;
        }
        
        .molecule-viewer-loading {
          height: 100%;
          display: flex;
          flex-direction: column;
          justify-content: center;
          align-items: center;
          color: #4a5568;
          text-align: center;
          padding: 0 20px;
        }
        
        .viewer-spinner {
          width: 40px;
          height: 40px;
          border: 4px solid rgba(0, 112, 243, 0.2);
          border-left-color: #0070f3;
          border-radius: 50%;
          animation: spin 1s linear infinite;
          margin-bottom: 10px;
        }
        
        @keyframes spin {
          to {
            transform: rotate(360deg);
          }
        }
        
        .molecule-viewer-error {
          height: 100%;
          display: flex;
          flex-direction: column;
          justify-content: center;
          align-items: center;
          color: #e53e3e;
          padding: 0 20px;
          text-align: center;
        }
        
        .retry-button {
          margin-top: 15px;
          background-color: #0070f3;
          color: white;
          border: none;
          border-radius: 4px;
          padding: 8px 16px;
          cursor: pointer;
          font-weight: 500;
          transition: background-color 0.2s;
        }
        
        .retry-button:hover {
          background-color: #0051b3;
        }
        
        .molecule-viewer-canvas {
          height: 100%;
          width: 100%;
          display: flex;
          justify-content: center;
          align-items: center;
        }
        
        .molecule-viewer-placeholder {
          text-align: center;
          display: flex;
          flex-direction: column;
          align-items: center;
        }
        
        .molecule-3d-model {
          width: 180px;
          height: 180px;
          position: relative;
          margin: 0 auto 20px;
          transform-style: preserve-3d;
          animation: rotate 15s infinite linear;
        }
        
        @keyframes rotate {
          from {
            transform: rotateY(0deg) rotateX(15deg);
          }
          to {
            transform: rotateY(360deg) rotateX(15deg);
          }
        }
        
        .atom {
          position: absolute;
          width: 30px;
          height: 30px;
          border-radius: 50%;
          opacity: 0.8;
        }
        
        .atom-1 {
          background-color: #ff4d4d;
          transform: translateZ(60px);
        }
        
        .atom-2 {
          background-color: #4dabff;
          transform: translateZ(-60px);
        }
        
        .atom-3 {
          background-color: #ffcc4d;
          transform: translateX(60px);
        }
        
        .atom-4 {
          background-color: #4dff77;
          transform: translateX(-60px);
        }
        
        .bond {
          position: absolute;
          background-color: #999;
          width: 10px;
          height: 120px;
          top: -45px;
          left: 85px;
        }
        
        .bond-1 {
          transform: rotateY(0deg);
        }
        
        .bond-2 {
          transform: rotateY(60deg);
        }
        
        .bond-3 {
          transform: rotateY(120deg);
        }
        
        .molecule-info {
          display: flex;
          flex-direction: column;
          align-items: center;
          gap: 5px;
        }
        
        .molecule-name {
          font-weight: bold;
          font-size: 1.2rem;
          margin-bottom: 2px;
        }
        
        .molecule-pubchem {
          margin: 0 0 5px 0;
        }
        
        .pubchem-link {
          color: #0070f3;
          text-decoration: none;
          font-size: 0.9rem;
          display: inline-flex;
          align-items: center;
          transition: color 0.2s;
        }
        
        .pubchem-link:hover {
          color: #0051b3;
          text-decoration: underline;
        }
        
        .molecule-smiles {
          font-family: monospace;
          background-color: #edf2f7;
          padding: 5px 10px;
          border-radius: 4px;
          display: inline-block;
          font-size: 0.9rem;
          margin: 5px 0;
          max-width: 100%;
          overflow-x: auto;
        }
        
        .viewer-mode-info {
          font-size: 0.9rem;
          color: #4a5568;
          margin-top: 5px;
        }
        
        .molecule-viewer-controls {
          display: flex;
          justify-content: center;
          gap: 10px;
          margin-top: 15px;
          flex-wrap: wrap;
        }
        
        .control-button {
          padding: 8px 16px;
          background-color: #0070f3;
          color: white;
          border: none;
          border-radius: 4px;
          cursor: pointer;
          font-weight: 500;
          transition: background-color 0.2s;
        }
        
        .control-button:hover:not(:disabled) {
          background-color: #0051b3;
        }
        
        .control-button:disabled {
          background-color: #a0aec0;
          cursor: not-allowed;
        }
        
        .zoom-controls {
          display: flex;
          gap: 4px;
        }
        
        .zoom-button {
          padding: 8px 12px;
        }
        
        .select-wrapper {
          position: relative;
        }
        
        .control-select {
          padding: 8px 16px;
          border: 1px solid #e2e8f0;
          border-radius: 4px;
          background-color: white;
          cursor: pointer;
          appearance: none;
          padding-right: 30px;
        }
        
        .control-select:disabled {
          background-color: #edf2f7;
          color: #a0aec0;
          cursor: not-allowed;
        }
        
        .select-wrapper::after {
          content: '';
          position: absolute;
          right: 12px;
          top: 50%;
          transform: translateY(-50%);
          width: 0;
          height: 0;
          border-left: 5px solid transparent;
          border-right: 5px solid transparent;
          border-top: 5px solid #718096;
          pointer-events: none;
        }
        
        @media (max-width: 768px) {
          .molecule-viewer-container {
            height: 300px;
          }
          
          .molecule-smiles {
            font-size: 0.8rem;
            max-width: 250px;
            white-space: nowrap;
            overflow: hidden;
            text-overflow: ellipsis;
          }
          
          .molecule-viewer-controls {
            flex-direction: column;
            align-items: center;
          }
          
          .control-button, .control-select {
            width: 100%;
            max-width: 250px;
          }
        }
      `}</style>
    </div>
  );
}
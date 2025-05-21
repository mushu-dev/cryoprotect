import React, { useEffect, useRef, useState } from 'react';

/**
 * A simple 3D molecule viewer component
 * @param {Object} props - Component props
 * @param {string} props.smiles - SMILES notation of the molecule
 * @param {string} props.name - Name of the molecule
 */
export default function MoleculeViewer3D({ smiles, name }) {
  const containerRef = useRef(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState(null);
  const [molData, setMolData] = useState(null);

  useEffect(() => {
    if (!smiles) {
      setError('No SMILES data provided');
      setLoading(false);
      return;
    }

    const fetchMolecule = async () => {
      try {
        setLoading(true);
        setError(null);
        
        // In a real app, we would use a library like 3DMol.js or RDKit to render the molecule
        // For this demo, we'll simulate loading the 3D structure
        
        // Simulate API call for 3D structure
        await new Promise(resolve => setTimeout(resolve, 1000));
        
        // Mock data for the molecule visualization
        setMolData({
          smiles,
          name,
          atoms: 'Sample 3D data would be loaded here'
        });
        
        setLoading(false);
      } catch (err) {
        console.error('Failed to load 3D structure:', err);
        setError('Failed to load 3D structure. Please try again later.');
        setLoading(false);
      }
    };

    fetchMolecule();
  }, [smiles, name]);

  return (
    <div className="molecule-viewer-3d">
      <div className="molecule-viewer-container" ref={containerRef}>
        {loading ? (
          <div className="molecule-viewer-loading">
            <div className="viewer-spinner"></div>
            <p>Loading 3D structure...</p>
          </div>
        ) : error ? (
          <div className="molecule-viewer-error">
            <p>{error}</p>
          </div>
        ) : (
          <div className="molecule-viewer-canvas">
            {/* This would be replaced with a real 3D viewer in a full implementation */}
            <div className="molecule-viewer-placeholder">
              <div className="molecule-3d-model">
                <div className="atom atom-1"></div>
                <div className="atom atom-2"></div>
                <div className="atom atom-3"></div>
                <div className="atom atom-4"></div>
                <div className="bond bond-1"></div>
                <div className="bond bond-2"></div>
                <div className="bond bond-3"></div>
              </div>
              <p className="molecule-name">{name || 'Molecule'}</p>
              <p className="molecule-smiles">{smiles}</p>
            </div>
          </div>
        )}
      </div>
      
      <div className="molecule-viewer-controls">
        <button className="control-button">Rotate</button>
        <button className="control-button">Zoom</button>
        <button className="control-button">Reset</button>
        <select className="control-select">
          <option value="ball-and-stick">Ball and Stick</option>
          <option value="space-filling">Space Filling</option>
          <option value="wireframe">Wireframe</option>
        </select>
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
          justify-content: center;
          align-items: center;
          color: #e53e3e;
          padding: 0 20px;
          text-align: center;
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
        
        .molecule-name {
          font-weight: bold;
          font-size: 1.2rem;
          margin-bottom: 5px;
        }
        
        .molecule-smiles {
          font-family: monospace;
          background-color: #edf2f7;
          padding: 5px 10px;
          border-radius: 4px;
          display: inline-block;
          font-size: 0.9rem;
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
        
        .control-button:hover {
          background-color: #0051b3;
        }
        
        .control-select {
          padding: 8px 16px;
          border: 1px solid #e2e8f0;
          border-radius: 4px;
          background-color: white;
          cursor: pointer;
        }
        
        @media (max-width: 768px) {
          .molecule-viewer-container {
            height: 300px;
          }
        }
      `}</style>
    </div>
  );
}
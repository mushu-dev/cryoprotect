import React, { useEffect, useState } from 'react';
import { getMoleculeDepiction } from '../utils/api';

/**
 * A 2D molecule structure viewer component
 * @param {Object} props - Component props
 * @param {string} props.smiles - SMILES notation of the molecule
 * @param {string} props.name - Name of the molecule
 */
export default function MoleculeViewer2D({ smiles, name }) {
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState(null);
  const [svgContent, setSvgContent] = useState('');
  
  useEffect(() => {
    if (!smiles) {
      setError('No SMILES data provided');
      setLoading(false);
      return;
    }

    const fetchMoleculeDepiction = async () => {
      try {
        setLoading(true);
        setError(null);
        
        // In a full implementation, we would get SVG from RDKit or similar service
        const svg = await getMoleculeDepiction(smiles);
        
        if (svg) {
          setSvgContent(svg);
        } else {
          // If no SVG is returned, we'll use a fallback placeholder
          setError('Could not generate molecule structure');
        }
        
        setLoading(false);
      } catch (err) {
        console.error('Failed to load 2D structure:', err);
        setError('Failed to load 2D structure');
        setLoading(false);
      }
    };

    fetchMoleculeDepiction();
  }, [smiles]);

  return (
    <div className="molecule-viewer-2d">
      <h3 className="viewer-title">2D Structure</h3>
      
      <div className="molecule-viewer-container">
        {loading ? (
          <div className="molecule-viewer-loading">
            <div className="viewer-spinner"></div>
            <p>Generating structure...</p>
          </div>
        ) : error ? (
          <div className="molecule-viewer-error">
            {/* Fallback placeholder for when we can't get a proper depiction */}
            <div className="molecule-fallback">
              <div className="molecule-placeholder-structure">
                <div className="atom-circle"></div>
                <div className="bond"></div>
                <div className="atom-circle"></div>
              </div>
              <p className="error-message">{error}</p>
            </div>
          </div>
        ) : (
          <div className="molecule-depiction">
            {svgContent ? (
              <div dangerouslySetInnerHTML={{ __html: svgContent }} />
            ) : (
              <div className="molecule-fallback">
                <div className="molecule-placeholder-structure">
                  <div className="atom-circle"></div>
                  <div className="bond"></div>
                  <div className="atom-circle"></div>
                </div>
              </div>
            )}
            <p className="molecule-formula">{name}</p>
          </div>
        )}
      </div>
      
      <style jsx>{`
        .molecule-viewer-2d {
          width: 100%;
          margin-top: 5px;
        }
        
        .viewer-title {
          margin-bottom: 10px;
          color: #4a5568;
          font-size: 1.1rem;
        }
        
        .molecule-viewer-container {
          min-height: 200px;
          background-color: #f8f9fa;
          border-radius: 8px;
          border: 1px solid #e2e8f0;
          padding: 15px;
          display: flex;
          justify-content: center;
          align-items: center;
        }
        
        .molecule-viewer-loading {
          display: flex;
          flex-direction: column;
          justify-content: center;
          align-items: center;
          color: #4a5568;
        }
        
        .viewer-spinner {
          width: 30px;
          height: 30px;
          border: 3px solid rgba(0, 112, 243, 0.2);
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
          width: 100%;
          text-align: center;
        }
        
        .molecule-fallback {
          display: flex;
          flex-direction: column;
          align-items: center;
          justify-content: center;
        }
        
        .molecule-placeholder-structure {
          display: flex;
          align-items: center;
          margin-bottom: 10px;
        }
        
        .atom-circle {
          width: 30px;
          height: 30px;
          border-radius: 50%;
          background-color: #0070f3;
          opacity: 0.6;
        }
        
        .bond {
          width: 60px;
          height: 6px;
          background-color: #718096;
        }
        
        .error-message {
          color: #e53e3e;
          font-size: 0.9rem;
        }
        
        .molecule-depiction {
          display: flex;
          flex-direction: column;
          align-items: center;
          max-width: 100%;
        }
        
        .molecule-depiction svg {
          max-width: 100%;
          height: auto;
        }
        
        .molecule-formula {
          margin-top: 10px;
          font-family: monospace;
          color: #4a5568;
        }
      `}</style>
    </div>
  );
}
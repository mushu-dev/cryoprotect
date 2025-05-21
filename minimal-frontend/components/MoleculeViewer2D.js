import React, { useEffect, useState, useRef } from 'react';
import { getMoleculeDepiction } from '../utils/api';

/**
 * A 2D molecule structure viewer component with improved accessibility and error handling
 * @param {Object} props - Component props
 * @param {string} props.smiles - SMILES notation of the molecule
 * @param {string} props.name - Name of the molecule
 * @param {string} props.formula - Chemical formula of the molecule
 * @param {string} props.pubchem_cid - PubChem CID for the molecule
 */
export default function MoleculeViewer2D({ smiles, name, formula, pubchem_cid }) {
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState(null);
  const [svgContent, setSvgContent] = useState('');
  const [retryCount, setRetryCount] = useState(0);
  const [loadingStatus, setLoadingStatus] = useState('');
  const containerRef = useRef(null);
  
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
        setLoadingStatus('Generating molecule structure...');
        
        // Get SVG representation from RDKit service via API utility
        const svg = await getMoleculeDepiction(smiles, {
          width: 300,
          height: 200,
          includeMoleculeDetails: true
        });
        
        if (svg) {
          // Process SVG to add ARIA attributes for better accessibility
          // This would be more complex in the actual implementation
          const processedSvg = svg
            .replace('<svg ', '<svg aria-label="Molecular structure" role="img" ')
            .replace(/id="mol/g, `id="mol-${name ? name.toLowerCase().replace(/\s+/g, '-') : 'molecule'}`);
          
          setSvgContent(processedSvg);
        } else {
          // If no SVG is returned, we'll use a fallback placeholder
          throw new Error('Failed to generate molecule structure');
        }
        
        setLoading(false);
        setLoadingStatus('');
      } catch (err) {
        console.error('Failed to load 2D structure:', err);
        
        // Implement retry logic with exponential backoff
        if (retryCount < 3) {
          const backoffTime = Math.pow(2, retryCount) * 1000; // Exponential backoff: 1s, 2s, 4s
          console.log(`Retrying in ${backoffTime/1000}s... (attempt ${retryCount + 1}/3)`);
          setLoadingStatus(`Connection issue. Retrying in ${backoffTime/1000}s... (attempt ${retryCount + 1}/3)`);
          
          setTimeout(() => {
            setRetryCount(prev => prev + 1);
          }, backoffTime);
        } else {
          setError('Failed to load 2D structure after multiple attempts');
          setLoading(false);
          setLoadingStatus('');
        }
      }
    };

    fetchMoleculeDepiction();
  }, [smiles, retryCount]);

  // Handle manual retry
  const handleRetry = () => {
    setRetryCount(0); // Reset retry count to trigger a new attempt
    setError(null);
  };

  // Simplified molecule structure for fallback
  const FallbackMolecule = () => (
    <div className="molecule-placeholder-structure" aria-hidden="true">
      <div className="atom-container">
        <div className="atom-circle atom-oxygen">O</div>
        <div className="bond bond-single"></div>
        <div className="atom-circle atom-carbon">C</div>
        <div className="bond bond-single vertical-bond"></div>
        <div className="atom-circle atom-oxygen">O</div>
      </div>
    </div>
  );

  return (
    <div className="molecule-viewer-2d">
      <h3 className="viewer-title">2D Structure</h3>
      
      <div 
        className="molecule-viewer-container" 
        ref={containerRef}
        role="img"
        aria-label={`2D structure of ${name || 'molecule'}`}
      >
        {loading ? (
          <div className="molecule-viewer-loading">
            <div className="viewer-spinner" aria-hidden="true"></div>
            <p>{loadingStatus || 'Generating structure...'}</p>
          </div>
        ) : error ? (
          <div className="molecule-viewer-error">
            {/* Fallback placeholder for when we can't get a proper depiction */}
            <div className="molecule-fallback">
              <FallbackMolecule />
              <p className="error-message">{error}</p>
              <button 
                className="retry-button"
                onClick={handleRetry}
                aria-label="Retry loading structure"
              >
                Try Again
              </button>
            </div>
          </div>
        ) : (
          <div className="molecule-depiction">
            {svgContent ? (
              <div 
                dangerouslySetInnerHTML={{ __html: svgContent }} 
                className="svg-container"
                aria-hidden="true" // The parent div already has the accessible label
              />
            ) : (
              <div className="molecule-fallback">
                <FallbackMolecule />
                <p className="fallback-message">Using simplified structure representation</p>
              </div>
            )}
            <div className="molecule-details">
              <p className="molecule-name" id={`molecule-name-${name ? name.toLowerCase().replace(/\s+/g, '-') : 'unknown'}`}>
                {name || 'Molecule'}
              </p>
              {formula && (
                <p className="molecule-formula" title={`Molecular formula: ${formula}`}>
                  {formula}
                </p>
              )}
              
              {pubchem_cid && (
                <p className="molecule-pubchem-link">
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
                <div className="molecule-smiles-container">
                  <p 
                    className="molecule-smiles" 
                    title={smiles}
                  >
                    {smiles.length > 30 ? `${smiles.substring(0, 30)}...` : smiles}
                  </p>
                  <button 
                    className="copy-button"
                    onClick={() => {
                      navigator.clipboard.writeText(smiles)
                        .then(() => alert('SMILES copied to clipboard'))
                        .catch(err => console.error('Failed to copy SMILES:', err));
                    }}
                    aria-label="Copy SMILES notation to clipboard"
                    title="Copy SMILES"
                  >
                    Copy
                  </button>
                </div>
              )}
            </div>
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
          text-align: center;
          padding: 0 20px;
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
          margin: 10px 0 20px;
        }
        
        .atom-container {
          display: flex;
          align-items: center;
          position: relative;
        }
        
        .atom-circle {
          width: 30px;
          height: 30px;
          border-radius: 50%;
          display: flex;
          justify-content: center;
          align-items: center;
          font-weight: bold;
          color: white;
          font-size: 14px;
        }
        
        .atom-oxygen {
          background-color: #ff4d4d;
        }
        
        .atom-carbon {
          background-color: #333333;
        }
        
        .bond {
          height: 6px;
          background-color: #718096;
        }
        
        .bond-single {
          width: 35px;
        }
        
        .vertical-bond {
          position: absolute;
          transform: rotate(90deg);
          width: 35px;
          top: 45px;
          left: 12px;
        }
        
        .error-message {
          color: #e53e3e;
          font-size: 0.9rem;
          margin: 0 0 10px 0;
        }
        
        .fallback-message {
          color: #718096;
          font-size: 0.85rem;
          font-style: italic;
          margin: 10px 0;
        }
        
        .retry-button {
          background-color: #0070f3;
          color: white;
          padding: 8px 12px;
          border: none;
          border-radius: 4px;
          cursor: pointer;
          font-size: 0.9rem;
          transition: background-color 0.2s;
        }
        
        .retry-button:hover {
          background-color: #0051b3;
        }
        
        .molecule-depiction {
          display: flex;
          flex-direction: column;
          align-items: center;
          max-width: 100%;
        }
        
        .svg-container {
          max-width: 100%;
          overflow: hidden;
        }
        
        .svg-container svg {
          max-width: 100%;
          height: auto;
          display: block;
        }
        
        .molecule-details {
          display: flex;
          flex-direction: column;
          align-items: center;
          margin-top: 10px;
          width: 100%;
        }
        
        .molecule-name {
          font-weight: bold;
          margin-bottom: 5px;
          text-align: center;
          color: #2d3748;
        }
        
        .molecule-formula {
          font-family: monospace;
          color: #4a5568;
          margin: 0 0 5px 0;
          background-color: #edf2f7;
          padding: 3px 6px;
          border-radius: 4px;
          font-size: 0.9rem;
        }
        
        .molecule-pubchem-link {
          margin-top: 5px;
        }
        
        .pubchem-link {
          color: #0070f3;
          font-size: 0.9rem;
          text-decoration: none;
          display: inline-flex;
          align-items: center;
          transition: color 0.2s;
        }
        
        .pubchem-link:hover {
          color: #0051b3;
          text-decoration: underline;
        }
        
        .molecule-smiles-container {
          display: flex;
          align-items: center;
          margin-top: 8px;
          max-width: 100%;
          flex-wrap: wrap;
          justify-content: center;
          gap: 8px;
        }
        
        .molecule-smiles {
          font-family: monospace;
          background-color: #edf2f7;
          padding: 5px 10px;
          border-radius: 4px;
          font-size: 0.9rem;
          margin: 0;
          white-space: nowrap;
          overflow: hidden;
          text-overflow: ellipsis;
          max-width: 250px;
        }
        
        .copy-button {
          background-color: #0070f3;
          color: white;
          border: none;
          border-radius: 4px;
          padding: 5px 8px;
          font-size: 0.8rem;
          cursor: pointer;
          transition: background-color 0.2s;
        }
        
        .copy-button:hover {
          background-color: #0051b3;
        }
        
        @media (max-width: 768px) {
          .molecule-viewer-container {
            padding: 10px;
          }
          
          .molecule-formula {
            font-size: 0.8rem;
          }
          
          .molecule-smiles {
            font-size: 0.8rem;
            max-width: 200px;
          }
        }
      `}</style>
    </div>
  );
}
import React, { useState, useEffect } from 'react';
import Layout from '../components/Layout';
import { getMoleculeDepiction } from '../utils/api';
import Loading from '../components/Loading';
import ErrorMessage from '../components/ErrorMessage';

export default function MoleculeViewer() {
  const [smiles, setSmiles] = useState('c1ccccc1CC');
  const [depiction, setDepiction] = useState('');
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);

  // Sample SMILES for users to try
  const sampleMolecules = [
    { name: 'Glycerol', smiles: 'C(C(CO)O)O' },
    { name: 'DMSO', smiles: 'CS(=O)C' },
    { name: 'Ethylene Glycol', smiles: 'C(CO)O' },
    { name: 'Propylene Glycol', smiles: 'CC(CO)O' },
    { name: 'Benzene', smiles: 'c1ccccc1' },
    { name: 'Caffeine', smiles: 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C' }
  ];

  const fetchDepiction = async () => {
    if (!smiles) return;

    try {
      setLoading(true);
      setError(null);
      const svgData = await getMoleculeDepiction(smiles);
      setDepiction(svgData);
    } catch (err) {
      console.error('Failed to get molecule depiction:', err);
      setError('Failed to generate molecule visualization. Please try again.');
    } finally {
      setLoading(false);
    }
  };

  // Generate depiction when SMILES changes
  useEffect(() => {
    fetchDepiction();
  }, [smiles]);

  return (
    <Layout title="Molecule Viewer">
      <div className="page-header">
        <h1 className="title">Molecule Viewer</h1>
        <p className="description">
          Visualize molecular structures using the RDKit service on Fly.io
        </p>
      </div>

      <div className="card molecule-viewer-card">
        <div className="molecule-input">
          <label htmlFor="smiles-input">SMILES Notation:</label>
          <input
            id="smiles-input"
            type="text"
            value={smiles}
            onChange={(e) => setSmiles(e.target.value)}
            placeholder="Enter SMILES notation"
            className="smiles-input"
          />
          <button onClick={fetchDepiction} className="button">
            View Molecule
          </button>
        </div>

        <div className="sample-molecules">
          <h3>Sample Molecules:</h3>
          <div className="sample-buttons">
            {sampleMolecules.map((molecule) => (
              <button
                key={molecule.name}
                onClick={() => setSmiles(molecule.smiles)}
                className="sample-button"
              >
                {molecule.name}
              </button>
            ))}
          </div>
        </div>

        <div className="molecule-depiction">
          {loading ? (
            <Loading />
          ) : error ? (
            <ErrorMessage message={error} onRetry={fetchDepiction} />
          ) : depiction ? (
            <div dangerouslySetInnerHTML={{ __html: depiction }} />
          ) : (
            <div className="molecule-placeholder">
              Enter a SMILES notation and click "View Molecule"
            </div>
          )}
        </div>
      </div>

      <style jsx>{`
        .molecule-viewer-card {
          max-width: 800px;
        }
        
        .molecule-input {
          display: flex;
          flex-direction: column;
          gap: 10px;
          margin-bottom: 20px;
        }
        
        .smiles-input {
          padding: 10px;
          border: 1px solid #ddd;
          border-radius: 4px;
          font-family: monospace;
          width: 100%;
        }
        
        .sample-molecules {
          margin-bottom: 20px;
        }
        
        .sample-buttons {
          display: flex;
          flex-wrap: wrap;
          gap: 10px;
          margin-top: 10px;
        }
        
        .sample-button {
          padding: 8px 12px;
          background-color: #f0f0f0;
          border: 1px solid #ddd;
          border-radius: 4px;
          cursor: pointer;
          transition: background-color 0.2s;
        }
        
        .sample-button:hover {
          background-color: #e0e0e0;
        }
        
        .molecule-depiction {
          border: 1px solid #ddd;
          border-radius: 8px;
          padding: 20px;
          min-height: 300px;
          display: flex;
          justify-content: center;
          align-items: center;
          background-color: white;
        }
        
        @media (max-width: 600px) {
          .sample-buttons {
            flex-direction: column;
          }
        }
      `}</style>
    </Layout>
  );
}
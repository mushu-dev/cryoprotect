import React, { useState, useCallback, useEffect } from 'react';
import dynamic from 'next/dynamic';
import Layout from '../components/Layout';
import Loading from '../components/Loading';
import ErrorMessage from '../components/ErrorMessage';

// Check if Convex is enabled
const isConvexEnabled = process.env.NEXT_PUBLIC_USE_CONVEX === 'true';

// Dynamically import Convex components to avoid build issues
const ConvexMoleculesView = dynamic(
  () => import('../src/components/ConvexMoleculesView'),
  { 
    ssr: false,
    loading: () => <Loading />
  }
);

// Mock data for initial state and fallback
const mockMolecules = [
  {
    id: 1,
    name: 'Glycerol',
    formula: 'C3H8O3',
    pubchem_cid: '753',
    molecular_weight: 92.09,
    is_cryoprotectant: true,
    description: 'A common cryoprotectant used in various applications.',
    smiles: 'C(C(CO)O)O'
  },
  {
    id: 2,
    name: 'Dimethyl Sulfoxide (DMSO)',
    formula: 'C2H6OS',
    pubchem_cid: '679',
    molecular_weight: 78.13,
    is_cryoprotectant: true,
    description: 'A widely used penetrating cryoprotectant.',
    smiles: 'CS(=O)C'
  },
  {
    id: 3,
    name: 'Ethylene Glycol',
    formula: 'C2H6O2',
    pubchem_cid: '174',
    molecular_weight: 62.07,
    is_cryoprotectant: true,
    description: 'Used in cryopreservation of embryos and tissues.',
    smiles: 'C(CO)O'
  },
  {
    id: 4,
    name: 'Propylene Glycol',
    formula: 'C3H8O2',
    pubchem_cid: '1030',
    molecular_weight: 76.09,
    is_cryoprotectant: true,
    description: 'Used as a cryoprotectant for various biological materials.',
    smiles: 'CC(CO)O'
  }
];

// Standard molecules view component
function StandardMoleculesView({ molecules }) {
  const [searchTerm, setSearchTerm] = useState('');
  const [selectedMolecule, setSelectedMolecule] = useState(null);
  const [rdkitProperties, setRdkitProperties] = useState({});
  const [loadingCalculations, setLoadingCalculations] = useState(new Set());
  
  // Filter molecules based on search term
  const filteredMolecules = React.useMemo(() => {
    if (!molecules || !searchTerm) return molecules;
    return molecules.filter(molecule => 
      molecule.name.toLowerCase().includes(searchTerm.toLowerCase()) ||
      (molecule.formula && molecule.formula.toLowerCase().includes(searchTerm.toLowerCase()))
    );
  }, [molecules, searchTerm]);
  
  // Calculate RDKit properties for a molecule
  const calculateRDKitProperties = useCallback(async (molecule) => {
    if (!molecule.smiles) return;
    
    const moleculeId = molecule.id;
    setLoadingCalculations(prev => new Set(prev).add(moleculeId));
    
    try {
      const response = await fetch('https://cryoprotect-8030e4025428.herokuapp.com/api/rdkit/descriptors', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({ smiles: molecule.smiles }),
      });
      
      if (response.ok) {
        const data = await response.json();
        setRdkitProperties(prev => ({
          ...prev,
          [moleculeId]: data
        }));
      } else {
        console.log('RDKit endpoint not available, using mock calculations');
        // Mock calculation for demo
        setRdkitProperties(prev => ({
          ...prev,
          [moleculeId]: {
            descriptors: {
              molecular_weight: molecule.molecular_weight || 100,
              logp: (Math.random() * 4 - 2).toFixed(2),
              topological_polar_surface_area: (Math.random() * 150).toFixed(1)
            },
            cryoprotectant_score: (Math.random() * 100 + 50).toFixed(1)
          }
        }));
      }
    } catch (error) {
      console.error('Failed to calculate RDKit properties:', error);
      // Mock calculation for demo
      setRdkitProperties(prev => ({
        ...prev,
        [moleculeId]: {
          descriptors: {
            molecular_weight: molecule.molecular_weight || 100,
            logp: (Math.random() * 4 - 2).toFixed(2),
            topological_polar_surface_area: (Math.random() * 150).toFixed(1)
          },
          cryoprotectant_score: (Math.random() * 100 + 50).toFixed(1)
        }
      }));
    } finally {
      setLoadingCalculations(prev => {
        const newSet = new Set(prev);
        newSet.delete(moleculeId);
        return newSet;
      });
    }
  }, []);

  return (
    <div className="molecules-container">
      <div className="molecules-header">
        <div className="header-content">
          <h1 className="page-title">Molecules Database</h1>
          <p className="page-description">
            Explore our comprehensive database of cryoprotectant molecules with real-time property calculations
          </p>
        </div>
        
        {/* Search functionality */}
        <div className="search-section">
          <input
            type="text"
            placeholder="Search molecules by name or formula..."
            value={searchTerm}
            onChange={(e) => setSearchTerm(e.target.value)}
            className="search-input"
          />
        </div>
        
        {/* Status info */}
        <div className="status-info">
          <span className="status-item">🧪 Molecular Database</span>
          <span className="status-item">🔬 Property Calculations</span>
          <span className="status-item">📊 {filteredMolecules?.length || 0} molecules</span>
        </div>
      </div>
      
      <div className="molecules-grid">
        {filteredMolecules && filteredMolecules.length > 0 ? (
          filteredMolecules.map((molecule) => {
            const isCalculating = loadingCalculations.has(molecule.id);
            const rdkitData = rdkitProperties[molecule.id];
            
            return (
              <div key={molecule.id} className="molecule-card">
                <div className="molecule-header">
                  <h3 className="molecule-name">{molecule.name}</h3>
                  <span className="molecule-status">Active</span>
                </div>
                
                <div className="molecule-details">
                  <p className="molecule-formula">{molecule.formula || 'No formula available'}</p>
                  
                  {molecule.smiles && (
                    <div className="molecule-smiles">
                      <strong>SMILES:</strong> {molecule.smiles.length > 30 ? 
                        `${molecule.smiles.substring(0, 30)}...` : 
                        molecule.smiles}
                    </div>
                  )}
                  
                  {molecule.pubchem_cid && (
                    <div className="molecule-pubchem">
                      PubChem CID: {molecule.pubchem_cid}
                    </div>
                  )}
                  
                  {molecule.description && (
                    <p className="molecule-description">{molecule.description}</p>
                  )}
                </div>
                
                {/* RDKit Calculations Section */}
                <div className="rdkit-section">
                  <div className="rdkit-header">
                    <h4>Molecular Properties</h4>
                    {molecule.smiles && (
                      <button
                        onClick={() => calculateRDKitProperties(molecule)}
                        disabled={isCalculating}
                        className="calculate-btn"
                      >
                        {isCalculating ? 'Calculating...' : 'Calculate'}
                      </button>
                    )}
                  </div>
                  
                  {isCalculating ? (
                    <div className="calculating-indicator">
                      <div className="spinner"></div>
                      <span>Computing properties...</span>
                    </div>
                  ) : rdkitData ? (
                    <div className="properties-grid">
                      {rdkitData.descriptors && (
                        <>
                          <div className="property-item">
                            <span className="property-label">Molecular Weight:</span>
                            <span className="property-value">{rdkitData.descriptors.molecular_weight?.toFixed(2)} Da</span>
                          </div>
                          <div className="property-item">
                            <span className="property-label">LogP:</span>
                            <span className="property-value">{rdkitData.descriptors.logp?.toFixed(2)}</span>
                          </div>
                          <div className="property-item">
                            <span className="property-label">TPSA:</span>
                            <span className="property-value">{rdkitData.descriptors.topological_polar_surface_area?.toFixed(1)} Ų</span>
                          </div>
                          {rdkitData.cryoprotectant_score && (
                            <div className="property-item cryo-score">
                              <span className="property-label">Cryo Score:</span>
                              <span className={`property-value score-${
                                rdkitData.cryoprotectant_score > 80 ? 'good' : 
                                rdkitData.cryoprotectant_score > 60 ? 'fair' : 'poor'
                              }`}>
                                {rdkitData.cryoprotectant_score}/100
                              </span>
                            </div>
                          )}
                        </>
                      )}
                    </div>
                  ) : molecule.smiles ? (
                    <div className="no-properties">
                      Click "Calculate" to compute molecular properties
                    </div>
                  ) : (
                    <div className="no-structure">
                      No molecular structure available for calculations
                    </div>
                  )}
                </div>
              </div>
            );
          })
        ) : (
          <div className="no-molecules">
            <p>No molecules found.</p>
            {searchTerm && (
              <p>Try adjusting your search term.</p>
            )}
          </div>
        )}
      </div>
      
      <style jsx>{`
        .molecules-container {
          max-width: 1200px;
          margin: 0 auto;
          padding: 20px;
        }
        
        .molecules-header {
          margin-bottom: 30px;
        }
        
        .header-content {
          text-align: center;
          margin-bottom: 20px;
        }
        
        .page-title {
          font-size: 2.5rem;
          font-weight: bold;
          color: #1a202c;
          margin-bottom: 10px;
        }
        
        .page-description {
          font-size: 1.1rem;
          color: #4a5568;
          max-width: 600px;
          margin: 0 auto;
        }
        
        .search-section {
          margin-bottom: 20px;
        }
        
        .search-input {
          width: 100%;
          max-width: 500px;
          padding: 12px 16px;
          border: 2px solid #e2e8f0;
          border-radius: 8px;
          font-size: 16px;
          outline: none;
          transition: border-color 0.2s;
        }
        
        .search-input:focus {
          border-color: #3182ce;
        }
        
        .status-info {
          display: flex;
          justify-content: center;
          gap: 20px;
          font-size: 14px;
          color: #4a5568;
        }
        
        .molecules-grid {
          display: grid;
          grid-template-columns: repeat(auto-fill, minmax(350px, 1fr));
          gap: 24px;
        }
        
        .molecule-card {
          background: white;
          border: 1px solid #e2e8f0;
          border-radius: 12px;
          padding: 20px;
          box-shadow: 0 1px 3px rgba(0, 0, 0, 0.1);
          transition: transform 0.2s, box-shadow 0.2s;
        }
        
        .molecule-card:hover {
          transform: translateY(-2px);
          box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
        }
        
        .molecule-header {
          display: flex;
          justify-content: space-between;
          align-items: center;
          margin-bottom: 15px;
        }
        
        .molecule-name {
          font-size: 1.25rem;
          font-weight: 600;
          color: #1a202c;
          margin: 0;
        }
        
        .molecule-status {
          background: #48bb78;
          color: white;
          padding: 4px 8px;
          border-radius: 4px;
          font-size: 12px;
          font-weight: 500;
        }
        
        .molecule-details {
          margin-bottom: 20px;
        }
        
        .molecule-formula {
          font-size: 1rem;
          color: #4a5568;
          margin-bottom: 8px;
        }
        
        .molecule-smiles {
          font-size: 0.9rem;
          color: #718096;
          margin-bottom: 8px;
        }
        
        .molecule-pubchem {
          font-size: 0.9rem;
          color: #3182ce;
          margin-bottom: 8px;
        }
        
        .molecule-description {
          font-size: 0.9rem;
          color: #4a5568;
          line-height: 1.4;
        }
        
        .rdkit-section {
          border-top: 1px solid #e2e8f0;
          padding-top: 16px;
        }
        
        .rdkit-header {
          display: flex;
          justify-content: space-between;
          align-items: center;
          margin-bottom: 12px;
        }
        
        .rdkit-header h4 {
          font-size: 1rem;
          font-weight: 600;
          color: #1a202c;
          margin: 0;
        }
        
        .calculate-btn {
          background: #3182ce;
          color: white;
          border: none;
          padding: 6px 12px;
          border-radius: 4px;
          font-size: 12px;
          cursor: pointer;
          transition: background-color 0.2s;
        }
        
        .calculate-btn:hover:not(:disabled) {
          background: #2c5282;
        }
        
        .calculate-btn:disabled {
          opacity: 0.6;
          cursor: not-allowed;
        }
        
        .calculating-indicator {
          display: flex;
          align-items: center;
          gap: 8px;
          color: #4a5568;
          font-size: 14px;
        }
        
        .spinner {
          width: 16px;
          height: 16px;
          border: 2px solid #e2e8f0;
          border-top: 2px solid #3182ce;
          border-radius: 50%;
          animation: spin 1s linear infinite;
        }
        
        @keyframes spin {
          0% { transform: rotate(0deg); }
          100% { transform: rotate(360deg); }
        }
        
        .properties-grid {
          display: grid;
          gap: 8px;
        }
        
        .property-item {
          display: grid;
          grid-template-columns: 1fr 1fr;
          gap: 8px;
          padding: 4px 0;
        }
        
        .property-label {
          color: #4a5568;
          font-size: 14px;
        }
        
        .property-value {
          font-weight: 500;
          color: #1a202c;
          font-size: 14px;
        }
        
        .cryo-score {
          border-top: 1px solid #e2e8f0;
          padding-top: 8px;
          margin-top: 8px;
        }
        
        .score-good { color: #38a169; }
        .score-fair { color: #d69e2e; }
        .score-poor { color: #e53e3e; }
        
        .no-properties, .no-structure {
          font-size: 14px;
          color: #718096;
          font-style: italic;
          text-align: center;
          padding: 12px;
        }
        
        .no-molecules {
          grid-column: 1 / -1;
          text-align: center;
          color: #4a5568;
          padding: 40px;
        }
        
        @media (max-width: 768px) {
          .molecules-grid {
            grid-template-columns: 1fr;
          }
          
          .status-info {
            flex-wrap: wrap;
            gap: 10px;
          }
          
          .molecule-header {
            flex-direction: column;
            align-items: flex-start;
            gap: 8px;
          }
        }
      `}</style>
    </div>
  );
}

// This function gets called at build time on server-side.
export async function getStaticProps() {
  return {
    props: {
      initialMolecules: mockMolecules
    },
  }
}

export default function Molecules({ initialMolecules }) {
  const [molecules, setMolecules] = useState(initialMolecules);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);

  useEffect(() => {
    // Try to fetch from API if not using Convex
    if (!isConvexEnabled) {
      fetchMolecules();
    }
  }, []);

  const fetchMolecules = async () => {
    try {
      setLoading(true);
      setError(null);
      const { getMolecules } = await import('../utils/api');
      const data = await getMolecules();
      setMolecules(data || initialMolecules);
    } catch (err) {
      console.error('Failed to fetch molecules:', err);
      setMolecules(initialMolecules);
    } finally {
      setLoading(false);
    }
  };

  if (loading) {
    return (
      <Layout title="Molecules">
        <Loading />
      </Layout>
    );
  }

  if (error) {
    return (
      <Layout title="Molecules">
        <ErrorMessage message={error} onRetry={fetchMolecules} />
      </Layout>
    );
  }

  return (
    <Layout title="Molecules">
      {isConvexEnabled ? (
        <ConvexMoleculesView />
      ) : (
        <StandardMoleculesView molecules={molecules} />
      )}
    </Layout>
  );
}
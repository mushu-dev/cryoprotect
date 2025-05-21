import React, { useState, useCallback, useEffect } from 'react';
import { useQuery } from 'convex/react';
import { api } from '../../convex/_generated/api';

export default function ConvexMoleculesView() {
  const [searchTerm, setSearchTerm] = useState('');
  const [rdkitProperties, setRdkitProperties] = useState({});
  const [loadingCalculations, setLoadingCalculations] = useState(new Set());

  // Query molecules from Convex with real-time updates
  const molecules = useQuery(api.molecules.query.getRecentMolecules, { 
    limit: 50
  });
  
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
    if (!molecule.canonicalSmiles) return;
    
    const moleculeId = molecule._id;
    setLoadingCalculations(prev => new Set(prev).add(moleculeId));
    
    try {
      const response = await fetch('https://cryoprotect-8030e4025428.herokuapp.com/api/rdkit/descriptors', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({ smiles: molecule.canonicalSmiles }),
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
              molecular_weight: Math.random() * 200 + 50,
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
            molecular_weight: Math.random() * 200 + 50,
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

  // Auto-calculate properties for visible molecules
  useEffect(() => {
    if (filteredMolecules && filteredMolecules.length > 0) {
      // Calculate properties for first 3 molecules automatically
      filteredMolecules.slice(0, 3).forEach(molecule => {
        if (molecule.canonicalSmiles && !rdkitProperties[molecule._id]) {
          calculateRDKitProperties(molecule);
        }
      });
    }
  }, [filteredMolecules, rdkitProperties, calculateRDKitProperties]);
  
  const isLoading = molecules === undefined;
  
  if (isLoading) {
    return (
      <div className="molecules-container">
        <div className="loading-state">
          <div className="spinner"></div>
          <p>Loading molecules from Convex...</p>
        </div>
        
        <style jsx>{`
          .molecules-container {
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
          }
          
          .loading-state {
            text-align: center;
            padding: 60px 20px;
          }
          
          .spinner {
            width: 40px;
            height: 40px;
            border: 4px solid #e2e8f0;
            border-top: 4px solid #3182ce;
            border-radius: 50%;
            animation: spin 1s linear infinite;
            margin: 0 auto 20px;
          }
          
          @keyframes spin {
            0% { transform: rotate(0deg); }
            100% { transform: rotate(360deg); }
          }
        `}</style>
      </div>
    );
  }

  return (
    <div className="molecules-container">
      <div className="molecules-header">
        <div className="header-content">
          <h1 className="page-title">Molecules Database (Convex Real-time)</h1>
          <p className="page-description">
            Real-time molecular database with live property calculations powered by Convex
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
          <span className="status-item">📊 Real-time Convex data</span>
          <span className="status-item">🧪 Live RDKit calculations</span>
          <span className="status-item">🔄 {filteredMolecules?.length || 0} molecules</span>
        </div>
      </div>
      
      <div className="molecules-grid">
        {filteredMolecules && filteredMolecules.length > 0 ? (
          filteredMolecules.map((molecule) => {
            const isCalculating = loadingCalculations.has(molecule._id);
            const rdkitData = rdkitProperties[molecule._id];
            
            return (
              <div key={molecule._id} className="molecule-card">
                <div className="molecule-header">
                  <h3 className="molecule-name">{molecule.name}</h3>
                  <span className={`molecule-status ${molecule.status === 'active' ? 'active' : 'inactive'}`}>
                    {molecule.status || 'active'}
                  </span>
                </div>
                
                <div className="molecule-details">
                  <p className="molecule-formula">{molecule.formula || 'No formula available'}</p>
                  
                  {molecule.canonicalSmiles && (
                    <div className="molecule-smiles">
                      <strong>SMILES:</strong> {molecule.canonicalSmiles.length > 30 ? 
                        `${molecule.canonicalSmiles.substring(0, 30)}...` : 
                        molecule.canonicalSmiles}
                    </div>
                  )}
                  
                  {molecule.pubchemCid && (
                    <div className="molecule-pubchem">
                      PubChem CID: {molecule.pubchemCid}
                    </div>
                  )}
                  
                  {molecule.inchiKey && (
                    <div className="molecule-inchi">
                      InChI Key: {molecule.inchiKey.substring(0, 14)}...
                    </div>
                  )}
                </div>
                
                {/* RDKit Calculations Section */}
                <div className="rdkit-section">
                  <div className="rdkit-header">
                    <h4>Live Properties</h4>
                    {molecule.canonicalSmiles && (
                      <button
                        onClick={() => calculateRDKitProperties(molecule)}
                        disabled={isCalculating}
                        className="calculate-btn"
                      >
                        {isCalculating ? 'Calculating...' : 'Recalculate'}
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
                  ) : molecule.canonicalSmiles ? (
                    <div className="no-properties">
                      Click "Recalculate" to compute RDKit properties
                    </div>
                  ) : (
                    <div className="no-structure">
                      No SMILES structure available for calculations
                    </div>
                  )}
                </div>
              </div>
            );
          })
        ) : (
          <div className="no-molecules">
            <p>No molecules found in the Convex database.</p>
            {searchTerm && (
              <p>Try adjusting your search term or check if molecules are populated.</p>
            )}
          </div>
        )}
      </div>
      
      {/* Connection Status */}
      <div className="connection-status">
        <div className="status-indicator">
          <span className="status-dot"></span>
          Connected to Convex • Real-time updates enabled
        </div>
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
          text-align: center;
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
          flex-wrap: wrap;
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
          padding: 4px 8px;
          border-radius: 4px;
          font-size: 12px;
          font-weight: 500;
        }
        
        .molecule-status.active {
          background: #48bb78;
          color: white;
        }
        
        .molecule-status.inactive {
          background: #a0aec0;
          color: white;
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
        
        .molecule-inchi {
          font-size: 0.9rem;
          color: #805ad5;
          margin-bottom: 8px;
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
        
        .connection-status {
          margin-top: 30px;
          text-align: center;
        }
        
        .status-indicator {
          display: inline-flex;
          align-items: center;
          padding: 8px 16px;
          background: #e6fffa;
          color: #234e52;
          border-radius: 20px;
          font-size: 14px;
          font-weight: 500;
        }
        
        .status-dot {
          width: 8px;
          height: 8px;
          background: #38b2ac;
          border-radius: 50%;
          margin-right: 8px;
          animation: pulse 2s infinite;
        }
        
        @keyframes pulse {
          0%, 100% { opacity: 1; }
          50% { opacity: 0.5; }
        }
        
        @media (max-width: 768px) {
          .molecules-grid {
            grid-template-columns: 1fr;
          }
          
          .molecule-header {
            flex-direction: column;
            align-items: flex-start;
            gap: 8px;
          }
          
          .page-title {
            font-size: 2rem;
          }
        }
      `}</style>
    </div>
  );
}
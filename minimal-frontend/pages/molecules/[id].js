import React, { useState, useEffect } from 'react';
import { useRouter } from 'next/router';
import Link from 'next/link';
import dynamic from 'next/dynamic';
import Layout from '../../components/Layout';
import Loading from '../../components/Loading';
import ErrorMessage from '../../components/ErrorMessage';
import { getMolecule } from '../../utils/api';

// Dynamically import components with SSR disabled for visualization
const MoleculeViewer3D = dynamic(() => import('../../components/MoleculeViewer3D'), {
  ssr: false,
  loading: () => <div className="dynamic-loading">Loading 3D viewer...</div>
});

const MoleculeViewer2D = dynamic(() => import('../../components/MoleculeViewer2D'), {
  ssr: false,
  loading: () => <div className="dynamic-loading">Loading 2D viewer...</div>
});

const PropertyDisplay = dynamic(() => import('../../components/PropertyDisplay'), {
  loading: () => <div className="dynamic-loading">Loading properties...</div>
});

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
    smiles: 'C(C(CO)O)O',
    properties: [
      { name: 'Melting Point', value: '18.2', units: '°C' },
      { name: 'Boiling Point', value: '290', units: '°C' },
      { name: 'Density', value: '1.26', units: 'g/mL' }
    ]
  },
  {
    id: 2,
    name: 'Dimethyl Sulfoxide (DMSO)',
    formula: 'C2H6OS',
    pubchem_cid: '679',
    molecular_weight: 78.13,
    is_cryoprotectant: true,
    description: 'A widely used penetrating cryoprotectant.',
    smiles: 'CS(=O)C',
    properties: [
      { name: 'Melting Point', value: '18.5', units: '°C' },
      { name: 'Boiling Point', value: '189', units: '°C' },
      { name: 'Density', value: '1.1', units: 'g/mL' }
    ]
  },
  {
    id: 3,
    name: 'Ethylene Glycol',
    formula: 'C2H6O2',
    pubchem_cid: '174',
    molecular_weight: 62.07,
    is_cryoprotectant: true,
    description: 'Used in cryopreservation of embryos and tissues.',
    smiles: 'C(CO)O',
    properties: [
      { name: 'Melting Point', value: '-12.9', units: '°C' },
      { name: 'Boiling Point', value: '197.3', units: '°C' },
      { name: 'Density', value: '1.11', units: 'g/mL' }
    ]
  },
  {
    id: 4,
    name: 'Propylene Glycol',
    formula: 'C3H8O2',
    pubchem_cid: '1030',
    molecular_weight: 76.09,
    is_cryoprotectant: true,
    description: 'Used as a cryoprotectant for various biological materials.',
    smiles: 'CC(CO)O',
    properties: [
      { name: 'Melting Point', value: '-59', units: '°C' },
      { name: 'Boiling Point', value: '188.2', units: '°C' },
      { name: 'Density', value: '1.04', units: 'g/mL' }
    ]
  }
];

export async function getStaticPaths() {
  return {
    paths: mockMolecules.map(molecule => ({
      params: { id: molecule.id.toString() }
    })),
    fallback: false // Only pre-render the molecules we know about
  };
}

export async function getStaticProps({ params }) {
  const { id } = params;
  const molecule = mockMolecules.find(m => m.id.toString() === id);
  
  return {
    props: {
      initialMolecule: molecule || null,
      id
    }
  };
}

export default function MoleculeDetail({ initialMolecule, id }) {
  const router = useRouter();
  
  const [molecule, setMolecule] = useState(initialMolecule);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  
  useEffect(() => {
    // If we have an initialMolecule from props, use it
    if (initialMolecule) {
      setMolecule(initialMolecule);
      setLoading(false);
    } else if (id) {
      // Only try to fetch from API if we don't have initial data
      fetchMolecule(id);
    }
  }, []);
  
  const fetchMolecule = async (moleculeId) => {
    try {
      setLoading(true);
      setError(null);
      const data = await getMolecule(moleculeId);
      setMolecule(data);
    } catch (err) {
      console.error('Failed to fetch molecule:', err);
      // Try to get mock data as fallback
      const mockMolecule = mockMolecules.find(m => m.id.toString() === moleculeId.toString());
      if (mockMolecule) {
        setMolecule(mockMolecule);
      } else {
        setError('Failed to fetch molecule details. Please try again later.');
      }
    } finally {
      setLoading(false);
    }
  };
  
  // This is unlikely to happen with getStaticProps/getStaticPaths, but keeping for safety
  if (!id && !initialMolecule) {
    return (
      <Layout title="Molecule Detail">
        <div className="page-header">
          <h1 className="title">Molecule Detail</h1>
        </div>
        <div className="card">
          <p>Loading molecule ID...</p>
        </div>
      </Layout>
    );
  }
  
  return (
    <Layout title={molecule ? `${molecule.name || 'Molecule'} Details` : 'Molecule Detail'}>
      <div className="page-header">
        <Link href="/molecules">
          <span className="back-link">← Back to Molecules</span>
        </Link>
        <h1 className="title">{molecule ? molecule.name || 'Unnamed Molecule' : 'Molecule Detail'}</h1>
      </div>
      
      {loading ? (
        <Loading />
      ) : error ? (
        <ErrorMessage message={error} onRetry={() => fetchMolecule(id)} />
      ) : molecule ? (
        <>
          <div className="molecule-detail-container">
            <div className="card molecule-detail-card">
              <h2 className="subtitle">Basic Information</h2>
              
              <div className="property-group">
                {molecule.pubchem_cid && (
                  <div className="property-row">
                    <span className="property-label">PubChem CID:</span>
                    <span className="property-value">{molecule.pubchem_cid}</span>
                  </div>
                )}
                
                {molecule.formula && (
                  <div className="property-row">
                    <span className="property-label">Formula:</span>
                    <span className="property-value molecule-formula">{molecule.formula}</span>
                  </div>
                )}
                
                {molecule.molecular_weight && (
                  <div className="property-row">
                    <span className="property-label">Molecular Weight:</span>
                    <span className="property-value">{molecule.molecular_weight} g/mol</span>
                  </div>
                )}
                
                {typeof molecule.is_cryoprotectant !== 'undefined' && (
                  <div className="property-row">
                    <span className="property-label">Cryoprotectant:</span>
                    <span className="property-value">{molecule.is_cryoprotectant ? 'Yes' : 'No'}</span>
                  </div>
                )}
                
                {molecule.smiles && (
                  <div className="property-row">
                    <span className="property-label">SMILES:</span>
                    <span className="property-value molecule-smiles">{molecule.smiles}</span>
                  </div>
                )}
              </div>
            </div>
            
            {molecule.description && (
              <div className="card molecule-detail-card">
                <h2 className="subtitle">Description</h2>
                <p>{molecule.description}</p>
              </div>
            )}
            
            <div className="card molecule-detail-card">
              <h2 className="subtitle">Molecular Visualization</h2>
              <div className="visualization-container">
                <div className="visualization-section">
                  {molecule.smiles ? (
                    <div className="visualization-3d">
                      <MoleculeViewer3D 
                        smiles={molecule.smiles} 
                        name={molecule.name} 
                      />
                    </div>
                  ) : (
                    <div className="molecule-placeholder">
                      SMILES notation required for 3D visualization
                    </div>
                  )}
                </div>
                
                <div className="visualization-section">
                  {molecule.smiles ? (
                    <div className="visualization-2d">
                      <MoleculeViewer2D 
                        smiles={molecule.smiles} 
                        name={molecule.name} 
                      />
                    </div>
                  ) : (
                    <div className="molecule-placeholder">
                      SMILES notation required for 2D visualization
                    </div>
                  )}
                </div>
              </div>
            </div>
            
            {molecule.properties && molecule.properties.length > 0 && (
              <div className="card molecule-detail-card">
                <h2 className="subtitle">Physical Properties</h2>
                <PropertyDisplay properties={molecule.properties} />
              </div>
            )}
          </div>
          
          <div className="actions-container">
            <Link href="/molecules">
              <span className="button">Back to Molecules</span>
            </Link>
          </div>
        </>
      ) : (
        <div className="card">
          <p>No molecule found with ID: {id}</p>
          <Link href="/molecules">
            <span className="button">Back to Molecules</span>
          </Link>
        </div>
      )}
      
      <style jsx>{`
        .molecule-detail-container {
          max-width: 900px;
        }
        
        .visualization-container {
          display: grid;
          grid-template-columns: 1fr 1fr;
          gap: 20px;
          margin-top: 15px;
        }
        
        .visualization-section {
          width: 100%;
        }
        
        .dynamic-loading {
          height: 200px;
          display: flex;
          justify-content: center;
          align-items: center;
          background-color: #f8f9fa;
          border-radius: 8px;
          color: #4a5568;
        }
        
        @media (max-width: 768px) {
          .visualization-container {
            grid-template-columns: 1fr;
          }
        }
      `}</style>
    </Layout>
  );
}
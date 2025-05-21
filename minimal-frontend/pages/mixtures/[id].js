import React, { useState, useEffect } from 'react';
import { useRouter } from 'next/router';
import Link from 'next/link';
import Layout from '../../components/Layout';
import Loading from '../../components/Loading';
import ErrorMessage from '../../components/ErrorMessage';
import { getMixture } from '../../utils/api';

// Mock data for initial state and fallback
const mockMixtures = [
  {
    id: 1,
    name: 'VS55',
    description: 'A cryoprotectant mixture used for organ preservation.',
    freezing_point: -80,
    created_at: '2024-05-01T12:00:00Z',
    components: [
      {
        molecule: {
          id: 1,
          name: 'Glycerol',
          formula: 'C3H8O3'
        },
        concentration: 15,
        concentration_unit: '%',
        role: 'Primary cryoprotectant'
      },
      {
        molecule: {
          id: 2,
          name: 'Dimethyl Sulfoxide (DMSO)',
          formula: 'C2H6OS'
        },
        concentration: 20,
        concentration_unit: '%',
        role: 'Penetrating cryoprotectant'
      },
      {
        molecule: {
          id: 3,
          name: 'Ethylene Glycol',
          formula: 'C2H6O2'
        },
        concentration: 20,
        concentration_unit: '%',
        role: 'Penetrating cryoprotectant'
      }
    ],
    properties: [
      { name: 'Glass Transition Temperature', value: '-123', units: '°C' },
      { name: 'Toxicity', value: 'Low', units: '' },
      { name: 'Viscosity', value: '4.2', units: 'cP' }
    ]
  },
  {
    id: 2,
    name: 'Standard Cell Freezing Medium',
    description: 'Common mixture for cell line preservation.',
    freezing_point: -20,
    created_at: '2024-05-02T10:30:00Z',
    components: [
      {
        molecule: {
          id: 2,
          name: 'Dimethyl Sulfoxide (DMSO)',
          formula: 'C2H6OS'
        },
        concentration: 10,
        concentration_unit: '%',
        role: 'Primary cryoprotectant'
      }
    ],
    properties: [
      { name: 'Glass Transition Temperature', value: '-132', units: '°C' },
      { name: 'Toxicity', value: 'Medium', units: '' }
    ]
  }
];

export async function getStaticPaths() {
  return {
    paths: mockMixtures.map(mixture => ({
      params: { id: mixture.id.toString() }
    })),
    fallback: false // Only pre-render the mixtures we know about
  };
}

export async function getStaticProps({ params }) {
  const { id } = params;
  const mixture = mockMixtures.find(m => m.id.toString() === id);
  
  return {
    props: {
      initialMixture: mixture || null,
      id
    }
  };
}

export default function MixtureDetail({ initialMixture, id }) {
  const router = useRouter();
  
  const [mixture, setMixture] = useState(initialMixture);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  
  useEffect(() => {
    // If we have an initialMixture from props, use it
    if (initialMixture) {
      setMixture(initialMixture);
      setLoading(false);
    } else if (id) {
      // Only try to fetch from API if we don't have initial data
      fetchMixture(id);
    }
  }, []);
  
  const fetchMixture = async (mixtureId) => {
    try {
      setLoading(true);
      setError(null);
      const data = await getMixture(mixtureId);
      setMixture(data);
    } catch (err) {
      console.error('Failed to fetch mixture:', err);
      // Try to get mock data as fallback
      const mockMixture = mockMixtures.find(m => m.id.toString() === mixtureId.toString());
      if (mockMixture) {
        setMixture(mockMixture);
      } else {
        setError('Failed to fetch mixture details. Please try again later.');
      }
    } finally {
      setLoading(false);
    }
  };
  
  // This is unlikely to happen with getStaticProps/getStaticPaths, but keeping for safety
  if (!id && !initialMixture) {
    return (
      <Layout title="Mixture Detail">
        <div className="page-header">
          <h1 className="title">Mixture Detail</h1>
        </div>
        <div className="card">
          <p>Loading mixture ID...</p>
        </div>
      </Layout>
    );
  }
  
  return (
    <Layout title={mixture ? `${mixture.name || 'Mixture'} Details` : 'Mixture Detail'}>
      <div className="page-header">
        <Link href="/mixtures">
          <span className="back-link">← Back to Mixtures</span>
        </Link>
        <h1 className="title">{mixture ? mixture.name || 'Unnamed Mixture' : 'Mixture Detail'}</h1>
      </div>
      
      {loading ? (
        <Loading />
      ) : error ? (
        <ErrorMessage message={error} onRetry={() => fetchMixture(id)} />
      ) : mixture ? (
        <>
          <div className="mixture-detail-container">
            <div className="card mixture-detail-card">
              <h2 className="subtitle">Basic Information</h2>
              
              {mixture.description && (
                <p className="mixture-description">{mixture.description}</p>
              )}
              
              <div className="property-group">
                {mixture.id && (
                  <div className="property-row">
                    <span className="property-label">Mixture ID:</span>
                    <span className="property-value">{mixture.id}</span>
                  </div>
                )}
                
                {mixture.freezing_point && (
                  <div className="property-row">
                    <span className="property-label">Freezing Point:</span>
                    <span className="property-value">{mixture.freezing_point} °C</span>
                  </div>
                )}
                
                {mixture.created_at && (
                  <div className="property-row">
                    <span className="property-label">Created:</span>
                    <span className="property-value">
                      {new Date(mixture.created_at).toLocaleDateString()}
                    </span>
                  </div>
                )}
              </div>
            </div>
            
            {mixture.components && mixture.components.length > 0 && (
              <div className="card mixture-detail-card">
                <h2 className="subtitle">Components</h2>
                <div className="components-list">
                  {mixture.components.map((component, index) => (
                    <div className="component-item" key={index}>
                      <div className="component-header">
                        <strong>{component.molecule?.name || 'Unnamed Molecule'}</strong>
                        {component.concentration && (
                          <span className="component-concentration">
                            {component.concentration} {component.concentration_unit || '%'}
                          </span>
                        )}
                      </div>
                      {component.molecule?.formula && (
                        <div className="component-formula">
                          {component.molecule.formula}
                        </div>
                      )}
                      {component.role && (
                        <div className="component-role">
                          Role: {component.role}
                        </div>
                      )}
                      {component.molecule?.id && (
                        <div className="component-link">
                          <Link href={`/molecules/${component.molecule.id}`}>
                            View Molecule Details
                          </Link>
                        </div>
                      )}
                    </div>
                  ))}
                </div>
              </div>
            )}
            
            <div className="card mixture-detail-card">
              <h2 className="subtitle">Visualization</h2>
              <p className="description">
                Interactive composition visualization will be available soon.
              </p>
              <div className="mixture-placeholder">
                Mixture composition chart placeholder
              </div>
            </div>
            
            {mixture.properties && mixture.properties.length > 0 && (
              <div className="card mixture-detail-card">
                <h2 className="subtitle">Properties</h2>
                <div className="property-group">
                  {mixture.properties.map((prop, index) => (
                    <div className="property-row" key={index}>
                      <span className="property-label">{prop.name || 'Property'}:</span>
                      <span className="property-value">{prop.value} {prop.units || ''}</span>
                    </div>
                  ))}
                </div>
              </div>
            )}
          </div>
          
          <div className="actions-container">
            <Link href="/mixtures">
              <span className="button">Back to Mixtures</span>
            </Link>
          </div>
        </>
      ) : (
        <div className="card">
          <p>No mixture found with ID: {id}</p>
          <Link href="/mixtures">
            <span className="button">Back to Mixtures</span>
          </Link>
        </div>
      )}
    </Layout>
  );
}
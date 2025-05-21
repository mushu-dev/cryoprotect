import React, { useState, useEffect } from 'react';
import Layout from '../components/Layout';
import MixtureCard from '../components/MixtureCard';
import Loading from '../components/Loading';
import ErrorMessage from '../components/ErrorMessage';
import { getMixtures } from '../utils/api';

// Mock data for initial state and fallback
const mockMixtures = [
  {
    id: 1,
    name: 'VS55',
    description: 'A cryoprotectant mixture used for organ preservation.',
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
    ]
  },
  {
    id: 2,
    name: 'Standard Cell Freezing Medium',
    description: 'Common mixture for cell line preservation.',
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
    ]
  }
];

// This function gets called at build time on server-side.
export async function getStaticProps() {
  // For the static build, we'll use mock data
  return {
    props: {
      initialMixtures: mockMixtures
    },
  }
}

export default function Mixtures({ initialMixtures }) {
  const [mixtures, setMixtures] = useState(initialMixtures);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);

  useEffect(() => {
    // We already have initial data from props, but try to fetch from API anyway
    fetchMixtures();
  }, []);

  const fetchMixtures = async () => {
    try {
      setLoading(true);
      setError(null);
      const data = await getMixtures();
      setMixtures(data || initialMixtures);
    } catch (err) {
      console.error('Failed to fetch mixtures:', err);
      // Don't show error since we already have initial data
      setMixtures(initialMixtures);
    } finally {
      setLoading(false);
    }
  };

  return (
    <Layout title="Mixtures">
      <div className="page-header">
        <h1 className="title">Mixtures</h1>
        <p className="description">
          Explore cryoprotectant mixtures and their properties.
        </p>
      </div>

      {loading ? (
        <Loading />
      ) : error ? (
        <ErrorMessage message={error} onRetry={fetchMixtures} />
      ) : mixtures.length === 0 ? (
        <div className="card">
          <p>No mixtures found in the database.</p>
        </div>
      ) : (
        <>
          <div className="mixtures-grid">
            {mixtures.map((mixture) => (
              <MixtureCard key={mixture.id} mixture={mixture} />
            ))}
          </div>

          <div className="card" style={{ marginTop: '30px' }}>
            <h2 className="subtitle">Coming Features</h2>
            <ul className="list">
              <li>Mixture Composition Analysis</li>
              <li>Performance Comparison Tools</li>
              <li>Interactive Composition Charts</li>
              <li>Freezing Point Depression Data</li>
              <li>Application-specific Recommendations</li>
            </ul>
          </div>
        </>
      )}
    </Layout>
  );
}
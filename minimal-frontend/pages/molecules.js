import React, { useState, useEffect } from 'react';
import Layout from '../components/Layout';
import MoleculeCard from '../components/MoleculeCard';
import Loading from '../components/Loading';
import ErrorMessage from '../components/ErrorMessage';
import { getMolecules } from '../utils/api';

// Mock data for initial state and fallback
const mockMolecules = [
  {
    id: 1,
    name: 'Glycerol',
    formula: 'C3H8O3',
    pubchem_cid: '753',
    molecular_weight: 92.09,
    is_cryoprotectant: true,
    description: 'A common cryoprotectant used in various applications.'
  },
  {
    id: 2,
    name: 'Dimethyl Sulfoxide (DMSO)',
    formula: 'C2H6OS',
    pubchem_cid: '679',
    molecular_weight: 78.13,
    is_cryoprotectant: true,
    description: 'A widely used penetrating cryoprotectant.'
  },
  {
    id: 3,
    name: 'Ethylene Glycol',
    formula: 'C2H6O2',
    pubchem_cid: '174',
    molecular_weight: 62.07,
    is_cryoprotectant: true,
    description: 'Used in cryopreservation of embryos and tissues.'
  },
  {
    id: 4,
    name: 'Propylene Glycol',
    formula: 'C3H8O2',
    pubchem_cid: '1030',
    molecular_weight: 76.09,
    is_cryoprotectant: true,
    description: 'Used as a cryoprotectant for various biological materials.'
  }
];

// This function gets called at build time on server-side.
// It won't be called on client-side, so you can even do
// direct database queries in it.
export async function getStaticProps() {
  // For the static build, we'll use mock data
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
    // We already have initial data from props, but try to fetch from API anyway
    fetchMolecules();
  }, []);

  const fetchMolecules = async () => {
    try {
      setLoading(true);
      setError(null);
      const data = await getMolecules();
      setMolecules(data || initialMolecules);
    } catch (err) {
      console.error('Failed to fetch molecules:', err);
      // Don't show error since we already have initial data
      setMolecules(initialMolecules);
    } finally {
      setLoading(false);
    }
  };

  return (
    <Layout title="Molecules">
      <div className="page-header">
        <h1 className="title">Molecules</h1>
        <p className="description">
          View and explore molecules in the CryoProtect database.
        </p>
      </div>

      {loading ? (
        <Loading />
      ) : error ? (
        <ErrorMessage message={error} onRetry={fetchMolecules} />
      ) : molecules.length === 0 ? (
        <div className="card">
          <p>No molecules found in the database.</p>
        </div>
      ) : (
        <>
          <div className="molecules-grid">
            {molecules.map((molecule) => (
              <MoleculeCard key={molecule.id} molecule={molecule} />
            ))}
          </div>

          <div className="card" style={{ marginTop: '30px' }}>
            <h2 className="subtitle">Coming Features</h2>
            <ul className="list">
              <li>3D Molecular Visualization</li>
              <li>Chemical Property Exploration</li>
              <li>Structure Search</li>
              <li>Predictive Property Analysis</li>
              <li>Molecular Comparison Tools</li>
            </ul>
          </div>
        </>
      )}
    </Layout>
  );
}
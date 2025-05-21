import { useState } from 'react';
import Link from 'next/link';
import { useMolecules } from '../../convex/hooks';
import MoleculeCard from '../../../components/MoleculeCard';
import Loading from '../../../components/Loading';
import ErrorMessage from '../../../components/ErrorMessage';

/**
 * ConvexMoleculesList component to display molecules from Convex
 * This component adapts Convex data to the format expected by MoleculeCard
 */
export function ConvexMoleculesList() {
  const molecules = useMolecules();
  const [error, setError] = useState(null);

  if (molecules === undefined) {
    return <Loading />;
  }

  if (error) {
    return <ErrorMessage message={error} onRetry={() => setError(null)} />;
  }

  if (molecules.length === 0) {
    return (
      <div className="card">
        <h2 className="subtitle">No Molecules Found</h2>
        <p>No molecules were found in the Convex database.</p>
      </div>
    );
  }

  // Adapt Convex data to the format expected by MoleculeCard
  const adaptedMolecules = molecules.map(molecule => ({
    id: molecule._id,
    name: molecule.name,
    pubchem_cid: molecule.pubchemCid,
    formula: molecule.formula,
    molecular_weight: molecule.molecularWeight,
    is_cryoprotectant: molecule.isCryoprotectant,
    // Add any other fields needed by MoleculeCard
  }));

  return (
    <div className="molecules-grid">
      {adaptedMolecules.map((molecule) => (
        <MoleculeCard key={molecule.id} molecule={molecule} />
      ))}
    </div>
  );
}
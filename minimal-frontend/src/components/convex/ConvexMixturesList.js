import { useState } from 'react';
import Link from 'next/link';
import { useMixtures } from '../../convex/hooks';
import MixtureCard from '../../../components/MixtureCard';
import Loading from '../../../components/Loading';
import ErrorMessage from '../../../components/ErrorMessage';

/**
 * ConvexMixturesList component to display mixtures from Convex
 * This component adapts Convex data to the format expected by MixtureCard
 */
export function ConvexMixturesList() {
  const mixtures = useMixtures();
  const [error, setError] = useState(null);

  if (mixtures === undefined) {
    return <Loading />;
  }

  if (error) {
    return <ErrorMessage message={error} onRetry={() => setError(null)} />;
  }

  if (mixtures.length === 0) {
    return (
      <div className="card">
        <h2 className="subtitle">No Mixtures Found</h2>
        <p>No mixtures were found in the Convex database.</p>
      </div>
    );
  }

  // Adapt Convex data to the format expected by MixtureCard
  const adaptedMixtures = mixtures.map(mixture => ({
    id: mixture._id,
    name: mixture.name,
    description: mixture.description,
    freezing_point: mixture.freezingPoint,
    components: mixture.components || [],
    // Add any other fields needed by MixtureCard
  }));

  return (
    <div className="mixtures-grid">
      {adaptedMixtures.map((mixture) => (
        <MixtureCard key={mixture.id} mixture={mixture} />
      ))}
    </div>
  );
}
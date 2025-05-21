import React from 'react';
import type { Metadata } from 'next';
import ExperimentsList from '../../features/experiments/components/ExperimentsList';

export const metadata: Metadata = {
  title: 'Experiments | CryoProtect',
  description: 'Design, track, and analyze cryopreservation experiments with detailed protocols',
};

export default function ExperimentsPage() {
  return (
    <div className="container mx-auto px-4 py-8">
      <ExperimentsList />
    </div>
  );
}
import React from 'react';
import Head from 'next/head';
import ExperimentsList from '../../features/experiments/components/ExperimentsList';

export default function ExperimentsPage() {
  return (
    <>
      <Head>
        <title>Experiments | CryoProtect</title>
        <meta name="description" content="Design, track, and analyze cryopreservation experiments with detailed protocols" />
      </Head>
    
      <div className="container mx-auto px-4 py-8">
        <ExperimentsList />
      </div>
    </>
  );
}
import React from 'react';
import { Metadata } from 'next';
import { ExperimentsList } from '../../features/experiments/components/experiments-list';
import { Button } from '../../components/ui/button';
import Link from 'next/link';

export const metadata: Metadata = {
  title: 'Experiments | CryoProtect',
  description: 'Design, track, and analyze cryopreservation experiments with detailed protocols',
};

export default function ExperimentsPage() {
  return (
    <div className="container mx-auto px-4 py-8">
      <div className="flex flex-col md:flex-row justify-between items-start md:items-center mb-8">
        <div>
          <h1 className="text-3xl font-bold mb-2">Experiments</h1>
          <p className="text-muted-foreground">
            Design, track, and analyze cryopreservation experiments
          </p>
        </div>
        
        <div className="mt-4 md:mt-0">
          <Link href="/experiments/create" passHref>
            <Button>Create New Experiment</Button>
          </Link>
        </div>
      </div>

      <ExperimentsList 
        initialParams={{ page: 1, per_page: 12 }}
        onExperimentClick={(id) => `/experiments/${id}`}
        onCreateNew={() => "/experiments/create"}
      />
    </div>
  );
}
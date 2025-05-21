import React from 'react';
import type { Metadata } from 'next';
import Link from 'next/link';
import { ExperimentCreationWizard } from '../../../features/experiments/components/experiment-creation-wizard';

export const metadata: Metadata = {
  title: 'Create Experiment | CryoProtect',
  description: 'Create a new cryopreservation experiment with detailed protocols and parameters',
};

export default function CreateExperimentPage() {
  return (
    <div className="container mx-auto px-4 py-8">
      <div className="mb-6">
        <Link 
          href="/experiments"
          className="text-primary hover:underline flex items-center"
        >
          <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" className="mr-2">
            <path d="M19 12H5"/>
            <path d="M12 19l-7-7 7-7"/>
          </svg>
          Back to Experiments
        </Link>
      </div>
      
      <div className="mb-6">
        <h1 className="text-3xl font-bold mb-2">Create New Experiment</h1>
        <p className="text-muted-foreground">
          Design a new cryopreservation experiment with protocol selection, parameters, and success criteria
        </p>
      </div>
      
      <ExperimentCreationWizard />
    </div>
  );
}
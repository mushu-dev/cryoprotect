import React from 'react';
import Link from 'next/link';
import { notFound } from 'next/navigation';
import type { Metadata } from 'next';
import ExperimentDetail from '../../../features/experiments/components/ExperimentDetail';
import { getExperimentById } from '../../../features/experiments/actions/experiment-actions';

interface ExperimentDetailPageProps {
  params: {
    id: string;
  };
}

export async function generateMetadata({ params }: ExperimentDetailPageProps): Promise<Metadata> {
  const experiment = await getExperimentById(params.id);
  
  if (!experiment) {
    return {
      title: 'Experiment Not Found | CryoProtect',
      description: 'The requested experiment could not be found',
    };
  }
  
  return {
    title: `${experiment.title} | CryoProtect`,
    description: experiment.description,
  };
}

export default async function ExperimentDetailPage({ params }: ExperimentDetailPageProps) {
  const experiment = await getExperimentById(params.id);
  
  if (!experiment) {
    notFound();
  }
  
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
      
      <ExperimentDetail experiment={experiment} />
    </div>
  );
}
import React, { useEffect, useState } from 'react';
import Link from 'next/link';
import Head from 'next/head';
import { useRouter } from 'next/router';
import ExperimentDetail from '../../features/experiments/components/ExperimentDetail';
import useExperimentData from '../../features/experiments/hooks/useExperimentData';

export default function ExperimentDetailPage() {
  const router = useRouter();
  const { id } = router.query;
  const { fetchExperiment } = useExperimentData();
  const [experiment, setExperiment] = useState(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState(null);

  // Fetch experiment data when ID is available
  useEffect(() => {
    async function loadExperiment() {
      if (id) {
        try {
          setLoading(true);
          const data = await fetchExperiment(id);
          setExperiment(data);
          setError(null);
        } catch (err) {
          setError('Failed to load experiment details');
          console.error('Error loading experiment:', err);
        } finally {
          setLoading(false);
        }
      }
    }

    loadExperiment();
  }, [id, fetchExperiment]);

  if (router.isFallback || loading) {
    return (
      <div className="container mx-auto px-4 py-8">
        <div className="py-12 text-center">
          <div className="animate-spin rounded-full h-12 w-12 border-b-2 border-primary mx-auto"></div>
          <p className="mt-4 text-gray-500">Loading experiment details...</p>
        </div>
      </div>
    );
  }

  if (error || !experiment) {
    return (
      <div className="container mx-auto px-4 py-8">
        <div className="py-12 text-center">
          <div className="text-red-500 mb-2">Error loading experiment</div>
          <p className="text-gray-500">The requested experiment could not be found</p>
          <Link href="/experiments">
            <a className="mt-4 inline-block text-blue-600 hover:underline">
              Return to Experiments
            </a>
          </Link>
        </div>
      </div>
    );
  }

  return (
    <>
      <Head>
        <title>{experiment.title} | CryoProtect</title>
        <meta name="description" content={experiment.description} />
      </Head>
    
      <div className="container mx-auto px-4 py-8">
        <div className="mb-6">
          <Link href="/experiments">
            <a className="text-primary hover:underline flex items-center">
              <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" className="mr-2">
                <path d="M19 12H5"/>
                <path d="M12 19l-7-7 7-7"/>
              </svg>
              Back to Experiments
            </a>
          </Link>
        </div>
        
        <ExperimentDetail experiment={experiment} />
      </div>
    </>
  );
}
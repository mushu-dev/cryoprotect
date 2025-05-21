import React from 'react';
import Link from 'next/link';

export default function ExperimentNotFound() {
  return (
    <div className="container mx-auto px-4 py-8">
      <div className="py-12 text-center">
        <div className="text-red-500 mb-2">Experiment Not Found</div>
        <p className="text-gray-500">The requested experiment could not be found</p>
        <Link 
          href="/experiments"
          className="mt-4 inline-block text-blue-600 hover:underline"
        >
          Return to Experiments
        </Link>
      </div>
    </div>
  );
}
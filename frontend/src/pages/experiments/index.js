import React from 'react';
import Link from 'next/link';
import Head from 'next/head';

export default function ExperimentsPage() {
  // This is a simplified implementation for Next.js 12
  return (
    <>
      <Head>
        <title>Experiments | CryoProtect</title>
        <meta name="description" content="Design, track, and analyze cryopreservation experiments with detailed protocols" />
      </Head>
    
      <div className="container mx-auto px-4 py-8">
        <div className="flex flex-col md:flex-row justify-between items-start md:items-center mb-8">
          <div>
            <h1 className="text-3xl font-bold mb-2">Experiments</h1>
            <p className="text-muted-foreground">
              Design, track, and analyze cryopreservation experiments
            </p>
          </div>
          
          <div className="mt-4 md:mt-0">
            <Link href="/experiments/create">
              <a className="inline-flex items-center justify-center rounded-md bg-primary px-4 py-2 text-sm font-medium text-primary-foreground shadow transition-colors hover:bg-primary/90">
                Create New Experiment
              </a>
            </Link>
          </div>
        </div>

        {/* Sample experiments list for demonstration */}
        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
          {[1, 2, 3, 4, 5, 6].map((id) => (
            <div key={id} className="rounded-lg border bg-card p-6 shadow-sm">
              <div className="flex items-center justify-between mb-4">
                <h3 className="text-lg font-semibold">Experiment #{id}</h3>
                <span className="inline-flex items-center rounded-full bg-primary/10 px-2.5 py-0.5 text-xs font-medium text-primary">
                  {id % 2 === 0 ? 'Completed' : 'In Progress'}
                </span>
              </div>
              <p className="text-muted-foreground mb-4">
                Sample experiment with cryoprotectant testing and viability analysis.
              </p>
              <div className="flex justify-between items-center mt-4">
                <span className="text-sm text-muted-foreground">May {id + 10}, 2025</span>
                <Link href={`/experiments/${id}`}>
                  <a className="text-primary hover:underline">View Details</a>
                </Link>
              </div>
            </div>
          ))}
        </div>
      </div>
    </>
  );
}
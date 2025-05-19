import React from 'react';
import Link from 'next/link';
import Head from 'next/head';

export default function Home() {
  return (
    <>
      <Head>
        <title>CryoProtect - A platform for cryoprotectant analysis</title>
        <meta name="description" content="A platform for cryoprotectant analysis and experiment management" />
      </Head>
      
      <div className="container mx-auto px-4 py-8">
        <section className="mb-12">
          <h1 className="text-4xl font-bold mb-4">CryoProtect</h1>
          <p className="text-xl text-gray-700 dark:text-gray-300">
            A platform for cryoprotectant analysis and experiment management
          </p>
        </section>
        
        <section className="mb-12">
          <h2 className="text-3xl font-bold mb-4">Molecule Database</h2>
          <p className="mb-6 text-gray-700 dark:text-gray-300">
            Explore our comprehensive database of cryoprotectant molecules and their properties.
          </p>
          <Link href="/molecules">
            <a className="text-blue-600 hover:underline font-medium">View Molecules</a>
          </Link>
        </section>
        
        <section className="mb-12">
          <h2 className="text-3xl font-bold mb-4">Mixtures</h2>
          <p className="mb-6 text-gray-700 dark:text-gray-300">
            Discover optimized cryoprotectant mixtures and their performance characteristics.
          </p>
          <Link href="/mixtures">
            <a className="text-blue-600 hover:underline font-medium">View Mixtures</a>
          </Link>
        </section>
        
        <section className="mb-12">
          <h2 className="text-3xl font-bold mb-4">Experiments</h2>
          <p className="mb-6 text-gray-700 dark:text-gray-300">
            Design, track, and analyze cryopreservation experiments with detailed protocols.
          </p>
          <Link href="/experiments">
            <a className="text-blue-600 hover:underline font-medium">Manage Experiments</a>
          </Link>
        </section>
        
        <section className="mb-12">
          <h2 className="text-3xl font-bold mb-4">Protocols</h2>
          <p className="mb-6 text-gray-700 dark:text-gray-300">
            Create and manage standardized protocols for reproducible cryopreservation procedures.
          </p>
          <Link href="/protocols">
            <a className="text-blue-600 hover:underline font-medium">Browse Protocols</a>
          </Link>
        </section>
        
        <footer className="text-sm text-gray-500 dark:text-gray-400 mt-12">
          © 2025 CryoProtect. All rights reserved.
        </footer>
      </div>
    </>
  );
}
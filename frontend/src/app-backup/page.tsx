import React from 'react';
import Link from 'next/link';

export default function HomePage() {
  return (
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
        <Link href="/molecules" className="text-blue-600 hover:underline font-medium">
          <button className="px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700">
            View Molecules
          </button>
        </Link>
      </section>
      
      <section className="mb-12">
        <h2 className="text-3xl font-bold mb-4">Mixtures</h2>
        <p className="mb-6 text-gray-700 dark:text-gray-300">
          Discover optimized cryoprotectant mixtures and their performance characteristics.
        </p>
        <Link href="/mixtures" className="text-blue-600 hover:underline font-medium">
          <button className="px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700">
            View Mixtures
          </button>
        </Link>
      </section>
      
      <section className="mb-12">
        <h2 className="text-3xl font-bold mb-4">Experiments</h2>
        <p className="mb-6 text-gray-700 dark:text-gray-300">
          Design, track, and analyze cryopreservation experiments with detailed protocols.
        </p>
        <Link href="/experiments" className="text-blue-600 hover:underline font-medium">
          <button className="px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700">
            Manage Experiments
          </button>
        </Link>
      </section>
      
      <section className="mb-12">
        <h2 className="text-3xl font-bold mb-4">Protocols</h2>
        <p className="mb-6 text-gray-700 dark:text-gray-300">
          Create and manage standardized protocols for reproducible cryopreservation procedures.
        </p>
        <Link href="/protocols" className="text-blue-600 hover:underline font-medium">
          <button className="px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700">
            Browse Protocols
          </button>
        </Link>
      </section>
      
      <footer className="text-sm text-gray-500 dark:text-gray-400 mt-12">
        © 2025 CryoProtect. All rights reserved.
      </footer>
    </div>
  );
}
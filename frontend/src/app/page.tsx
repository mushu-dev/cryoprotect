import React from 'react';
import Link from 'next/link';
import { Button } from '../components/ui/button';

export default function Home() {
  return (
    <div className="container mx-auto px-4 py-12">
      <section className="flex flex-col items-center justify-center text-center py-12">
        <h1 className="text-4xl font-bold mb-4">CryoProtect</h1>
        <p className="text-xl mb-8">A platform for cryoprotectant analysis and experiment management</p>
        
        <div className="grid grid-cols-1 md:grid-cols-2 gap-6 mt-8 max-w-4xl w-full">
          <div className="bg-card rounded-lg p-6 shadow-md">
            <h2 className="text-2xl font-semibold mb-4">Molecule Database</h2>
            <p className="mb-4">Explore our comprehensive database of cryoprotectant molecules and their properties.</p>
            <Link href="/molecules" passHref>
              <Button className="w-full">View Molecules</Button>
            </Link>
          </div>
          
          <div className="bg-card rounded-lg p-6 shadow-md">
            <h2 className="text-2xl font-semibold mb-4">Mixtures</h2>
            <p className="mb-4">Discover optimized cryoprotectant mixtures and their performance characteristics.</p>
            <Link href="/mixtures" passHref>
              <Button className="w-full">View Mixtures</Button>
            </Link>
          </div>
          
          <div className="bg-card rounded-lg p-6 shadow-md">
            <h2 className="text-2xl font-semibold mb-4">Experiments</h2>
            <p className="mb-4">Design, track, and analyze cryopreservation experiments with detailed protocols.</p>
            <Link href="/experiments" passHref>
              <Button className="w-full">Manage Experiments</Button>
            </Link>
          </div>
          
          <div className="bg-card rounded-lg p-6 shadow-md">
            <h2 className="text-2xl font-semibold mb-4">Protocols</h2>
            <p className="mb-4">Create and manage standardized protocols for reproducible cryopreservation procedures.</p>
            <Link href="/protocols" passHref>
              <Button className="w-full">Browse Protocols</Button>
            </Link>
          </div>
        </div>
      </section>
    </div>
  );
}
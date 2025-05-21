import React from 'react';
import Head from 'next/head';
import { useEffect, useState } from 'react';

export default function MoleculesPage() {
  const [isLoading, setIsLoading] = useState(true);
  const [molecules, setMolecules] = useState([]);
  const [error, setError] = useState(null);

  useEffect(() => {
    async function fetchMolecules() {
      try {
        setIsLoading(true);
        const response = await fetch(`${process.env.NEXT_PUBLIC_API_URL}/molecules`);
        if (!response.ok) {
          throw new Error(`API error: ${response.status}`);
        }
        const data = await response.json();
        setMolecules(data);
      } catch (err) {
        console.error('Error fetching molecules:', err);
        setError('Unable to load molecules. Please try again later.');
      } finally {
        setIsLoading(false);
      }
    }

    fetchMolecules();
  }, []);

  return (
    <>
      <Head>
        <title>Molecules - CryoProtect</title>
        <meta name="description" content="Browse and search for cryoprotectant molecules" />
      </Head>

      <div className="container mx-auto px-4 py-8">
        <h1 className="text-3xl font-bold mb-6">Molecules</h1>
        
        {error && (
          <div className="bg-red-100 border border-red-400 text-red-700 px-4 py-3 rounded mb-4">
            {error}
          </div>
        )}
        
        {isLoading ? (
          <div className="text-center py-8">
            <p>Loading molecules...</p>
          </div>
        ) : (
          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
            {molecules && molecules.length > 0 ? (
              molecules.map((molecule) => (
                <div key={molecule.id} className="border rounded-lg p-4 shadow-sm hover:shadow-md transition-shadow">
                  <h2 className="text-xl font-semibold mb-2">{molecule.name}</h2>
                  <p className="text-gray-600 mb-3">{molecule.formula || 'No formula available'}</p>
                  {molecule.smiles && (
                    <div className="text-sm text-gray-500 mb-2">
                      SMILES: {molecule.smiles.length > 30 ? `${molecule.smiles.substring(0, 30)}...` : molecule.smiles}
                    </div>
                  )}
                  {molecule.properties && (
                    <div className="text-sm text-gray-500">
                      Properties: {Object.keys(molecule.properties).length}
                    </div>
                  )}
                </div>
              ))
            ) : (
              <div className="col-span-3 text-center py-8">
                <p>No molecules found. Please try again later.</p>
              </div>
            )}
          </div>
        )}
      </div>
    </>
  );
}
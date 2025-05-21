import React from 'react';
import Head from 'next/head';
import { useQuery } from 'convex/react';
import { api } from '../../../convex/_generated/api';

export default function MoleculesPageConvex() {
  // Use Convex query to fetch molecules
  const molecules = useQuery(api.molecules.getAllMolecules, { limit: 50 });
  
  // Determine if we're still loading
  const isLoading = molecules === undefined;
  
  return (
    <>
      <Head>
        <title>Molecules - CryoProtect (Convex)</title>
        <meta name="description" content="Browse and search for cryoprotectant molecules using Convex" />
      </Head>

      <div className="container mx-auto px-4 py-8">
        <h1 className="text-3xl font-bold mb-6">Molecules (Convex-powered)</h1>
        
        {isLoading ? (
          <div className="text-center py-8">
            <p>Loading molecules from Convex...</p>
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
                  
                  {/* Show properties */}
                  {molecule.properties && Object.keys(molecule.properties).length > 0 && (
                    <div className="mt-4">
                      <h3 className="text-lg font-medium mb-2">Properties</h3>
                      <ul className="text-sm">
                        {Object.entries(molecule.properties).map(([key, prop]) => (
                          <li key={key} className="grid grid-cols-2 gap-2 mb-1">
                            <span className="text-gray-600">{prop.displayName || key}:</span>
                            <span className="font-medium">
                              {prop.value} {prop.units || ''}
                            </span>
                          </li>
                        ))}
                      </ul>
                    </div>
                  )}
                </div>
              ))
            ) : (
              <div className="col-span-3 text-center py-8">
                <p>No molecules found in the Convex database.</p>
              </div>
            )}
          </div>
        )}
      </div>
    </>
  );
}
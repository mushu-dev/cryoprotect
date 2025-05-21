import React from 'react';
import Head from 'next/head';
import { useEffect, useState } from 'react';

export default function MixturesPage() {
  const [isLoading, setIsLoading] = useState(true);
  const [mixtures, setMixtures] = useState([]);
  const [error, setError] = useState(null);

  useEffect(() => {
    async function fetchMixtures() {
      try {
        setIsLoading(true);
        const response = await fetch(`${process.env.NEXT_PUBLIC_API_URL}/mixtures`);
        if (!response.ok) {
          throw new Error(`API error: ${response.status}`);
        }
        const data = await response.json();
        setMixtures(data);
      } catch (err) {
        console.error('Error fetching mixtures:', err);
        setError('Unable to load mixtures. Please try again later.');
      } finally {
        setIsLoading(false);
      }
    }

    fetchMixtures();
  }, []);

  return (
    <>
      <Head>
        <title>Mixtures - CryoProtect</title>
        <meta name="description" content="Browse and search for cryoprotectant mixtures" />
      </Head>

      <div className="container mx-auto px-4 py-8">
        <h1 className="text-3xl font-bold mb-6">Mixtures</h1>
        
        {error && (
          <div className="bg-red-100 border border-red-400 text-red-700 px-4 py-3 rounded mb-4">
            {error}
          </div>
        )}
        
        {isLoading ? (
          <div className="text-center py-8">
            <p>Loading mixtures...</p>
          </div>
        ) : (
          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
            {mixtures && mixtures.length > 0 ? (
              mixtures.map((mixture) => (
                <div key={mixture.id} className="border rounded-lg p-4 shadow-sm hover:shadow-md transition-shadow">
                  <h2 className="text-xl font-semibold mb-2">{mixture.name}</h2>
                  <p className="text-gray-600 mb-3">{mixture.description || 'No description available'}</p>
                  <div className="text-sm text-gray-500">
                    Components: {mixture.components?.length || 0}
                  </div>
                </div>
              ))
            ) : (
              <div className="col-span-3 text-center py-8">
                <p>No mixtures found. Please try again later.</p>
              </div>
            )}
          </div>
        )}
      </div>
    </>
  );
}
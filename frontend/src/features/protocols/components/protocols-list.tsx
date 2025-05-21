import React, { useState, useEffect } from 'react';
import Link from 'next/link';
import { getAllProtocols } from '../actions/protocol-actions';

/**
 * Component for displaying a list of protocols with filtering capability
 */
export default function ProtocolsList() {
  const [activeTab, setActiveTab] = useState('all');
  const [protocols, setProtocols] = useState<any[]>([]);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [search, setSearch] = useState('');

  // Fetch protocols on component mount
  useEffect(() => {
    const loadProtocols = async () => {
      try {
        setLoading(true);
        const data = await getAllProtocols();
        setProtocols(data);
        setError(null);
      } catch (err) {
        console.error('Error loading protocols:', err);
        setError('Failed to load protocols. Please try again later.');
      } finally {
        setLoading(false);
      }
    };

    loadProtocols();
  }, []);

  // Filter protocols based on tab and search query
  const filteredProtocols = protocols.filter(protocol => {
    // Apply tab filter
    if (activeTab === 'my' && protocol.author.id !== 'user1') {
      return false;
    }
    if (activeTab === 'templates' && !protocol.is_template) {
      return false;
    }

    // Apply search filter (if search is not empty)
    if (search.trim() !== '') {
      const searchLower = search.toLowerCase();
      return (
        protocol.name.toLowerCase().includes(searchLower) ||
        protocol.description.toLowerCase().includes(searchLower)
      );
    }

    return true;
  });

  if (loading) {
    return (
      <div className="py-8 text-center">
        <div className="animate-spin rounded-full h-12 w-12 border-b-2 border-primary mx-auto"></div>
        <p className="mt-4 text-gray-500">Loading protocols...</p>
      </div>
    );
  }

  if (error) {
    return (
      <div className="py-8 text-center">
        <div className="text-red-500 mb-2">Error loading protocols</div>
        <p className="text-gray-500">{error}</p>
      </div>
    );
  }

  return (
    <>
      <div className="flex flex-col md:flex-row justify-between items-start md:items-center mb-8">
        <div>
          <h1 className="text-3xl font-bold mb-2">Protocols</h1>
          <p className="text-muted-foreground">
            Standardized protocols for reproducible cryopreservation procedures
          </p>
        </div>
        
        <div className="mt-4 md:mt-0">
          <Link 
            href="/protocols/create"
            className="inline-flex items-center justify-center rounded-md bg-primary px-4 py-2 text-sm font-medium text-primary-foreground shadow transition-colors hover:bg-primary/90"
          >
            Create New Protocol
          </Link>
        </div>
      </div>
    
      <div className="mb-8">
        <div className="max-w-md mx-auto mb-8">
          <input 
            type="search" 
            placeholder="Search protocols..." 
            className="w-full px-4 py-2 border rounded-md"
            value={search}
            onChange={(e) => setSearch(e.target.value)}
          />
        </div>

        <div className="border-b mb-6">
          <div className="flex flex-wrap -mb-px">
            <button
              className={`inline-block p-4 border-b-2 ${
                activeTab === 'all' 
                  ? 'border-primary text-primary' 
                  : 'border-transparent text-muted-foreground hover:text-foreground hover:border-muted'
              }`}
              onClick={() => setActiveTab('all')}
            >
              All Protocols
            </button>
            <button
              className={`inline-block p-4 border-b-2 ${
                activeTab === 'my' 
                  ? 'border-primary text-primary' 
                  : 'border-transparent text-muted-foreground hover:text-foreground hover:border-muted'
              }`}
              onClick={() => setActiveTab('my')}
            >
              My Protocols
            </button>
            <button
              className={`inline-block p-4 border-b-2 ${
                activeTab === 'templates' 
                  ? 'border-primary text-primary' 
                  : 'border-transparent text-muted-foreground hover:text-foreground hover:border-muted'
              }`}
              onClick={() => setActiveTab('templates')}
            >
              Templates
            </button>
          </div>
        </div>
        
        {filteredProtocols.length === 0 ? (
          <div className="py-12 text-center border rounded-lg bg-gray-50">
            <p className="text-gray-500 mb-2">No protocols found</p>
            <button
              onClick={() => { setSearch(''); setActiveTab('all'); }}
              className="text-blue-600 hover:text-blue-800 hover:underline text-sm"
            >
              Clear filters
            </button>
          </div>
        ) : (
          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
            {filteredProtocols.map(protocol => (
              <div key={protocol.id} className="rounded-lg border bg-card shadow-sm overflow-hidden">
                <div className="p-6">
                  <div className="flex flex-wrap items-center gap-2 mb-2">
                    <h3 className="text-lg font-semibold">{protocol.name}</h3>
                    <span className="inline-flex items-center rounded-full bg-muted px-2.5 py-0.5 text-xs font-medium">
                      v{protocol.version}
                    </span>
                    {protocol.is_template && (
                      <span className="inline-flex items-center rounded-full border px-2.5 py-0.5 text-xs font-medium">
                        Template
                      </span>
                    )}
                  </div>
                  <p className="text-muted-foreground mb-4 text-sm">{protocol.description}</p>
                  <div className="flex text-sm text-muted-foreground mb-4">
                    <div className="mr-4">{protocol.steps.length} steps</div>
                    <div>{protocol.duration}</div>
                  </div>
                  <div className="flex justify-between items-center">
                    <div className="text-xs text-muted-foreground">
                      By {protocol.author.name}
                    </div>
                    {protocol.is_template ? (
                      <button className="text-primary text-sm hover:underline">
                        Use Template
                      </button>
                    ) : (
                      <Link href={`/protocols/${protocol.id}`}>
                        <span className="text-primary text-sm hover:underline">View Details</span>
                      </Link>
                    )}
                  </div>
                </div>
              </div>
            ))}
          </div>
        )}
      </div>
    </>
  );
}
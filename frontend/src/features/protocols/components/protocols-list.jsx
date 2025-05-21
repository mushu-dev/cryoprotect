import React, { useState, useEffect } from 'react';
import Link from 'next/link';
import { useProtocols } from '../hooks/use-protocols';

/**
 * Protocols List Component
 * Displays a list of available protocols
 */
export default function ProtocolsList() {
  const [activeTab, setActiveTab] = useState('all');
  const [protocols, setProtocols] = useState([]);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState(null);
  const [search, setSearch] = useState('');
  
  // Generate mock protocols for demonstration
  useEffect(() => {
    const mockProtocols = [
      {
        id: '1',
        name: 'Standard Cryopreservation Protocol',
        description: 'A standard protocol for cryopreservation of biological samples.',
        category: 'Cryopreservation',
        steps: [
          { title: 'Sample Preparation', description: 'Prepare the sample for cryopreservation' },
          { title: 'Add Cryoprotectant', description: 'Add appropriate cryoprotectant solution' },
          { title: 'Controlled Cooling', description: 'Cool the sample at controlled rate' },
          { title: 'Storage', description: 'Store in appropriate conditions' }
        ],
        createdAt: new Date(Date.now() - 60 * 24 * 60 * 60 * 1000),
        author: 'Dr. Smith'
      },
      {
        id: '2',
        name: 'Vitrification Protocol',
        description: 'Rapid cooling protocol to achieve vitrification of samples.',
        category: 'Vitrification',
        steps: [
          { title: 'Sample Preparation', description: 'Prepare the sample for vitrification' },
          { title: 'Add Cryoprotectant Mixture', description: 'Add vitrification solution' },
          { title: 'Rapid Cooling', description: 'Plunge into liquid nitrogen' },
          { title: 'Storage', description: 'Store in liquid nitrogen' }
        ],
        createdAt: new Date(Date.now() - 45 * 24 * 60 * 60 * 1000),
        author: 'Dr. Johnson'
      },
      {
        id: '3',
        name: 'Slow Freezing Protocol',
        description: 'Controlled slow cooling protocol for sensitive samples.',
        category: 'Slow Freezing',
        steps: [
          { title: 'Prepare Samples', description: 'Prepare samples for freezing' },
          { title: 'Add Cryoprotectant', description: 'Add cryoprotectant solution' },
          { title: 'Load Freezing Container', description: 'Place samples in freezing container' },
          { title: 'Program Freezer', description: 'Set appropriate cooling rate' },
          { title: 'Initiate Freezing', description: 'Start freezing process' },
          { title: 'Transfer to Storage', description: 'Transfer to long-term storage' }
        ],
        createdAt: new Date(Date.now() - 30 * 24 * 60 * 60 * 1000),
        author: 'Dr. Williams'
      },
      {
        id: '4',
        name: 'Cell Line Cryopreservation',
        description: 'Specialized protocol for cryopreservation of cell lines.',
        category: 'Cell Preservation',
        steps: [
          { title: 'Harvest Cells', description: 'Collect cells at optimal density' },
          { title: 'Centrifugation', description: 'Pellet cells by centrifugation' },
          { title: 'Resuspend in Media', description: 'Resuspend in freezing media' },
          { title: 'Aliquot', description: 'Distribute into cryovials' },
          { title: 'Controlled Freezing', description: 'Freeze at -1°C/minute' },
          { title: 'Transfer to Storage', description: 'Move to liquid nitrogen storage' }
        ],
        createdAt: new Date(Date.now() - 15 * 24 * 60 * 60 * 1000),
        author: 'Dr. Brown'
      }
    ];
    
    setTimeout(() => {
      setProtocols(mockProtocols);
      setLoading(false);
    }, 1000);
  }, []);
  
  // Filter protocols based on active tab and search term
  const filteredProtocols = protocols.filter(protocol => {
    // Filter by tab
    if (activeTab !== 'all' && protocol.category.toLowerCase() !== activeTab.toLowerCase()) {
      return false;
    }
    
    // Filter by search term
    if (search && !protocol.name.toLowerCase().includes(search.toLowerCase()) && 
        !protocol.description.toLowerCase().includes(search.toLowerCase())) {
      return false;
    }
    
    return true;
  });
  
  if (loading) {
    return (
      <div className="flex justify-center items-center h-64">
        <div className="text-center">
          <div className="animate-spin rounded-full h-12 w-12 border-t-2 border-b-2 border-primary mx-auto mb-4"></div>
          <p className="text-muted-foreground">Loading protocols...</p>
        </div>
      </div>
    );
  }
  
  if (error) {
    return (
      <div className="p-4 border border-red-200 bg-red-50 dark:bg-red-900/20 dark:border-red-800 rounded-md">
        <p className="text-red-600 dark:text-red-400">Error loading protocols: {error}</p>
      </div>
    );
  }
  
  return (
    <div>
      {/* Search and Filters */}
      <div className="mb-6">
        <div className="flex flex-col sm:flex-row gap-4">
          <div className="relative flex-1">
            <svg className="absolute top-3 left-3 h-5 w-5 text-gray-400" xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" stroke="currentColor">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M21 21l-6-6m2-5a7 7 0 11-14 0 7 7 0 0114 0z" />
            </svg>
            <input
              type="text"
              placeholder="Search protocols..."
              className="pl-10 pr-4 py-2 w-full border border-gray-300 dark:border-gray-600 rounded-md shadow-sm focus:ring-primary focus:border-primary dark:bg-gray-700"
              value={search}
              onChange={(e) => setSearch(e.target.value)}
            />
          </div>
          
          <div className="flex space-x-2">
            <button
              onClick={() => setActiveTab('all')}
              className={`px-3 py-2 text-sm font-medium rounded-md ${
                activeTab === 'all'
                  ? 'bg-primary text-white'
                  : 'bg-gray-100 dark:bg-gray-800 text-gray-700 dark:text-gray-300 hover:bg-gray-200 dark:hover:bg-gray-700'
              }`}
            >
              All
            </button>
            <button
              onClick={() => setActiveTab('cryopreservation')}
              className={`px-3 py-2 text-sm font-medium rounded-md ${
                activeTab === 'cryopreservation'
                  ? 'bg-primary text-white'
                  : 'bg-gray-100 dark:bg-gray-800 text-gray-700 dark:text-gray-300 hover:bg-gray-200 dark:hover:bg-gray-700'
              }`}
            >
              Cryopreservation
            </button>
            <button
              onClick={() => setActiveTab('vitrification')}
              className={`px-3 py-2 text-sm font-medium rounded-md ${
                activeTab === 'vitrification'
                  ? 'bg-primary text-white'
                  : 'bg-gray-100 dark:bg-gray-800 text-gray-700 dark:text-gray-300 hover:bg-gray-200 dark:hover:bg-gray-700'
              }`}
            >
              Vitrification
            </button>
          </div>
        </div>
      </div>
      
      {/* Protocols List */}
      {filteredProtocols.length === 0 ? (
        <div className="text-center py-10">
          <p className="text-gray-500 dark:text-gray-400">No protocols found matching your criteria</p>
        </div>
      ) : (
        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
          {filteredProtocols.map(protocol => (
            <Link href={`/protocols/${protocol.id}`} key={protocol.id} className="group">
              <div className="h-full bg-white dark:bg-gray-800 border border-gray-200 dark:border-gray-700 rounded-lg shadow-sm hover:shadow-md transition-shadow">
                <div className="p-5">
                  <div className="flex justify-between items-start mb-3">
                    <div>
                      <span className="inline-block px-2 py-1 text-xs font-medium rounded-full bg-blue-100 text-blue-800 dark:bg-blue-900/30 dark:text-blue-300">
                        {protocol.category}
                      </span>
                    </div>
                    <span className="text-xs text-gray-500 dark:text-gray-400">
                      {protocol.createdAt.toLocaleDateString()}
                    </span>
                  </div>
                  
                  <h3 className="text-lg font-semibold mb-2 group-hover:text-primary transition-colors">
                    {protocol.name}
                  </h3>
                  
                  <p className="text-sm text-gray-600 dark:text-gray-300 line-clamp-2 mb-3">
                    {protocol.description}
                  </p>
                  
                  <div className="flex items-center justify-between text-sm">
                    <div className="flex items-center">
                      <svg className="h-4 w-4 text-gray-400 mr-1" xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 21V5a2 2 0 00-2-2H7a2 2 0 00-2 2v16m14 0h2m-2 0h-5m-9 0H3m2 0h5M9 7h1m-1 4h1m4-4h1m-1 4h1m-5 10v-5a1 1 0 011-1h2a1 1 0 011 1v5m-4 0h4" />
                      </svg>
                      <span>{protocol.steps?.length || 0} steps</span>
                    </div>
                    
                    <span className="text-gray-500 dark:text-gray-400">{protocol.author}</span>
                  </div>
                </div>
              </div>
            </Link>
          ))}
        </div>
      )}
    </div>
  );
}
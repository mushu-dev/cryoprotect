import React from 'react';
import Link from 'next/link';
import Head from 'next/head';

export default function ProtocolsPage() {
  const [activeTab, setActiveTab] = React.useState('all');
  
  // Simplified mock data for demonstration
  const myProtocols = [
    {
      id: '1',
      name: 'Standard Cell Freezing',
      description: 'Basic protocol for freezing mammalian cells using DMSO',
      version: '1.2.0',
      created_at: '2023-09-01',
      updated_at: '2023-10-15',
      step_count: 8,
      duration: '2 hours',
      is_template: false,
      author: { name: 'John Doe', id: 'user1' }
    },
    {
      id: '2',
      name: 'Vitrification Protocol',
      description: 'Advanced vitrification method with controlled cooling rates',
      version: '2.0.1',
      created_at: '2023-08-12',
      updated_at: '2023-10-18',
      step_count: 12,
      duration: '3.5 hours',
      is_template: false,
      author: { name: 'John Doe', id: 'user1' }
    }
  ];
  
  const templateProtocols = [
    {
      id: '3',
      name: 'Standard Slow Freezing',
      description: 'Industry standard slow freezing protocol for various cell types',
      version: '3.1.0',
      created_at: '2023-06-15',
      updated_at: '2023-09-20',
      step_count: 10,
      duration: '2.5 hours',
      is_template: true,
      author: { name: 'System', id: 'system' }
    },
    {
      id: '4',
      name: 'Rapid Cooling Method',
      description: 'Quick cooling protocol for small sample volumes',
      version: '1.5.2',
      created_at: '2023-07-22',
      updated_at: '2023-08-30',
      step_count: 6,
      duration: '1 hour',
      is_template: true,
      author: { name: 'System', id: 'system' }
    },
    {
      id: '5',
      name: 'Plant Tissue Cryopreservation',
      description: 'Specialized protocol for preserving plant tissue samples',
      version: '2.2.0',
      created_at: '2023-05-18',
      updated_at: '2023-09-10',
      step_count: 14,
      duration: '4 hours',
      is_template: true,
      author: { name: 'Jane Smith', id: 'user2' }
    }
  ];

  return (
    <>
      <Head>
        <title>Protocols | CryoProtect</title>
        <meta name="description" content="Browse and manage standardized protocols for cryopreservation experiments" />
      </Head>
    
      <div className="container mx-auto px-4 py-8">
        <div className="flex flex-col md:flex-row justify-between items-start md:items-center mb-8">
          <div>
            <h1 className="text-3xl font-bold mb-2">Protocols</h1>
            <p className="text-muted-foreground">
              Standardized protocols for reproducible cryopreservation procedures
            </p>
          </div>
          
          <div className="mt-4 md:mt-0">
            <Link href="/protocols/create">
              <a className="inline-flex items-center justify-center rounded-md bg-primary px-4 py-2 text-sm font-medium text-primary-foreground shadow transition-colors hover:bg-primary/90">
                Create New Protocol
              </a>
            </Link>
          </div>
        </div>
      
        <div className="mb-8">
          <div className="max-w-md mx-auto mb-8">
            <input 
              type="search" 
              placeholder="Search protocols..." 
              className="w-full px-4 py-2 border rounded-md"
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
          
          {activeTab === 'all' && (
            <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
              {[...myProtocols, ...templateProtocols].map(protocol => (
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
                      <div className="mr-4">{protocol.step_count} steps</div>
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
                          <a className="text-primary text-sm hover:underline">View Details</a>
                        </Link>
                      )}
                    </div>
                  </div>
                </div>
              ))}
            </div>
          )}
          
          {activeTab === 'my' && (
            <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
              {myProtocols.map(protocol => (
                <div key={protocol.id} className="rounded-lg border bg-card shadow-sm overflow-hidden">
                  <div className="p-6">
                    <div className="flex flex-wrap items-center gap-2 mb-2">
                      <h3 className="text-lg font-semibold">{protocol.name}</h3>
                      <span className="inline-flex items-center rounded-full bg-muted px-2.5 py-0.5 text-xs font-medium">
                        v{protocol.version}
                      </span>
                    </div>
                    <p className="text-muted-foreground mb-4 text-sm">{protocol.description}</p>
                    <div className="flex text-sm text-muted-foreground mb-4">
                      <div className="mr-4">{protocol.step_count} steps</div>
                      <div>{protocol.duration}</div>
                    </div>
                    <div className="flex justify-between items-center">
                      <div className="text-xs text-muted-foreground">
                        Updated {new Date(protocol.updated_at).toLocaleDateString()}
                      </div>
                      <Link href={`/protocols/${protocol.id}`}>
                        <a className="text-primary text-sm hover:underline">View Details</a>
                      </Link>
                    </div>
                  </div>
                </div>
              ))}
            </div>
          )}
          
          {activeTab === 'templates' && (
            <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
              {templateProtocols.map(protocol => (
                <div key={protocol.id} className="rounded-lg border bg-card shadow-sm overflow-hidden">
                  <div className="p-6">
                    <div className="flex flex-wrap items-center gap-2 mb-2">
                      <h3 className="text-lg font-semibold">{protocol.name}</h3>
                      <span className="inline-flex items-center rounded-full bg-muted px-2.5 py-0.5 text-xs font-medium">
                        v{protocol.version}
                      </span>
                      <span className="inline-flex items-center rounded-full border px-2.5 py-0.5 text-xs font-medium">
                        Template
                      </span>
                    </div>
                    <p className="text-muted-foreground mb-4 text-sm">{protocol.description}</p>
                    <div className="flex text-sm text-muted-foreground mb-4">
                      <div className="mr-4">{protocol.step_count} steps</div>
                      <div>{protocol.duration}</div>
                    </div>
                    <div className="flex justify-between items-center">
                      <div className="text-xs text-muted-foreground">
                        By {protocol.author.name}
                      </div>
                      <button className="text-primary text-sm hover:underline">
                        Use Template
                      </button>
                    </div>
                  </div>
                </div>
              ))}
            </div>
          )}
        </div>
      </div>
    </>
  );
}
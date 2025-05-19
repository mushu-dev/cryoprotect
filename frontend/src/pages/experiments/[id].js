import React from 'react';
import Link from 'next/link';
import Head from 'next/head';
import { useRouter } from 'next/router';

export default function ExperimentDetailPage() {
  const router = useRouter();
  const { id } = router.query;
  
  // This would typically come from an API call
  const experiment = {
    id: id,
    title: `Experiment #${id}`,
    status: id % 2 === 0 ? 'Completed' : 'In Progress',
    description: 'Comprehensive analysis of DMSO and trehalose combinations for mammalian cell cryopreservation.',
    date: `May ${Number(id) + 10}, 2025`,
    createdBy: 'Dr. Jane Smith',
    protocol: {
      id: '1',
      name: 'Standard Cell Freezing Protocol',
      version: '2.1.0'
    },
    cryoprotectants: [
      { id: '1', name: 'DMSO', concentration: '10%' },
      { id: '2', name: 'Trehalose', concentration: '100mM' }
    ],
    cellType: 'Human dermal fibroblasts',
    freezingRate: '-1°C/min controlled rate',
    storageTemperature: '-196°C (liquid nitrogen)',
    thawingMethod: 'Rapid thawing in 37°C water bath',
    results: {
      viability: '87%',
      recovery: '92%',
      functionality: '85%',
      notes: 'Cells exhibited normal morphology and growth characteristics following thawing.'
    },
    attachments: [
      { id: '1', name: 'Raw Data.xlsx', type: 'spreadsheet', size: '1.2 MB' },
      { id: '2', name: 'Microscopy Images.zip', type: 'archive', size: '15.8 MB' },
      { id: '3', name: 'Analysis Report.pdf', type: 'document', size: '3.5 MB' }
    ]
  };

  if (router.isFallback || !id) {
    return <div className="container mx-auto px-4 py-8">Loading experiment details...</div>;
  }

  return (
    <>
      <Head>
        <title>{experiment.title} | CryoProtect</title>
        <meta name="description" content={experiment.description} />
      </Head>
    
      <div className="container mx-auto px-4 py-8">
        <div className="mb-6">
          <Link href="/experiments">
            <a className="text-primary hover:underline flex items-center">
              <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" className="mr-2">
                <path d="M19 12H5"/>
                <path d="M12 19l-7-7 7-7"/>
              </svg>
              Back to Experiments
            </a>
          </Link>
        </div>
        
        <div className="flex flex-col md:flex-row justify-between items-start md:items-center mb-6">
          <div>
            <h1 className="text-3xl font-bold mb-2">{experiment.title}</h1>
            <p className="text-muted-foreground">{experiment.description}</p>
          </div>
          <div className="mt-4 md:mt-0 flex flex-wrap gap-2">
            <span className={`inline-flex items-center rounded-full px-2.5 py-0.5 text-xs font-medium ${
              experiment.status === 'Completed' 
                ? 'bg-green-100 text-green-800' 
                : 'bg-blue-100 text-blue-800'
            }`}>
              {experiment.status}
            </span>
            <button className="inline-flex items-center justify-center rounded-md bg-primary px-4 py-2 text-sm font-medium text-primary-foreground shadow transition-colors hover:bg-primary/90">
              Edit Experiment
            </button>
            <button className="inline-flex items-center justify-center rounded-md bg-destructive px-4 py-2 text-sm font-medium text-destructive-foreground shadow transition-colors hover:bg-destructive/90">
              Delete
            </button>
          </div>
        </div>
        
        <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
          <div className="lg:col-span-2 space-y-6">
            <div className="bg-card rounded-lg border p-6 shadow-sm">
              <h2 className="text-xl font-semibold mb-4">Experiment Overview</h2>
              <dl className="grid grid-cols-1 md:grid-cols-2 gap-x-4 gap-y-6">
                <div>
                  <dt className="text-sm font-medium text-muted-foreground">Date</dt>
                  <dd className="mt-1">{experiment.date}</dd>
                </div>
                <div>
                  <dt className="text-sm font-medium text-muted-foreground">Created By</dt>
                  <dd className="mt-1">{experiment.createdBy}</dd>
                </div>
                <div>
                  <dt className="text-sm font-medium text-muted-foreground">Protocol</dt>
                  <dd className="mt-1">
                    <Link href={`/protocols/${experiment.protocol.id}`}>
                      <a className="text-primary hover:underline">
                        {experiment.protocol.name} (v{experiment.protocol.version})
                      </a>
                    </Link>
                  </dd>
                </div>
                <div>
                  <dt className="text-sm font-medium text-muted-foreground">Cell Type</dt>
                  <dd className="mt-1">{experiment.cellType}</dd>
                </div>
              </dl>
            </div>
            
            <div className="bg-card rounded-lg border p-6 shadow-sm">
              <h2 className="text-xl font-semibold mb-4">Cryopreservation Conditions</h2>
              <dl className="grid grid-cols-1 md:grid-cols-2 gap-x-4 gap-y-6">
                <div>
                  <dt className="text-sm font-medium text-muted-foreground">Cryoprotectants</dt>
                  <dd className="mt-1">
                    <ul className="list-disc list-inside">
                      {experiment.cryoprotectants.map(cp => (
                        <li key={cp.id}>
                          <Link href={`/molecules/${cp.id}`}>
                            <a className="text-primary hover:underline">{cp.name}</a>
                          </Link> ({cp.concentration})
                        </li>
                      ))}
                    </ul>
                  </dd>
                </div>
                <div>
                  <dt className="text-sm font-medium text-muted-foreground">Freezing Rate</dt>
                  <dd className="mt-1">{experiment.freezingRate}</dd>
                </div>
                <div>
                  <dt className="text-sm font-medium text-muted-foreground">Storage Temperature</dt>
                  <dd className="mt-1">{experiment.storageTemperature}</dd>
                </div>
                <div>
                  <dt className="text-sm font-medium text-muted-foreground">Thawing Method</dt>
                  <dd className="mt-1">{experiment.thawingMethod}</dd>
                </div>
              </dl>
            </div>
            
            <div className="bg-card rounded-lg border p-6 shadow-sm">
              <h2 className="text-xl font-semibold mb-4">Results</h2>
              <dl className="grid grid-cols-1 md:grid-cols-3 gap-x-4 gap-y-6 mb-4">
                <div>
                  <dt className="text-sm font-medium text-muted-foreground">Viability</dt>
                  <dd className="mt-1 text-lg font-medium">{experiment.results.viability}</dd>
                </div>
                <div>
                  <dt className="text-sm font-medium text-muted-foreground">Recovery</dt>
                  <dd className="mt-1 text-lg font-medium">{experiment.results.recovery}</dd>
                </div>
                <div>
                  <dt className="text-sm font-medium text-muted-foreground">Functionality</dt>
                  <dd className="mt-1 text-lg font-medium">{experiment.results.functionality}</dd>
                </div>
              </dl>
              <div>
                <h3 className="text-sm font-medium text-muted-foreground mb-2">Notes</h3>
                <p>{experiment.results.notes}</p>
              </div>
            </div>
            
            <div className="bg-card rounded-lg border p-6 shadow-sm">
              <h2 className="text-xl font-semibold mb-4">Attachments</h2>
              <ul className="divide-y">
                {experiment.attachments.map(file => (
                  <li key={file.id} className="py-3 flex justify-between items-center">
                    <div className="flex items-center">
                      <svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" className="mr-3 text-muted-foreground">
                        <path d="M14 2H6a2 2 0 0 0-2 2v16a2 2 0 0 0 2 2h12a2 2 0 0 0 2-2V8z"/>
                        <polyline points="14 2 14 8 20 8"/>
                        <line x1="16" y1="13" x2="8" y2="13"/>
                        <line x1="16" y1="17" x2="8" y2="17"/>
                        <polyline points="10 9 9 9 8 9"/>
                      </svg>
                      <div>
                        <p className="font-medium">{file.name}</p>
                        <p className="text-sm text-muted-foreground">{file.type} • {file.size}</p>
                      </div>
                    </div>
                    <button className="inline-flex items-center justify-center rounded-md bg-secondary px-3 py-1.5 text-sm font-medium text-secondary-foreground shadow-sm hover:bg-secondary/80">
                      Download
                    </button>
                  </li>
                ))}
              </ul>
            </div>
          </div>
          
          <div className="space-y-6">
            <div className="bg-card rounded-lg border p-6 shadow-sm">
              <h2 className="text-xl font-semibold mb-4">Activity Log</h2>
              <div className="space-y-4">
                <div className="border-l-2 border-muted pl-4 pb-4">
                  <div className="text-sm font-medium">{experiment.date}</div>
                  <div className="mt-1">Experiment completed</div>
                </div>
                <div className="border-l-2 border-muted pl-4 pb-4">
                  <div className="text-sm font-medium">{experiment.date.replace(`${Number(id) + 10}`, `${Number(id) + 9}`)}</div>
                  <div className="mt-1">Thawing and recovery analysis performed</div>
                </div>
                <div className="border-l-2 border-muted pl-4 pb-4">
                  <div className="text-sm font-medium">{experiment.date.replace(`${Number(id) + 10}`, `${Number(id) + 5}`)}</div>
                  <div className="mt-1">Samples transferred to long-term storage</div>
                </div>
                <div className="border-l-2 border-muted pl-4 pb-4">
                  <div className="text-sm font-medium">{experiment.date.replace(`${Number(id) + 10}`, `${Number(id) + 3}`)}</div>
                  <div className="mt-1">Freezing protocol initiated</div>
                </div>
                <div className="border-l-2 border-muted pl-4">
                  <div className="text-sm font-medium">{experiment.date.replace(`${Number(id) + 10}`, `${Number(id) + 1}`)}</div>
                  <div className="mt-1">Experiment created</div>
                </div>
              </div>
            </div>
            
            <div className="bg-card rounded-lg border p-6 shadow-sm">
              <h2 className="text-xl font-semibold mb-4">Related Experiments</h2>
              <div className="space-y-3">
                <div className="rounded border p-3">
                  <h3 className="font-medium mb-1">
                    <Link href={`/experiments/${Number(id) + 1}`}>
                      <a className="text-primary hover:underline">Experiment #{Number(id) + 1}</a>
                    </Link>
                  </h3>
                  <p className="text-sm text-muted-foreground">Similar protocol with different cell type</p>
                </div>
                <div className="rounded border p-3">
                  <h3 className="font-medium mb-1">
                    <Link href={`/experiments/${Number(id) - 1 > 0 ? Number(id) - 1 : Number(id) + 2}`}>
                      <a className="text-primary hover:underline">Experiment #{Number(id) - 1 > 0 ? Number(id) - 1 : Number(id) + 2}</a>
                    </Link>
                  </h3>
                  <p className="text-sm text-muted-foreground">Different cryoprotectant combination</p>
                </div>
              </div>
            </div>
          </div>
        </div>
      </div>
    </>
  );
}
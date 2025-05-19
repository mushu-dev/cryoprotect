import React from 'react';
import Link from 'next/link';
import Head from 'next/head';
import { useRouter } from 'next/router';

export default function ProtocolDetailPage() {
  const router = useRouter();
  const { id } = router.query;
  
  // This would typically come from an API call
  const protocol = {
    id: id,
    name: id === '1' ? 'Standard Cell Freezing' : 'Vitrification Protocol',
    version: id === '1' ? '1.2.0' : '2.0.1',
    description: id === '1' 
      ? 'Basic protocol for freezing mammalian cells using DMSO as a cryoprotectant with controlled rate cooling.'
      : 'Advanced vitrification method with controlled cooling rates and combination cryoprotectant solution.',
    author: { name: 'Dr. John Doe', id: 'user1' },
    created_at: '2023-09-01',
    updated_at: '2023-10-15',
    duration: id === '1' ? '2 hours' : '3.5 hours',
    is_template: false,
    cell_types: ['Human embryonic stem cells', 'Primary fibroblasts', 'CHO cells'],
    equipment: [
      'Controlled-rate freezer',
      'Cryovials',
      'Liquid nitrogen storage tank',
      'Water bath (37°C)',
      'Cell culture materials'
    ],
    cryoprotectants: [
      { id: '1', name: 'DMSO', concentration: '10%' },
      { id: '2', name: 'Fetal bovine serum', concentration: '20%' }
    ],
    freezing_rate: '-1°C/min',
    storage_temperature: '-196°C (liquid nitrogen)',
    thawing_method: 'Rapid thawing in 37°C water bath',
    steps: [
      {
        order: 1,
        title: 'Prepare cryopreservation medium',
        description: 'Mix cell culture medium with DMSO to achieve 10% final concentration and FBS to 20% final concentration. Keep on ice.',
        duration: '15 minutes',
        is_critical: true
      },
      {
        order: 2,
        title: 'Prepare cells',
        description: 'Harvest cells in log phase growth (70-80% confluency). Count cells and check viability.',
        duration: '20 minutes',
        is_critical: true
      },
      {
        order: 3,
        title: 'Centrifuge cell suspension',
        description: 'Centrifuge at 200 x g for 5 minutes at room temperature.',
        duration: '10 minutes',
        is_critical: false
      },
      {
        order: 4,
        title: 'Resuspend cells in cryopreservation medium',
        description: 'Gently resuspend cell pellet in cold cryopreservation medium to achieve 1-5 x 10^6 cells/mL.',
        duration: '10 minutes',
        is_critical: true
      },
      {
        order: 5,
        title: 'Aliquot into cryovials',
        description: 'Transfer 1 mL of cell suspension to each pre-labeled cryovial.',
        duration: '15 minutes',
        is_critical: false
      },
      {
        order: 6,
        title: 'Controlled-rate freezing',
        description: 'Place cryovials in controlled-rate freezer and run program: hold at 4°C for 10 minutes, then cool at -1°C/min to -80°C.',
        duration: '90 minutes',
        is_critical: true
      },
      {
        order: 7,
        title: 'Transfer to long-term storage',
        description: 'Quickly transfer frozen vials to liquid nitrogen storage tank.',
        duration: '10 minutes',
        is_critical: true
      }
    ],
    notes: 'For optimal results, cells should be in log phase of growth and have viability >90% before freezing. DMSO is toxic to cells at room temperature, so work quickly when adding cells to freezing medium.',
    references: [
      { title: 'Optimization of cryopreservation procedures for mammalian cells', authors: 'Smith et al.', journal: 'Cryobiology', year: 2020 },
      { title: 'Best practices for cell banking', authors: 'Johnson & Williams', journal: 'Nature Methods', year: 2019 }
    ],
    used_in_experiments: [
      { id: '1', title: 'DMSO concentration optimization', date: '2023-11-10' },
      { id: '3', title: 'Fibroblast freezing comparison', date: '2023-10-20' }
    ]
  };

  if (router.isFallback || !id) {
    return <div className="container mx-auto px-4 py-8">Loading protocol details...</div>;
  }

  return (
    <>
      <Head>
        <title>{protocol.name} | CryoProtect</title>
        <meta name="description" content={protocol.description} />
      </Head>
    
      <div className="container mx-auto px-4 py-8">
        <div className="mb-6">
          <Link href="/protocols">
            <a className="text-primary hover:underline flex items-center">
              <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" className="mr-2">
                <path d="M19 12H5"/>
                <path d="M12 19l-7-7 7-7"/>
              </svg>
              Back to Protocols
            </a>
          </Link>
        </div>
        
        <div className="flex flex-col md:flex-row justify-between items-start md:items-center mb-6">
          <div>
            <div className="flex flex-wrap items-center gap-2 mb-2">
              <h1 className="text-3xl font-bold">{protocol.name}</h1>
              <span className="inline-flex items-center rounded-full bg-muted px-2.5 py-0.5 text-xs font-medium">
                v{protocol.version}
              </span>
              {protocol.is_template && (
                <span className="inline-flex items-center rounded-full border px-2.5 py-0.5 text-xs font-medium">
                  Template
                </span>
              )}
            </div>
            <p className="text-muted-foreground">{protocol.description}</p>
          </div>
          <div className="mt-4 md:mt-0 flex flex-wrap gap-2">
            <button className="inline-flex items-center justify-center rounded-md bg-primary px-4 py-2 text-sm font-medium text-primary-foreground shadow transition-colors hover:bg-primary/90">
              Use Protocol
            </button>
            <button className="inline-flex items-center justify-center rounded-md bg-secondary px-4 py-2 text-sm font-medium text-secondary-foreground shadow transition-colors hover:bg-secondary/80">
              Duplicate
            </button>
            <button className="inline-flex items-center justify-center rounded-md bg-secondary px-4 py-2 text-sm font-medium text-secondary-foreground shadow transition-colors hover:bg-secondary/80">
              Export
            </button>
          </div>
        </div>
        
        <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
          <div className="lg:col-span-2 space-y-6">
            <div className="bg-card rounded-lg border p-6 shadow-sm">
              <h2 className="text-xl font-semibold mb-4">Protocol Steps</h2>
              <ol className="space-y-5">
                {protocol.steps.map((step) => (
                  <li key={step.order} className="relative pl-8 pb-5 border-l border-muted last:border-l-0">
                    <div className="absolute left-0 top-0 flex items-center justify-center h-6 w-6 rounded-full bg-primary text-white text-sm font-medium">
                      {step.order}
                    </div>
                    <div className="ml-2">
                      <h3 className="font-medium">
                        {step.title}
                        {step.is_critical && (
                          <span className="ml-2 inline-flex items-center rounded-full bg-destructive/10 px-2 py-0.5 text-xs font-medium text-destructive">
                            Critical Step
                          </span>
                        )}
                      </h3>
                      <p className="mt-1 text-muted-foreground text-sm">{step.description}</p>
                      <p className="mt-1 text-xs font-medium">Duration: {step.duration}</p>
                    </div>
                  </li>
                ))}
              </ol>
            </div>
            
            <div className="bg-card rounded-lg border p-6 shadow-sm">
              <h2 className="text-xl font-semibold mb-4">Materials and Equipment</h2>
              <div className="mb-6">
                <h3 className="text-base font-medium mb-2">Cryoprotectants</h3>
                <ul className="list-disc list-inside space-y-1">
                  {protocol.cryoprotectants.map((cp, index) => (
                    <li key={index}>
                      <Link href={`/molecules/${cp.id}`}>
                        <a className="text-primary hover:underline">{cp.name}</a>
                      </Link> ({cp.concentration})
                    </li>
                  ))}
                </ul>
              </div>
              <div>
                <h3 className="text-base font-medium mb-2">Equipment</h3>
                <ul className="list-disc list-inside space-y-1">
                  {protocol.equipment.map((item, index) => (
                    <li key={index}>{item}</li>
                  ))}
                </ul>
              </div>
            </div>
            
            <div className="bg-card rounded-lg border p-6 shadow-sm">
              <h2 className="text-xl font-semibold mb-4">Parameters</h2>
              <dl className="grid grid-cols-1 md:grid-cols-2 gap-x-4 gap-y-6">
                <div>
                  <dt className="text-sm font-medium text-muted-foreground">Compatible Cell Types</dt>
                  <dd className="mt-1">
                    <ul className="list-disc list-inside">
                      {protocol.cell_types.map((cell, index) => (
                        <li key={index}>{cell}</li>
                      ))}
                    </ul>
                  </dd>
                </div>
                <div>
                  <dt className="text-sm font-medium text-muted-foreground">Freezing Rate</dt>
                  <dd className="mt-1">{protocol.freezing_rate}</dd>
                </div>
                <div>
                  <dt className="text-sm font-medium text-muted-foreground">Storage Temperature</dt>
                  <dd className="mt-1">{protocol.storage_temperature}</dd>
                </div>
                <div>
                  <dt className="text-sm font-medium text-muted-foreground">Thawing Method</dt>
                  <dd className="mt-1">{protocol.thawing_method}</dd>
                </div>
                <div>
                  <dt className="text-sm font-medium text-muted-foreground">Total Duration</dt>
                  <dd className="mt-1">{protocol.duration}</dd>
                </div>
              </dl>
            </div>
            
            <div className="bg-card rounded-lg border p-6 shadow-sm">
              <h2 className="text-xl font-semibold mb-4">Notes</h2>
              <div className="text-sm prose max-w-none">
                <p>{protocol.notes}</p>
              </div>
            </div>
            
            {protocol.references.length > 0 && (
              <div className="bg-card rounded-lg border p-6 shadow-sm">
                <h2 className="text-xl font-semibold mb-4">References</h2>
                <ul className="space-y-3">
                  {protocol.references.map((ref, index) => (
                    <li key={index} className="text-sm">
                      <span className="font-medium">{ref.title}</span>
                      <div className="text-muted-foreground">
                        {ref.authors}. {ref.journal}, {ref.year}
                      </div>
                    </li>
                  ))}
                </ul>
              </div>
            )}
          </div>
          
          <div className="space-y-6">
            <div className="bg-card rounded-lg border p-6 shadow-sm">
              <h2 className="text-lg font-semibold mb-4">Protocol Information</h2>
              <dl className="divide-y">
                <div className="py-3 flex justify-between">
                  <dt className="text-sm font-medium text-muted-foreground">Author</dt>
                  <dd className="text-sm">{protocol.author.name}</dd>
                </div>
                <div className="py-3 flex justify-between">
                  <dt className="text-sm font-medium text-muted-foreground">Created</dt>
                  <dd className="text-sm">{new Date(protocol.created_at).toLocaleDateString()}</dd>
                </div>
                <div className="py-3 flex justify-between">
                  <dt className="text-sm font-medium text-muted-foreground">Last Updated</dt>
                  <dd className="text-sm">{new Date(protocol.updated_at).toLocaleDateString()}</dd>
                </div>
                <div className="py-3 flex justify-between">
                  <dt className="text-sm font-medium text-muted-foreground">Version</dt>
                  <dd className="text-sm">{protocol.version}</dd>
                </div>
                <div className="py-3 flex justify-between">
                  <dt className="text-sm font-medium text-muted-foreground">Steps</dt>
                  <dd className="text-sm">{protocol.steps.length}</dd>
                </div>
              </dl>
            </div>
            
            <div className="bg-card rounded-lg border p-6 shadow-sm">
              <h2 className="text-lg font-semibold mb-4">Used in Experiments</h2>
              {protocol.used_in_experiments.length > 0 ? (
                <ul className="space-y-3">
                  {protocol.used_in_experiments.map((exp) => (
                    <li key={exp.id} className="border-b pb-3 last:border-0 last:pb-0">
                      <Link href={`/experiments/${exp.id}`}>
                        <a className="text-primary hover:underline font-medium">{exp.title}</a>
                      </Link>
                      <div className="text-sm text-muted-foreground">
                        {new Date(exp.date).toLocaleDateString()}
                      </div>
                    </li>
                  ))}
                </ul>
              ) : (
                <p className="text-sm text-muted-foreground">
                  This protocol hasn't been used in experiments yet.
                </p>
              )}
            </div>
            
            <div className="bg-card rounded-lg border p-6 shadow-sm">
              <h2 className="text-lg font-semibold mb-4">Quick Actions</h2>
              <div className="space-y-2">
                <Link href={`/protocols/${id}/edit`}>
                  <a className="w-full inline-flex items-center justify-center rounded-md text-sm font-medium transition-colors focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-ring focus-visible:ring-offset-2 disabled:opacity-50 disabled:pointer-events-none ring-offset-background bg-secondary text-secondary-foreground hover:bg-secondary/80 h-10 py-2 px-4">
                    Edit Protocol
                  </a>
                </Link>
                <button className="w-full inline-flex items-center justify-center rounded-md text-sm font-medium transition-colors focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-ring focus-visible:ring-offset-2 disabled:opacity-50 disabled:pointer-events-none ring-offset-background bg-secondary text-secondary-foreground hover:bg-secondary/80 h-10 py-2 px-4">
                  Download PDF
                </button>
                <button className="w-full inline-flex items-center justify-center rounded-md text-sm font-medium transition-colors focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-ring focus-visible:ring-offset-2 disabled:opacity-50 disabled:pointer-events-none ring-offset-background border border-input bg-background hover:bg-accent hover:text-accent-foreground h-10 py-2 px-4">
                  Share Protocol
                </button>
              </div>
            </div>
          </div>
        </div>
      </div>
    </>
  );
}
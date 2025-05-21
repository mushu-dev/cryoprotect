import React, { useState, useEffect } from 'react';
import Link from 'next/link';
import Head from 'next/head';
import { useRouter } from 'next/router';
import { ProtocolVersionComparison } from '../../features/protocols/components/protocol-version-comparison';
import { useProtocol, useProtocolVersions } from '../../features/protocols/hooks/use-protocols';
import { 
  ArrowLeft, 
  ChevronDown, 
  ChevronUp, 
  Clock, 
  Download, 
  FileClock, 
  FileEdit, 
  FileText, 
  History, 
  Share, 
  Users,
  FlaskConical
} from 'lucide-react';

export default function ProtocolDetailPage() {
  const router = useRouter();
  const { id } = router.query;
  const [showVersionHistory, setShowVersionHistory] = useState(false);
  const [selectedVersionId, setSelectedVersionId] = useState(null);
  const [comparisonProtocol, setComparisonProtocol] = useState(null);
  const [showComparison, setShowComparison] = useState(false);
  const [activeTab, setActiveTab] = useState('steps');
  
  // Fetch the current protocol and its versions
  const { protocol, loading: protocolLoading } = useProtocol(id);
  const { versions, loading: versionsLoading, compareVersions } = useProtocolVersions(id);
  
  // If we have a selected version different from the current one, fetch it for comparison
  useEffect(() => {
    if (selectedVersionId && selectedVersionId !== id) {
      // In a real app, we would fetch the selected version protocol here
      // For now, we'll use dummy data
      setComparisonProtocol(generateDummyProtocol(selectedVersionId));
    } else {
      setComparisonProtocol(null);
      setShowComparison(false);
    }
  }, [selectedVersionId, id]);
  
  const handleVersionSelect = (versionId) => {
    if (versionId === selectedVersionId) {
      setSelectedVersionId(null);
      setComparisonProtocol(null);
      setShowComparison(false);
    } else {
      setSelectedVersionId(versionId);
    }
  };
  
  const handleCompare = () => {
    if (selectedVersionId && protocol) {
      setShowComparison(true);
    }
  };
  
  const handleViewVersion = (selectedProtocol) => {
    if (selectedProtocol.id !== id) {
      router.push(`/protocols/${selectedProtocol.id}`);
    }
  };
  
  if (router.isFallback || !id || protocolLoading) {
    return (
      <div className="container mx-auto px-4 py-8">
        <div className="flex justify-center items-center h-64">
          <div className="text-center">
            <div className="animate-spin rounded-full h-12 w-12 border-t-2 border-b-2 border-primary mx-auto mb-4"></div>
            <p className="text-muted-foreground">Loading protocol details...</p>
          </div>
        </div>
      </div>
    );
  }
  
  // This would be replaced with real data in a production application
  const protocolData = protocol || {
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
        id: 's1',
        order: 1,
        title: 'Prepare cryopreservation medium',
        description: 'Mix cell culture medium with DMSO to achieve 10% final concentration and FBS to 20% final concentration. Keep on ice.',
        duration: '15 minutes',
        is_critical: true
      },
      {
        id: 's2',
        order: 2,
        title: 'Prepare cells',
        description: 'Harvest cells in log phase growth (70-80% confluency). Count cells and check viability.',
        duration: '20 minutes',
        is_critical: true
      },
      {
        id: 's3',
        order: 3,
        title: 'Centrifuge cell suspension',
        description: 'Centrifuge at 200 x g for 5 minutes at room temperature.',
        duration: '10 minutes',
        is_critical: false
      },
      {
        id: 's4',
        order: 4,
        title: 'Resuspend cells in cryopreservation medium',
        description: 'Gently resuspend cell pellet in cold cryopreservation medium to achieve 1-5 x 10^6 cells/mL.',
        duration: '10 minutes',
        is_critical: true
      },
      {
        id: 's5',
        order: 5,
        title: 'Aliquot into cryovials',
        description: 'Transfer 1 mL of cell suspension to each pre-labeled cryovial.',
        duration: '15 minutes',
        is_critical: false
      },
      {
        id: 's6',
        order: 6,
        title: 'Controlled-rate freezing',
        description: 'Place cryovials in controlled-rate freezer and run program: hold at 4°C for 10 minutes, then cool at -1°C/min to -80°C.',
        duration: '90 minutes',
        is_critical: true
      },
      {
        id: 's7',
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
  
  // Generate dummy version history for the protocol
  const versionHistoryData = versions?.length > 0 ? versions : [
    {
      id: id,
      version: protocolData.version,
      created_at: protocolData.created_at,
      created_by: protocolData.author.name,
      description: 'Current version'
    },
    {
      id: id === '1' ? '3' : '4',  // Different IDs for previous versions
      version: id === '1' ? '1.1.0' : '1.9.0',
      created_at: '2023-08-15',
      created_by: protocolData.author.name,
      description: 'Added improved cooling rates'
    },
    {
      id: id === '1' ? '4' : '5',
      version: id === '1' ? '1.0.0' : '1.8.0',
      created_at: '2023-07-01',
      created_by: protocolData.author.name,
      description: 'Initial version'
    }
  ];
  
  return (
    <>
      <Head>
        <title>{protocolData.name} | CryoProtect</title>
        <meta name="description" content={protocolData.description} />
      </Head>
    
      <div className="container mx-auto px-4 py-8">
        {/* Breadcrumb Navigation */}
        <nav className="flex mb-6" aria-label="Breadcrumb">
          <ol className="inline-flex items-center space-x-1 md:space-x-3">
            <li className="inline-flex items-center">
              <Link href="/">
                <a className="inline-flex items-center text-sm font-medium text-muted-foreground hover:text-foreground">
                  <svg className="w-3 h-3 mr-2.5" aria-hidden="true" xmlns="http://www.w3.org/2000/svg" fill="currentColor" viewBox="0 0 20 20">
                    <path d="m19.707 9.293-9-9a1 1 0 0 0-1.414 0l-9 9a1 1 0 0 0 0 1.414l9 9a1 1 0 0 0 1.414 0l9-9a1 1 0 0 0 0-1.414Z"/>
                  </svg>
                  Home
                </a>
              </Link>
            </li>
            <li>
              <div className="flex items-center">
                <svg className="w-3 h-3 text-muted-foreground mx-1" aria-hidden="true" xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 6 10">
                  <path stroke="currentColor" strokeLinecap="round" strokeLinejoin="round" strokeWidth="2" d="m1 9 4-4-4-4"/>
                </svg>
                <Link href="/protocols">
                  <a className="ml-1 text-sm font-medium text-muted-foreground hover:text-foreground md:ml-2">Protocols</a>
                </Link>
              </div>
            </li>
            <li aria-current="page">
              <div className="flex items-center">
                <svg className="w-3 h-3 text-muted-foreground mx-1" aria-hidden="true" xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 6 10">
                  <path stroke="currentColor" strokeLinecap="round" strokeLinejoin="round" strokeWidth="2" d="m1 9 4-4-4-4"/>
                </svg>
                <span className="ml-1 text-sm font-medium text-foreground md:ml-2 truncate max-w-[200px]">
                  {protocolData.name}
                </span>
              </div>
            </li>
          </ol>
        </nav>
        
        {/* Protocol Header */}
        <div className="bg-card rounded-lg border shadow-sm p-6 mb-6">
          <div className="flex flex-col md:flex-row justify-between items-start md:items-center gap-4">
            <div>
              <div className="flex flex-wrap items-center gap-2 mb-2">
                <h1 className="text-3xl font-bold">{protocolData.name}</h1>
                <span className="inline-flex items-center rounded-full bg-muted px-2.5 py-0.5 text-xs font-medium">
                  v{protocolData.version}
                </span>
                {protocolData.is_template && (
                  <span className="inline-flex items-center rounded-full border px-2.5 py-0.5 text-xs font-medium">
                    Template
                  </span>
                )}
              </div>
              <p className="text-muted-foreground">{protocolData.description}</p>
              
              <div className="flex flex-wrap items-center gap-x-4 gap-y-1 mt-2 text-sm text-muted-foreground">
                <span className="flex items-center">
                  <Users className="h-4 w-4 mr-1" />
                  Author: {protocolData.author.name}
                </span>
                <span className="flex items-center">
                  <Clock className="h-4 w-4 mr-1" />
                  Duration: {protocolData.duration}
                </span>
                <span className="flex items-center">
                  <FileClock className="h-4 w-4 mr-1" />
                  Updated: {new Date(protocolData.updated_at).toLocaleDateString()}
                </span>
              </div>
            </div>
            
            <div className="flex flex-wrap gap-2 mt-4 md:mt-0">
              <Link href={`/protocols/${id}/experiment-integration`}>
                <a className="inline-flex items-center justify-center rounded-md bg-primary px-4 py-2 text-sm font-medium text-primary-foreground shadow transition-colors hover:bg-primary/90">
                  <FlaskConical className="mr-2 h-4 w-4" />
                  Create Experiment
                </a>
              </Link>
              <button className="inline-flex items-center justify-center rounded-md bg-secondary px-4 py-2 text-sm font-medium text-secondary-foreground shadow transition-colors hover:bg-secondary/80">
                <FileEdit className="mr-2 h-4 w-4" />
                Edit
              </button>
              <button className="inline-flex items-center justify-center rounded-md bg-secondary px-4 py-2 text-sm font-medium text-secondary-foreground shadow transition-colors hover:bg-secondary/80">
                <FileText className="mr-2 h-4 w-4" />
                Duplicate
              </button>
              <button className="inline-flex items-center justify-center rounded-md bg-secondary px-4 py-2 text-sm font-medium text-secondary-foreground shadow transition-colors hover:bg-secondary/80">
                <Download className="mr-2 h-4 w-4" />
                Export
              </button>
            </div>
          </div>
        </div>
        
        {/* Version History Toggle */}
        <div className="mb-6">
          <button 
            onClick={() => setShowVersionHistory(!showVersionHistory)}
            className="flex items-center text-sm font-medium text-primary hover:underline"
          >
            <History className="mr-1 h-4 w-4" />
            {showVersionHistory ? "Hide Version History" : "Show Version History"}
            {showVersionHistory ? <ChevronUp className="ml-1 h-4 w-4" /> : <ChevronDown className="ml-1 h-4 w-4" />}
          </button>
        </div>
        
        {/* Version History Panel */}
        {showVersionHistory && (
          <div className="mb-6 border rounded-lg p-4 bg-card shadow-sm">
            <div className="flex justify-between items-center mb-4">
              <h2 className="text-xl font-semibold">Version History</h2>
              <button
                onClick={handleCompare}
                disabled={!selectedVersionId}
                className={`inline-flex items-center px-3 py-1 rounded-md text-sm font-medium ${
                  selectedVersionId 
                    ? 'bg-primary text-primary-foreground' 
                    : 'bg-muted text-muted-foreground cursor-not-allowed'
                }`}
              >
                Compare with Current
              </button>
            </div>
            
            <div className="overflow-x-auto">
              <table className="w-full">
                <thead>
                  <tr className="border-b">
                    <th className="text-left py-2 px-3 font-medium text-sm"></th>
                    <th className="text-left py-2 px-3 font-medium text-sm">Version</th>
                    <th className="text-left py-2 px-3 font-medium text-sm">Date</th>
                    <th className="text-left py-2 px-3 font-medium text-sm">Author</th>
                    <th className="text-left py-2 px-3 font-medium text-sm">Notes</th>
                    <th className="text-left py-2 px-3 font-medium text-sm">Actions</th>
                  </tr>
                </thead>
                <tbody>
                  {versionHistoryData.map((version) => (
                    <tr key={version.id} className="border-b hover:bg-muted/50">
                      <td className="py-2 px-3">
                        <input 
                          type="radio" 
                          name="versionSelect" 
                          checked={selectedVersionId === version.id}
                          disabled={version.id === id}
                          onChange={() => handleVersionSelect(version.id)}
                          className={version.id === id ? 'opacity-50 cursor-not-allowed' : ''}
                        />
                      </td>
                      <td className="py-2 px-3">
                        <span className={version.id === id ? 'font-semibold' : ''}>
                          v{version.version}
                          {version.id === id && ' (current)'}
                        </span>
                      </td>
                      <td className="py-2 px-3">{new Date(version.created_at).toLocaleDateString()}</td>
                      <td className="py-2 px-3">{version.created_by}</td>
                      <td className="py-2 px-3">{version.description}</td>
                      <td className="py-2 px-3">
                        {version.id !== id && (
                          <Link href={`/protocols/${version.id}`}>
                            <a className="text-primary hover:underline text-sm">View</a>
                          </Link>
                        )}
                      </td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </div>
        )}
        
        {/* Version Comparison */}
        {showComparison && comparisonProtocol && (
          <div className="mb-6">
            <ProtocolVersionComparison 
              protocolA={comparisonProtocol} 
              protocolB={protocolData}
              onViewVersion={handleViewVersion}
            />
          </div>
        )}
        
        {/* Protocol Content */}
        <div className="mb-6">
          <div className="border-b">
            <nav className="-mb-px flex space-x-8">
              <button
                onClick={() => setActiveTab('steps')}
                className={`py-4 px-1 inline-flex items-center border-b-2 font-medium text-sm ${
                  activeTab === 'steps' 
                    ? 'border-primary text-primary' 
                    : 'border-transparent text-muted-foreground hover:text-foreground hover:border-muted'
                }`}
              >
                Steps
              </button>
              <button
                onClick={() => setActiveTab('materials')}
                className={`py-4 px-1 inline-flex items-center border-b-2 font-medium text-sm ${
                  activeTab === 'materials' 
                    ? 'border-primary text-primary' 
                    : 'border-transparent text-muted-foreground hover:text-foreground hover:border-muted'
                }`}
              >
                Materials & Equipment
              </button>
              <button
                onClick={() => setActiveTab('parameters')}
                className={`py-4 px-1 inline-flex items-center border-b-2 font-medium text-sm ${
                  activeTab === 'parameters' 
                    ? 'border-primary text-primary' 
                    : 'border-transparent text-muted-foreground hover:text-foreground hover:border-muted'
                }`}
              >
                Parameters
              </button>
              <button
                onClick={() => setActiveTab('notes')}
                className={`py-4 px-1 inline-flex items-center border-b-2 font-medium text-sm ${
                  activeTab === 'notes' 
                    ? 'border-primary text-primary' 
                    : 'border-transparent text-muted-foreground hover:text-foreground hover:border-muted'
                }`}
              >
                Notes & References
              </button>
            </nav>
          </div>
        </div>
        
        <div className="grid grid-cols-1 lg:grid-cols-4 gap-6">
          {/* Main Content */}
          <div className="lg:col-span-3 space-y-6">
            {/* Steps Tab */}
            {activeTab === 'steps' && (
              <div className="bg-card rounded-lg border p-6 shadow-sm">
                <h2 className="text-xl font-semibold mb-4">Protocol Steps</h2>
                <ol className="space-y-5">
                  {protocolData.steps.map((step) => (
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
            )}
            
            {/* Materials Tab */}
            {activeTab === 'materials' && (
              <div className="bg-card rounded-lg border p-6 shadow-sm">
                <h2 className="text-xl font-semibold mb-4">Materials and Equipment</h2>
                <div className="mb-6">
                  <h3 className="text-base font-medium mb-2">Cryoprotectants</h3>
                  <ul className="list-disc list-inside space-y-1">
                    {protocolData.cryoprotectants.map((cp, index) => (
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
                    {protocolData.equipment.map((item, index) => (
                      <li key={index}>{item}</li>
                    ))}
                  </ul>
                </div>
              </div>
            )}
            
            {/* Parameters Tab */}
            {activeTab === 'parameters' && (
              <div className="bg-card rounded-lg border p-6 shadow-sm">
                <h2 className="text-xl font-semibold mb-4">Parameters</h2>
                <dl className="grid grid-cols-1 md:grid-cols-2 gap-x-4 gap-y-6">
                  <div>
                    <dt className="text-sm font-medium text-muted-foreground">Compatible Cell Types</dt>
                    <dd className="mt-1">
                      <ul className="list-disc list-inside">
                        {protocolData.cell_types.map((cell, index) => (
                          <li key={index}>{cell}</li>
                        ))}
                      </ul>
                    </dd>
                  </div>
                  <div>
                    <dt className="text-sm font-medium text-muted-foreground">Freezing Rate</dt>
                    <dd className="mt-1">{protocolData.freezing_rate}</dd>
                  </div>
                  <div>
                    <dt className="text-sm font-medium text-muted-foreground">Storage Temperature</dt>
                    <dd className="mt-1">{protocolData.storage_temperature}</dd>
                  </div>
                  <div>
                    <dt className="text-sm font-medium text-muted-foreground">Thawing Method</dt>
                    <dd className="mt-1">{protocolData.thawing_method}</dd>
                  </div>
                  <div>
                    <dt className="text-sm font-medium text-muted-foreground">Total Duration</dt>
                    <dd className="mt-1">{protocolData.duration}</dd>
                  </div>
                </dl>
              </div>
            )}
            
            {/* Notes Tab */}
            {activeTab === 'notes' && (
              <>
                <div className="bg-card rounded-lg border p-6 shadow-sm">
                  <h2 className="text-xl font-semibold mb-4">Notes</h2>
                  <div className="text-sm prose max-w-none">
                    <p>{protocolData.notes}</p>
                  </div>
                </div>
                
                {protocolData.references.length > 0 && (
                  <div className="bg-card rounded-lg border p-6 shadow-sm">
                    <h2 className="text-xl font-semibold mb-4">References</h2>
                    <ul className="space-y-3">
                      {protocolData.references.map((ref, index) => (
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
              </>
            )}
          </div>
          
          {/* Sidebar */}
          <div className="space-y-6">
            <div className="bg-card rounded-lg border p-6 shadow-sm sticky top-6">
              <h2 className="text-lg font-semibold mb-4">Quick Actions</h2>
              <div className="space-y-2">
                <Link href={`/protocols/${id}/edit`}>
                  <a className="w-full inline-flex items-center justify-center rounded-md text-sm font-medium transition-colors focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-ring focus-visible:ring-offset-2 disabled:opacity-50 disabled:pointer-events-none ring-offset-background bg-secondary text-secondary-foreground hover:bg-secondary/80 h-10 py-2 px-4">
                    <FileEdit className="mr-2 h-4 w-4" />
                    Edit Protocol
                  </a>
                </Link>
                <button className="w-full inline-flex items-center justify-center rounded-md text-sm font-medium transition-colors focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-ring focus-visible:ring-offset-2 disabled:opacity-50 disabled:pointer-events-none ring-offset-background bg-secondary text-secondary-foreground hover:bg-secondary/80 h-10 py-2 px-4">
                  <Download className="mr-2 h-4 w-4" />
                  Download PDF
                </button>
                <button className="w-full inline-flex items-center justify-center rounded-md text-sm font-medium transition-colors focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-ring focus-visible:ring-offset-2 disabled:opacity-50 disabled:pointer-events-none ring-offset-background border border-input bg-background hover:bg-accent hover:text-accent-foreground h-10 py-2 px-4">
                  <Share className="mr-2 h-4 w-4" />
                  Share Protocol
                </button>
                <Link href={`/protocols/${id}/experiment-integration`}>
                  <a className="w-full inline-flex items-center justify-center rounded-md text-sm font-medium transition-colors focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-ring focus-visible:ring-offset-2 disabled:opacity-50 disabled:pointer-events-none ring-offset-background bg-primary text-primary-foreground hover:bg-primary/90 h-10 py-2 px-4">
                    <FlaskConical className="mr-2 h-4 w-4" />
                    Experiment Integration
                  </a>
                </Link>
                
                <button 
                  onClick={() => setShowVersionHistory(!showVersionHistory)}
                  className="w-full inline-flex items-center justify-center rounded-md text-sm font-medium transition-colors focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-ring focus-visible:ring-offset-2 disabled:opacity-50 disabled:pointer-events-none ring-offset-background border border-input bg-background hover:bg-accent hover:text-accent-foreground h-10 py-2 px-4"
                >
                  <History className="mr-2 h-4 w-4" />
                  {showVersionHistory ? "Hide Version History" : "View Version History"}
                </button>
              </div>
            </div>
            
            <div className="bg-card rounded-lg border p-6 shadow-sm">
              <h2 className="text-lg font-semibold mb-4">Protocol Information</h2>
              <dl className="divide-y">
                <div className="py-3 flex justify-between">
                  <dt className="text-sm font-medium text-muted-foreground">Author</dt>
                  <dd className="text-sm">{protocolData.author.name}</dd>
                </div>
                <div className="py-3 flex justify-between">
                  <dt className="text-sm font-medium text-muted-foreground">Created</dt>
                  <dd className="text-sm">{new Date(protocolData.created_at).toLocaleDateString()}</dd>
                </div>
                <div className="py-3 flex justify-between">
                  <dt className="text-sm font-medium text-muted-foreground">Last Updated</dt>
                  <dd className="text-sm">{new Date(protocolData.updated_at).toLocaleDateString()}</dd>
                </div>
                <div className="py-3 flex justify-between">
                  <dt className="text-sm font-medium text-muted-foreground">Version</dt>
                  <dd className="text-sm">{protocolData.version}</dd>
                </div>
                <div className="py-3 flex justify-between">
                  <dt className="text-sm font-medium text-muted-foreground">Steps</dt>
                  <dd className="text-sm">{protocolData.steps.length}</dd>
                </div>
              </dl>
            </div>
            
            <div className="bg-card rounded-lg border p-6 shadow-sm">
              <h2 className="text-lg font-semibold mb-4">Used in Experiments</h2>
              {protocolData.used_in_experiments.length > 0 ? (
                <ul className="space-y-3">
                  {protocolData.used_in_experiments.map((exp) => (
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
          </div>
        </div>
      </div>
    </>
  );
}

// Helper function to generate a dummy protocol for comparison
function generateDummyProtocol(id) {
  const isPreviousVersion = true;
  
  return {
    id: id,
    name: 'Standard Cell Freezing Protocol',
    version: '1.1.0', // Previous version
    description: 'Basic protocol for freezing mammalian cells using DMSO as a cryoprotectant.',
    author: { name: 'Dr. John Doe', id: 'user1' },
    created_at: '2023-08-15',
    updated_at: '2023-08-15',
    duration: '1.5 hours',
    is_template: false,
    tags: ['cryopreservation', 'mammalian cells'],
    cell_types: ['Human embryonic stem cells', 'Primary fibroblasts'],
    equipment: [
      'Controlled-rate freezer',
      'Cryovials',
      'Liquid nitrogen storage tank',
      'Water bath (37°C)'
    ],
    cryoprotectants: [
      { id: '1', name: 'DMSO', concentration: '10%' }
    ],
    freezing_rate: '-0.5°C/min',
    storage_temperature: '-196°C (liquid nitrogen)',
    thawing_method: 'Rapid thawing in 37°C water bath',
    steps: [
      {
        id: 's1',
        order: 1,
        name: 'Prepare cryopreservation medium',
        description: 'Mix cell culture medium with DMSO to achieve 10% final concentration. Keep on ice.',
        duration: 15,
        duration_unit: 'minutes',
        is_critical: true
      },
      {
        id: 's2',
        order: 2,
        name: 'Prepare cells',
        description: 'Harvest cells in log phase growth (70-80% confluency). Count cells and check viability.',
        duration: 20,
        duration_unit: 'minutes',
        is_critical: true
      },
      {
        id: 's3',
        order: 3,
        name: 'Centrifuge cell suspension',
        description: 'Centrifuge at 200 x g for 5 minutes at room temperature.',
        duration: 10,
        duration_unit: 'minutes',
        is_critical: false
      },
      {
        id: 's4',
        order: 4,
        name: 'Resuspend cells in cryopreservation medium',
        description: 'Gently resuspend cell pellet in cold cryopreservation medium to achieve 1 x 10^6 cells/mL.',
        duration: 10,
        duration_unit: 'minutes',
        is_critical: true
      },
      {
        id: 's5',
        order: 5,
        name: 'Aliquot into cryovials',
        description: 'Transfer 1 mL of cell suspension to each pre-labeled cryovial.',
        duration: 15,
        duration_unit: 'minutes',
        is_critical: false
      },
      {
        id: 's6',
        order: 6,
        name: 'Controlled-rate freezing',
        description: 'Place cryovials in controlled-rate freezer and run program: cool at -0.5°C/min to -80°C.',
        duration: 60,
        duration_unit: 'minutes',
        is_critical: true
      }
    ],
    notes: 'For optimal results, cells should be in log phase of growth and have viability >90% before freezing.',
    references: [
      { title: 'Best practices for cell banking', authors: 'Johnson & Williams', journal: 'Nature Methods', year: 2019 }
    ],
    used_in_experiments: []
  };
}
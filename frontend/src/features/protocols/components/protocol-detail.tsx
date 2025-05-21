import React, { useState } from 'react';
import Link from 'next/link';
import { useRouter } from 'next/navigation';
import { getProtocolVersions } from '../actions/protocol-actions';

/**
 * Component for displaying detailed protocol information
 */
export default function ProtocolDetail({ protocol }: { protocol: any }) {
  const router = useRouter();
  const [activeTab, setActiveTab] = useState('steps');
  const [showVersionHistory, setShowVersionHistory] = useState(false);
  const [versions, setVersions] = useState<any[]>([]);
  const [versionsLoading, setVersionsLoading] = useState(false);
  const [selectedVersionId, setSelectedVersionId] = useState<string | null>(null);
  const [showComparison, setShowComparison] = useState(false);
  
  // Load version history
  const loadVersionHistory = async () => {
    try {
      setVersionsLoading(true);
      const versionData = await getProtocolVersions(protocol.id);
      setVersions(versionData);
      setShowVersionHistory(true);
    } catch (error) {
      console.error('Error loading version history:', error);
    } finally {
      setVersionsLoading(false);
    }
  };
  
  // Handle version selection
  const handleVersionSelect = (versionId: string) => {
    if (versionId === selectedVersionId) {
      setSelectedVersionId(null);
      setShowComparison(false);
    } else {
      setSelectedVersionId(versionId);
    }
  };
  
  // Toggle version history display
  const toggleVersionHistory = async () => {
    if (!showVersionHistory && versions.length === 0) {
      await loadVersionHistory();
    } else {
      setShowVersionHistory(!showVersionHistory);
    }
  };
  
  // Handle comparison
  const handleCompare = () => {
    if (selectedVersionId) {
      setShowComparison(true);
    }
  };
  
  // Render steps tab content
  const renderStepsTab = () => (
    <div className="bg-card rounded-lg border p-6 shadow-sm">
      <h2 className="text-xl font-semibold mb-4">Protocol Steps</h2>
      <ol className="space-y-5">
        {protocol.steps.map((step: any) => (
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
  );
  
  // Render materials tab content
  const renderMaterialsTab = () => (
    <div className="bg-card rounded-lg border p-6 shadow-sm">
      <h2 className="text-xl font-semibold mb-4">Materials and Equipment</h2>
      <div className="mb-6">
        <h3 className="text-base font-medium mb-2">Cryoprotectants</h3>
        <ul className="list-disc list-inside space-y-1">
          {protocol.cryoprotectants.map((cp: any, index: number) => (
            <li key={index}>
              {cp.id ? (
                <Link href={`/molecules/${cp.id}`}>
                  <span className="text-primary hover:underline">{cp.name}</span>
                </Link>
              ) : (
                <span>{cp.name}</span>
              )} ({cp.concentration})
            </li>
          ))}
        </ul>
      </div>
      <div>
        <h3 className="text-base font-medium mb-2">Equipment</h3>
        <ul className="list-disc list-inside space-y-1">
          {protocol.equipment.map((item: string, index: number) => (
            <li key={index}>{item}</li>
          ))}
        </ul>
      </div>
    </div>
  );
  
  // Render parameters tab content
  const renderParametersTab = () => (
    <div className="bg-card rounded-lg border p-6 shadow-sm">
      <h2 className="text-xl font-semibold mb-4">Parameters</h2>
      <dl className="grid grid-cols-1 md:grid-cols-2 gap-x-4 gap-y-6">
        <div>
          <dt className="text-sm font-medium text-muted-foreground">Compatible Cell Types</dt>
          <dd className="mt-1">
            <ul className="list-disc list-inside">
              {protocol.cell_types.map((cell: string, index: number) => (
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
  );
  
  // Render notes tab content
  const renderNotesTab = () => (
    <>
      <div className="bg-card rounded-lg border p-6 shadow-sm">
        <h2 className="text-xl font-semibold mb-4">Notes</h2>
        <div className="text-sm prose max-w-none">
          <p>{protocol.notes}</p>
        </div>
      </div>
      
      {protocol.references && protocol.references.length > 0 && (
        <div className="bg-card rounded-lg border p-6 shadow-sm mt-6">
          <h2 className="text-xl font-semibold mb-4">References</h2>
          <ul className="space-y-3">
            {protocol.references.map((ref: any, index: number) => (
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
  );
  
  // Render tab content based on active tab
  const renderTabContent = () => {
    switch (activeTab) {
      case 'steps':
        return renderStepsTab();
      case 'materials':
        return renderMaterialsTab();
      case 'parameters':
        return renderParametersTab();
      case 'notes':
        return renderNotesTab();
      default:
        return renderStepsTab();
    }
  };
  
  return (
    <div>
      {/* Protocol Header */}
      <div className="bg-card rounded-lg border shadow-sm p-6 mb-6">
        <div className="flex flex-col md:flex-row justify-between items-start md:items-center gap-4">
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
            
            <div className="flex flex-wrap items-center gap-x-4 gap-y-1 mt-2 text-sm text-muted-foreground">
              <span className="flex items-center">
                <svg className="w-4 h-4 mr-1" fill="none" stroke="currentColor" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M16 7a4 4 0 11-8 0 4 4 0 018 0zM12 14a7 7 0 00-7 7h14a7 7 0 00-7-7z" />
                </svg>
                Author: {protocol.author.name}
              </span>
              <span className="flex items-center">
                <svg className="w-4 h-4 mr-1" fill="none" stroke="currentColor" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 8v4l3 3m6-3a9 9 0 11-18 0 9 9 0 0118 0z" />
                </svg>
                Duration: {protocol.duration}
              </span>
              <span className="flex items-center">
                <svg className="w-4 h-4 mr-1" fill="none" stroke="currentColor" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M8 7V3m8 4V3m-9 8h10M5 21h14a2 2 0 002-2V7a2 2 0 00-2-2H5a2 2 0 00-2 2v12a2 2 0 002 2z" />
                </svg>
                Updated: {new Date(protocol.updated_at).toLocaleDateString()}
              </span>
            </div>
          </div>
          
          <div className="flex flex-wrap gap-2 mt-4 md:mt-0">
            <Link 
              href={`/experiments/create?protocol=${protocol.id}`}
              className="inline-flex items-center justify-center rounded-md bg-primary px-4 py-2 text-sm font-medium text-primary-foreground shadow transition-colors hover:bg-primary/90"
            >
              <svg className="w-4 h-4 mr-2" fill="none" stroke="currentColor" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19.428 15.428a2 2 0 00-1.022-.547l-2.387-.477a6 6 0 00-3.86.517l-.318.158a6 6 0 01-3.86.517L6.05 15.21a2 2 0 00-1.806.547M8 4h8l-1 1v5.172a2 2 0 00.586 1.414l5 5c1.26 1.26.367 3.414-1.415 3.414H4.828c-1.782 0-2.674-2.154-1.414-3.414l5-5A2 2 0 009 10.172V5L8 4z" />
              </svg>
              Create Experiment
            </Link>
            <Link 
              href={`/protocols/${protocol.id}/edit`}
              className="inline-flex items-center justify-center rounded-md bg-secondary px-4 py-2 text-sm font-medium text-secondary-foreground shadow transition-colors hover:bg-secondary/80"
            >
              <svg className="w-4 h-4 mr-2" fill="none" stroke="currentColor" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M11 5H6a2 2 0 00-2 2v11a2 2 0 002 2h11a2 2 0 002-2v-5m-1.414-9.414a2 2 0 112.828 2.828L11.828 15H9v-2.828l8.586-8.586z" />
              </svg>
              Edit
            </Link>
            <button 
              className="inline-flex items-center justify-center rounded-md bg-secondary px-4 py-2 text-sm font-medium text-secondary-foreground shadow transition-colors hover:bg-secondary/80"
            >
              <svg className="w-4 h-4 mr-2" fill="none" stroke="currentColor" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 12h6m-6 4h6m2 5H7a2 2 0 01-2-2V5a2 2 0 012-2h5.586a1 1 0 01.707.293l5.414 5.414a1 1 0 01.293.707V19a2 2 0 01-2 2z" />
              </svg>
              Duplicate
            </button>
          </div>
        </div>
      </div>
      
      {/* Version History Toggle */}
      <div className="mb-6">
        <button 
          onClick={toggleVersionHistory}
          className="flex items-center text-sm font-medium text-primary hover:underline"
        >
          <svg className="w-4 h-4 mr-1" fill="none" stroke="currentColor" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 8v4l3 3m6-3a9 9 0 11-18 0 9 9 0 0118 0z" />
          </svg>
          {versionsLoading 
            ? 'Loading version history...' 
            : showVersionHistory 
              ? 'Hide Version History' 
              : 'Show Version History'
          }
          <svg 
            className={`w-4 h-4 ml-1 transform ${showVersionHistory ? 'rotate-180' : ''}`} 
            fill="none" 
            stroke="currentColor" 
            viewBox="0 0 24 24" 
            xmlns="http://www.w3.org/2000/svg"
          >
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
          </svg>
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
                {versions.map((version) => (
                  <tr key={version.id} className="border-b hover:bg-muted/50">
                    <td className="py-2 px-3">
                      <input 
                        type="radio" 
                        name="versionSelect" 
                        checked={selectedVersionId === version.id}
                        disabled={version.id === protocol.id}
                        onChange={() => handleVersionSelect(version.id)}
                        className={version.id === protocol.id ? 'opacity-50 cursor-not-allowed' : ''}
                      />
                    </td>
                    <td className="py-2 px-3">
                      <span className={version.id === protocol.id ? 'font-semibold' : ''}>
                        v{version.version}
                        {version.id === protocol.id && ' (current)'}
                      </span>
                    </td>
                    <td className="py-2 px-3">{new Date(version.created_at).toLocaleDateString()}</td>
                    <td className="py-2 px-3">{version.created_by}</td>
                    <td className="py-2 px-3">{version.description}</td>
                    <td className="py-2 px-3">
                      {version.id !== protocol.id && (
                        <Link href={`/protocols/${version.id}`}>
                          <span className="text-primary hover:underline text-sm">View</span>
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
          {renderTabContent()}
        </div>
        
        {/* Sidebar */}
        <div className="space-y-6">
          <div className="bg-card rounded-lg border p-6 shadow-sm sticky top-6">
            <h2 className="text-lg font-semibold mb-4">Quick Actions</h2>
            <div className="space-y-2">
              <Link
                href={`/protocols/${protocol.id}/edit`}
                className="w-full inline-flex items-center justify-center rounded-md text-sm font-medium transition-colors focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-ring focus-visible:ring-offset-2 disabled:opacity-50 disabled:pointer-events-none ring-offset-background bg-secondary text-secondary-foreground hover:bg-secondary/80 h-10 py-2 px-4"
              >
                <svg className="w-4 h-4 mr-2" fill="none" stroke="currentColor" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M11 5H6a2 2 0 00-2 2v11a2 2 0 002 2h11a2 2 0 002-2v-5m-1.414-9.414a2 2 0 112.828 2.828L11.828 15H9v-2.828l8.586-8.586z" />
                </svg>
                Edit Protocol
              </Link>
              <button className="w-full inline-flex items-center justify-center rounded-md text-sm font-medium transition-colors focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-ring focus-visible:ring-offset-2 disabled:opacity-50 disabled:pointer-events-none ring-offset-background bg-secondary text-secondary-foreground hover:bg-secondary/80 h-10 py-2 px-4">
                <svg className="w-4 h-4 mr-2" fill="none" stroke="currentColor" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 16v1a3 3 0 003 3h10a3 3 0 003-3v-1m-4-4l-4 4m0 0l-4-4m4 4V4" />
                </svg>
                Download PDF
              </button>
              <button className="w-full inline-flex items-center justify-center rounded-md text-sm font-medium transition-colors focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-ring focus-visible:ring-offset-2 disabled:opacity-50 disabled:pointer-events-none ring-offset-background border border-input bg-background hover:bg-accent hover:text-accent-foreground h-10 py-2 px-4">
                <svg className="w-4 h-4 mr-2" fill="none" stroke="currentColor" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M8.684 13.342C8.886 12.938 9 12.482 9 12c0-.482-.114-.938-.316-1.342m0 2.684a3 3 0 110-2.684m0 2.684l6.632 3.316m-6.632-6l6.632-3.316m0 0a3 3 0 105.367-2.684 3 3 0 00-5.367 2.684zm0 9.316a3 3 0 105.368 2.684 3 3 0 00-5.368-2.684z" />
                </svg>
                Share Protocol
              </Link>
              <Link 
                href={`/experiments/create?protocol=${protocol.id}`}
                className="w-full inline-flex items-center justify-center rounded-md text-sm font-medium transition-colors focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-ring focus-visible:ring-offset-2 disabled:opacity-50 disabled:pointer-events-none ring-offset-background bg-primary text-primary-foreground hover:bg-primary/90 h-10 py-2 px-4"
              >
                <svg className="w-4 h-4 mr-2" fill="none" stroke="currentColor" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19.428 15.428a2 2 0 00-1.022-.547l-2.387-.477a6 6 0 00-3.86.517l-.318.158a6 6 0 01-3.86.517L6.05 15.21a2 2 0 00-1.806.547M8 4h8l-1 1v5.172a2 2 0 00.586 1.414l5 5c1.26 1.26.367 3.414-1.415 3.414H4.828c-1.782 0-2.674-2.154-1.414-3.414l5-5A2 2 0 009 10.172V5L8 4z" />
                </svg>
                Create Experiment
              </Link>
              
              <button 
                onClick={toggleVersionHistory}
                className="w-full inline-flex items-center justify-center rounded-md text-sm font-medium transition-colors focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-ring focus-visible:ring-offset-2 disabled:opacity-50 disabled:pointer-events-none ring-offset-background border border-input bg-background hover:bg-accent hover:text-accent-foreground h-10 py-2 px-4"
              >
                <svg className="w-4 h-4 mr-2" fill="none" stroke="currentColor" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 8v4l3 3m6-3a9 9 0 11-18 0 9 9 0 0118 0z" />
                </svg>
                {versionsLoading 
                  ? 'Loading...' 
                  : showVersionHistory 
                    ? 'Hide Version History' 
                    : 'View Version History'
                }
              </button>
            </div>
          </div>
          
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
            {protocol.used_in_experiments && protocol.used_in_experiments.length > 0 ? (
              <ul className="space-y-3">
                {protocol.used_in_experiments.map((exp: any) => (
                  <li key={exp.id} className="border-b pb-3 last:border-0 last:pb-0">
                    <Link href={`/experiments/${exp.id}`}>
                      <span className="text-primary hover:underline font-medium">{exp.title}</span>
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
  );
}
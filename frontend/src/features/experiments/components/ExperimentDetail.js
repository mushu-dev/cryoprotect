import React, { useState } from 'react';
import Link from 'next/link';
import ExperimentChart from './ExperimentChart';

/**
 * Enhanced experiment detail component with data visualization
 */
export default function ExperimentDetail({ experiment }) {
  const [activeTab, setActiveTab] = useState('overview');
  
  if (!experiment) {
    return <div className="p-6 text-center">Loading experiment details...</div>;
  }

  // Render tabs for different sections of detail view
  const renderTabContent = () => {
    switch (activeTab) {
      case 'overview':
        return renderOverviewTab();
      case 'results':
        return renderResultsTab();
      case 'methods':
        return renderMethodsTab();
      case 'files':
        return renderFilesTab();
      default:
        return renderOverviewTab();
    }
  };

  // Overview tab with basic experiment info
  const renderOverviewTab = () => (
    <div className="space-y-6">
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
    </div>
  );

  // Results tab with visualization
  const renderResultsTab = () => (
    <div className="space-y-6">
      <ExperimentChart experiment={experiment} />
      
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
        <h2 className="text-xl font-semibold mb-4">Statistical Analysis</h2>
        <div className="p-6 text-center text-gray-500">
          <p>Enhanced statistical analysis will be available in the next release.</p>
          <p className="text-sm mt-2">Features will include comparative analysis, trend identification, and significance testing.</p>
        </div>
      </div>
    </div>
  );

  // Methods tab with detailed protocol info
  const renderMethodsTab = () => (
    <div className="space-y-6">
      <div className="bg-card rounded-lg border p-6 shadow-sm">
        <h2 className="text-xl font-semibold mb-4">Protocol Details</h2>
        <div className="space-y-4">
          <div>
            <h3 className="text-md font-medium mb-2">
              {experiment.protocol.name} (v{experiment.protocol.version})
            </h3>
            <p className="text-gray-600">
              Standard protocol for controlled rate freezing of mammalian cells using cryoprotective agents.
            </p>
          </div>
          
          <div>
            <h3 className="text-md font-medium mb-2">Materials</h3>
            <ul className="list-disc list-inside space-y-1 text-gray-600">
              <li>Cell culture in exponential growth phase</li>
              <li>Cryopreservation medium</li>
              <li>Cryovials</li>
              <li>Controlled-rate freezing container</li>
              <li>Liquid nitrogen storage system</li>
            </ul>
          </div>
          
          <div>
            <h3 className="text-md font-medium mb-2">Procedure</h3>
            <ol className="list-decimal list-inside space-y-2 text-gray-600">
              <li>Harvest cells in exponential growth phase</li>
              <li>Prepare cryopreservation medium containing specified cryoprotectants</li>
              <li>Resuspend cells in cryopreservation medium at 1-5 × 10^6 cells/mL</li>
              <li>Transfer 1 mL aliquots to cryovials</li>
              <li>Place cryovials in controlled-rate freezing container</li>
              <li>Transfer container to -80°C freezer for 24 hours</li>
              <li>Transfer cryovials to liquid nitrogen storage</li>
            </ol>
          </div>
        </div>
      </div>
    </div>
  );

  // Files tab with attachments
  const renderFilesTab = () => (
    <div className="space-y-6">
      <div className="bg-card rounded-lg border p-6 shadow-sm">
        <h2 className="text-xl font-semibold mb-4">Attachments</h2>
        <ul className="divide-y">
          {experiment.attachments?.map(file => (
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
  );

  // Main component render
  return (
    <div>
      <div className="mb-6">
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
          </div>
        </div>
        
        {/* Tab navigation */}
        <div className="border-b">
          <nav className="flex space-x-8" aria-label="Tabs">
            {['overview', 'results', 'methods', 'files'].map((tab) => (
              <button
                key={tab}
                onClick={() => setActiveTab(tab)}
                className={`
                  py-4 px-1 border-b-2 font-medium text-sm whitespace-nowrap
                  ${activeTab === tab 
                    ? 'border-blue-500 text-blue-600' 
                    : 'border-transparent text-gray-500 hover:text-gray-700 hover:border-gray-300'}
                `}
              >
                {tab.charAt(0).toUpperCase() + tab.slice(1)}
              </button>
            ))}
          </nav>
        </div>
      </div>
      
      {/* Tab content */}
      <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
        <div className="lg:col-span-2">
          {renderTabContent()}
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
                <div className="text-sm font-medium">{experiment.date.replace(/\d+/, (m) => Number(m) - 1)}</div>
                <div className="mt-1">Thawing and recovery analysis performed</div>
              </div>
              <div className="border-l-2 border-muted pl-4 pb-4">
                <div className="text-sm font-medium">{experiment.date.replace(/\d+/, (m) => Number(m) - 5)}</div>
                <div className="mt-1">Samples transferred to long-term storage</div>
              </div>
              <div className="border-l-2 border-muted pl-4">
                <div className="text-sm font-medium">{experiment.date.replace(/\d+/, (m) => Number(m) - 7)}</div>
                <div className="mt-1">Experiment created</div>
              </div>
            </div>
          </div>
          
          <div className="bg-card rounded-lg border p-6 shadow-sm">
            <h2 className="text-xl font-semibold mb-4">Related Experiments</h2>
            <div className="space-y-3">
              <div className="rounded border p-3">
                <h3 className="font-medium mb-1">
                  <Link href={`/experiments/${Number(experiment.id) + 1}`}>
                    <a className="text-primary hover:underline">Experiment #{Number(experiment.id) + 1}</a>
                  </Link>
                </h3>
                <p className="text-sm text-muted-foreground">Similar protocol with different cell type</p>
              </div>
              <div className="rounded border p-3">
                <h3 className="font-medium mb-1">
                  <Link href={`/experiments/${Number(experiment.id) - 1 > 0 ? Number(experiment.id) - 1 : Number(experiment.id) + 2}`}>
                    <a className="text-primary hover:underline">Experiment #{Number(experiment.id) - 1 > 0 ? Number(experiment.id) - 1 : Number(experiment.id) + 2}</a>
                  </Link>
                </h3>
                <p className="text-sm text-muted-foreground">Different cryoprotectant combination</p>
              </div>
            </div>
          </div>
        </div>
      </div>
    </div>
  );
}
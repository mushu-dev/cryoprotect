import React from 'react';
import { Metadata } from 'next';
import { Tabs, TabsContent, TabsList, TabsTrigger } from '../../../components/ui/tabs';
import { Button } from '../../../components/ui/button';
import { Badge } from '../../../components/ui/badge';
import { ProtocolStepsVisualization } from '../../../features/protocols/components/protocol-steps-visualization';
import { useProtocol } from '../../../features/protocols/hooks/use-protocols';
import Link from 'next/link';

export const metadata: Metadata = {
  title: 'Protocol Details | CryoProtect',
  description: 'View detailed protocol information, steps, and versioning',
};

function formatDate(dateString: string) {
  return new Date(dateString).toLocaleDateString('en-US', {
    year: 'numeric',
    month: 'long',
    day: 'numeric',
  });
}

type ProtocolParams = {
  params: {
    id: string;
  };
};

export default function ProtocolDetailPage({ params }: ProtocolParams) {
  // This would be replaced with server components data fetching
  // For now, we'll use the client component as a wrapper
  return (
    <div className="container mx-auto px-4 py-8">
      <ProtocolDetail id={params.id} />
    </div>
  );
}

// Client component wrapper
function ProtocolDetail({ id }: { id: string }) {
  // In a real implementation, this would use server components
  // and data would be passed as props
  const { protocol, loading, error, validateProtocol } = useProtocol(id);

  if (loading) {
    return <div className="text-center py-12">Loading protocol data...</div>;
  }

  if (error || !protocol) {
    return (
      <div className="text-center py-12">
        <p className="text-red-500 mb-4">Error loading protocol</p>
        <Link href="/protocols" passHref>
          <Button variant="outline">Back to Protocols</Button>
        </Link>
      </div>
    );
  }

  return (
    <>
      <div className="flex flex-col md:flex-row justify-between items-start md:items-center mb-8">
        <div>
          <div className="flex items-center gap-2 mb-2">
            <h1 className="text-3xl font-bold">{protocol.name}</h1>
            <Badge variant="secondary">v{protocol.version}</Badge>
            {protocol.is_template && <Badge variant="outline">Template</Badge>}
          </div>
          <p className="text-muted-foreground">
            Created by {protocol.author?.name || 'Unknown'} on {formatDate(protocol.created_at)}
            {protocol.created_at !== protocol.updated_at && 
              ` • Updated on ${formatDate(protocol.updated_at)}`}
          </p>
        </div>
        
        <div className="flex gap-2 mt-4 md:mt-0">
          {protocol.is_template ? (
            <Button>Use as Template</Button>
          ) : (
            <>
              <Button variant="outline">Export</Button>
              <Button>Edit Protocol</Button>
            </>
          )}
        </div>
      </div>

      <Tabs defaultValue="steps" className="w-full">
        <TabsList className="mb-8">
          <TabsTrigger value="steps">Steps</TabsTrigger>
          <TabsTrigger value="overview">Overview</TabsTrigger>
          <TabsTrigger value="versions">Versions</TabsTrigger>
          <TabsTrigger value="validation">Validation</TabsTrigger>
        </TabsList>
        
        <TabsContent value="steps">
          <div className="bg-card rounded-lg p-6 shadow-sm">
            <h3 className="text-xl font-medium mb-6">Protocol Steps</h3>
            {protocol.steps && protocol.steps.length > 0 ? (
              <ProtocolStepsVisualization protocol={protocol} />
            ) : (
              <div className="text-center py-12">
                <p>This protocol has no steps defined.</p>
              </div>
            )}
          </div>
        </TabsContent>
        
        <TabsContent value="overview">
          <div className="grid grid-cols-1 md:grid-cols-3 gap-6 mb-6">
            <div className="bg-card rounded-lg p-6 shadow-sm">
              <h3 className="text-lg font-medium mb-2">Protocol Type</h3>
              <p>{protocol.type || 'Standard'}</p>
            </div>
            <div className="bg-card rounded-lg p-6 shadow-sm">
              <h3 className="text-lg font-medium mb-2">Duration</h3>
              <p>{protocol.duration || 'Not specified'}</p>
            </div>
            <div className="bg-card rounded-lg p-6 shadow-sm">
              <h3 className="text-lg font-medium mb-2">Step Count</h3>
              <p>{protocol.steps?.length || 0} steps</p>
            </div>
          </div>
          
          <div className="bg-card rounded-lg p-6 shadow-sm mb-6">
            <h3 className="text-lg font-medium mb-4">Description</h3>
            <p className="whitespace-pre-line">{protocol.description}</p>
          </div>
          
          <div className="bg-card rounded-lg p-6 shadow-sm">
            <h3 className="text-lg font-medium mb-4">Compatible With</h3>
            <div className="grid grid-cols-1 sm:grid-cols-2 md:grid-cols-3 gap-4">
              {protocol.compatible_with ? 
                protocol.compatible_with.map((item, index) => (
                  <div key={index} className="border rounded-md p-3">
                    <p>{item}</p>
                  </div>
                )) : 
                <p>No compatibility information specified</p>
              }
            </div>
          </div>
        </TabsContent>
        
        <TabsContent value="versions">
          {protocol.versions && protocol.versions.length > 0 ? (
            <div className="bg-card rounded-lg p-6 shadow-sm">
              <h3 className="text-xl font-medium mb-6">Version History</h3>
              <div className="border rounded-md divide-y">
                {protocol.versions.map((version, index) => (
                  <div key={index} className="p-4 flex justify-between items-center">
                    <div>
                      <div className="flex items-center gap-2">
                        <h4 className="font-medium">Version {version.version}</h4>
                        {version.version === protocol.version && (
                          <Badge variant="secondary">Current</Badge>
                        )}
                      </div>
                      <p className="text-muted-foreground">
                        Updated on {formatDate(version.date)} by {version.author}
                      </p>
                      <p className="mt-1">{version.changes}</p>
                    </div>
                    <div className="flex gap-2">
                      <Button variant="outline" size="sm">View</Button>
                      {version.version !== protocol.version && (
                        <Button variant="outline" size="sm">Compare</Button>
                      )}
                    </div>
                  </div>
                ))}
              </div>
            </div>
          ) : (
            <div className="text-center py-12 bg-card rounded-lg">
              <p>No version history is available for this protocol.</p>
            </div>
          )}
        </TabsContent>
        
        <TabsContent value="validation">
          <div className="bg-card rounded-lg p-6 shadow-sm">
            <div className="flex justify-between items-center mb-6">
              <h3 className="text-xl font-medium">Protocol Validation</h3>
              <Button onClick={() => validateProtocol()}>Validate Now</Button>
            </div>
            
            {protocol.validation_results ? (
              <div>
                <div className="mb-6">
                  <div className="flex items-center gap-2">
                    <h4 className="font-medium">Last Validation:</h4>
                    <Badge 
                      variant={protocol.validation_results.passed ? "success" : "destructive"}
                    >
                      {protocol.validation_results.passed ? "Passed" : "Failed"}
                    </Badge>
                  </div>
                  <p className="text-muted-foreground">
                    Validated on {formatDate(protocol.validation_results.date)}
                  </p>
                </div>
                
                <div>
                  <h4 className="font-medium mb-4">Validation Results</h4>
                  <div className="border rounded-md divide-y">
                    {protocol.validation_results.checks.map((check, index) => (
                      <div key={index} className="p-4">
                        <div className="flex items-start gap-2">
                          <div className={`w-6 h-6 rounded-full flex items-center justify-center flex-shrink-0 ${
                            check.passed ? 'bg-green-100 text-green-600' : 'bg-red-100 text-red-600'
                          }`}>
                            {check.passed ? '✓' : '✗'}
                          </div>
                          <div>
                            <p className="font-medium">{check.name}</p>
                            <p className="text-muted-foreground">{check.message}</p>
                          </div>
                        </div>
                      </div>
                    ))}
                  </div>
                </div>
              </div>
            ) : (
              <div className="text-center py-8">
                <p className="mb-4">This protocol has not been validated yet.</p>
                <p className="text-muted-foreground">
                  Validation checks for scientific accuracy, completeness, and consistency.
                </p>
              </div>
            )}
          </div>
        </TabsContent>
      </Tabs>
    </>
  );
}
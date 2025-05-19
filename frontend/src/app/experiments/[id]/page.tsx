import React from 'react';
import { Metadata } from 'next';
import { Tabs, TabsContent, TabsList, TabsTrigger } from '../../../components/ui/tabs';
import { Button } from '../../../components/ui/button';
import { ExperimentResultsChart } from '../../../features/experiments/components/experiment-results-chart';
import { useExperiment } from '../../../features/experiments/hooks/use-experiments';
import { Badge } from '../../../components/ui/badge';
import Link from 'next/link';

export const metadata: Metadata = {
  title: 'Experiment Details | CryoProtect',
  description: 'View detailed experiment information, results, and analysis',
};

function formatDate(dateString: string) {
  return new Date(dateString).toLocaleDateString('en-US', {
    year: 'numeric',
    month: 'long',
    day: 'numeric',
  });
}

type ExperimentParams = {
  params: {
    id: string;
  };
};

export default function ExperimentDetailPage({ params }: ExperimentParams) {
  // This would be replaced with server components data fetching
  // For now, we'll use the client component as a wrapper
  return (
    <div className="container mx-auto px-4 py-8">
      <ExperimentDetail id={params.id} />
    </div>
  );
}

// Client component wrapper
function ExperimentDetail({ id }: { id: string }) {
  // In a real implementation, this would use server components
  // and data would be passed as props
  const { experiment, loading, error } = useExperiment(id);

  if (loading) {
    return <div className="text-center py-12">Loading experiment data...</div>;
  }

  if (error || !experiment) {
    return (
      <div className="text-center py-12">
        <p className="text-red-500 mb-4">Error loading experiment</p>
        <Link href="/experiments" passHref>
          <Button variant="outline">Back to Experiments</Button>
        </Link>
      </div>
    );
  }

  return (
    <>
      <div className="flex flex-col md:flex-row justify-between items-start md:items-center mb-8">
        <div>
          <div className="flex items-center gap-2 mb-2">
            <h1 className="text-3xl font-bold">{experiment.title}</h1>
            <Badge variant={
              experiment.status === 'completed' ? 'success' : 
              experiment.status === 'in_progress' ? 'warning' : 
              experiment.status === 'planned' ? 'secondary' : 'default'
            }>
              {experiment.status}
            </Badge>
          </div>
          <p className="text-muted-foreground">
            Created by {experiment.researcher?.name || 'Unknown'} on {formatDate(experiment.created_at)}
          </p>
        </div>
        
        <div className="flex gap-2 mt-4 md:mt-0">
          <Button variant="outline">Export</Button>
          <Button variant="outline">Add Results</Button>
          <Button>Edit Experiment</Button>
        </div>
      </div>

      <Tabs defaultValue="overview" className="w-full">
        <TabsList className="mb-8">
          <TabsTrigger value="overview">Overview</TabsTrigger>
          <TabsTrigger value="results">Results</TabsTrigger>
          <TabsTrigger value="protocol">Protocol</TabsTrigger>
          <TabsTrigger value="analysis">Analysis</TabsTrigger>
        </TabsList>
        
        <TabsContent value="overview">
          <div className="grid grid-cols-1 md:grid-cols-3 gap-6 mb-6">
            <div className="bg-card rounded-lg p-6 shadow-sm">
              <h3 className="text-lg font-medium mb-2">Experiment Type</h3>
              <p>{experiment.type}</p>
            </div>
            <div className="bg-card rounded-lg p-6 shadow-sm">
              <h3 className="text-lg font-medium mb-2">Sample Type</h3>
              <p>{experiment.sample_type}</p>
            </div>
            <div className="bg-card rounded-lg p-6 shadow-sm">
              <h3 className="text-lg font-medium mb-2">Date Conducted</h3>
              <p>{experiment.execution_date ? formatDate(experiment.execution_date) : 'Not yet conducted'}</p>
            </div>
          </div>
          
          <div className="bg-card rounded-lg p-6 shadow-sm mb-6">
            <h3 className="text-lg font-medium mb-4">Description</h3>
            <p className="whitespace-pre-line">{experiment.description}</p>
          </div>
          
          <div className="bg-card rounded-lg p-6 shadow-sm">
            <h3 className="text-lg font-medium mb-4">Cryoprotectants Used</h3>
            <div className="grid grid-cols-1 sm:grid-cols-2 md:grid-cols-3 gap-4">
              {experiment.cryoprotectants?.map((cp, index) => (
                <div key={index} className="border rounded-md p-3">
                  <p className="font-medium">{cp.name}</p>
                  <p className="text-muted-foreground text-sm">{cp.concentration}% concentration</p>
                </div>
              ))}
            </div>
          </div>
        </TabsContent>
        
        <TabsContent value="results">
          {experiment.results && experiment.results.length > 0 ? (
            <div className="bg-card rounded-lg p-6 shadow-sm">
              <h3 className="text-xl font-medium mb-6">Experiment Results</h3>
              <ExperimentResultsChart results={experiment.results} />
            </div>
          ) : (
            <div className="text-center py-12 bg-card rounded-lg">
              <p className="mb-4">No results have been entered for this experiment yet.</p>
              <Button>Add Results</Button>
            </div>
          )}
        </TabsContent>
        
        <TabsContent value="protocol">
          {experiment.protocol ? (
            <div className="bg-card rounded-lg p-6 shadow-sm">
              <div className="flex justify-between items-center mb-6">
                <h3 className="text-xl font-medium">Protocol: {experiment.protocol.name}</h3>
                <Button variant="outline" size="sm">View Full Protocol</Button>
              </div>
              
              <div className="mb-6">
                <h4 className="font-medium mb-2">Protocol Summary</h4>
                <p>{experiment.protocol.description}</p>
              </div>
              
              <div>
                <h4 className="font-medium mb-4">Key Steps</h4>
                <div className="border rounded-md divide-y">
                  {experiment.protocol.steps?.map((step, index) => (
                    <div key={index} className="p-4">
                      <div className="flex justify-between items-start">
                        <div>
                          <p className="font-medium">{index + 1}. {step.name}</p>
                          <p className="text-muted-foreground">{step.description}</p>
                        </div>
                        <Badge>{step.duration}</Badge>
                      </div>
                    </div>
                  ))}
                </div>
              </div>
            </div>
          ) : (
            <div className="text-center py-12 bg-card rounded-lg">
              <p className="mb-4">No protocol is attached to this experiment.</p>
              <Button>Attach Protocol</Button>
            </div>
          )}
        </TabsContent>
        
        <TabsContent value="analysis">
          <div className="text-center py-12 bg-card rounded-lg">
            <p className="mb-4">Advanced analysis tools will be available here in a future update.</p>
            <Button variant="outline">Compare with Other Experiments</Button>
          </div>
        </TabsContent>
      </Tabs>
    </>
  );
}
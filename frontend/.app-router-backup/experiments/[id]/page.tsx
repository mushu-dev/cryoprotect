import React from 'react';
import { Metadata } from 'next';
import { Tabs, TabsContent, TabsList, TabsTrigger } from '../../../components/ui/tabs';
import { Button } from '../../../components/ui/button';
import { ExperimentResultsChart } from '../../../features/experiments/components/experiment-results-chart';
import { ResultsEntryForm } from '@/features/experiments/components/results-entry-form';
import { AnalysisDashboard } from '@/features/experiments/components/analysis-dashboard';
import { useExperiment } from '../../../features/experiments/hooks/use-experiments';
import { Badge } from '../../../components/ui/badge';
import Link from 'next/link';
import { MolstarViewer, ExperimentalDataPoint } from '../../../components/protein-visualizer/MolstarViewer';
import { Card, CardContent, CardHeader, CardTitle } from '../../../components/ui/card';
import { 
  FileDown, 
  BarChart4, 
  Plus, 
  Pencil,
  ArrowLeft
} from 'lucide-react';

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
  const { experiment, loading, error, fetchResults, addResult } = useExperiment(id);
  const [activeTab, setActiveTab] = React.useState('overview');

  React.useEffect(() => {
    // Load results if they're not already loaded
    if (experiment && (!experiment.results || experiment.results.length === 0)) {
      fetchResults();
    }
  }, [experiment, fetchResults]);

  // Prepare experimental data for visualization
  const experimentalData = React.useMemo(() => {
    if (!experiment?.results) return [];
    
    // Convert experiment results to format needed by MolstarViewer
    return experiment.results
      .filter(result => result.viability_percentage !== undefined || result.recovery_rate !== undefined)
      .map(result => {
        // Use viability percentage if available, otherwise recovery rate
        const value = result.viability_percentage !== undefined 
          ? result.viability_percentage 
          : result.recovery_rate !== undefined 
            ? result.recovery_rate 
            : 0;
            
        // Create entity ID from result information
        // This is a simplified version - in a real app, you would need to map to actual residue IDs
        const entityId = `${result.tissue_type_id}:${result.mixture_id || result.molecule_id || '0'}`;
        
        return {
          entityId,
          value,
          label: `${value.toFixed(1)}% ${result.viability_percentage !== undefined ? 'viability' : 'recovery'}`
        } as ExperimentalDataPoint;
      });
  }, [experiment?.results]);
  
  // Handle result selection in visualization
  const handleEntitySelect = (entityId: string, value?: number) => {
    if (value !== undefined) {
      console.log(`Selected entity ${entityId} with value ${value}`);
      // Here you could show a details panel or highlight the corresponding result
    }
  };

  if (loading) {
    return <div className="text-center py-12">Loading experiment data...</div>;
  }

  if (error || !experiment) {
    return (
      <div className="text-center py-12">
        <p className="text-red-500 mb-4">Error loading experiment: {error?.message || 'Experiment not found'}</p>
        <Link href="/experiments" passHref>
          <Button variant="outline">Back to Experiments</Button>
        </Link>
      </div>
    );
  }

  return (
    <>
      <div className="flex flex-col gap-4 mb-8">
        <div className="flex items-center gap-2">
          <Link href="/experiments" passHref>
            <Button variant="ghost" size="sm" className="h-8 w-8 p-0">
              <ArrowLeft className="h-4 w-4" />
            </Button>
          </Link>
          <h1 className="text-3xl font-bold">{experiment.name}</h1>
          <Badge variant={
            experiment.status === 'completed' ? 'default' : 
            experiment.status === 'in_progress' ? 'secondary' : 
            experiment.status === 'planned' ? 'outline' : 
            experiment.status === 'aborted' ? 'destructive' : 
            'default'
          }>
            {experiment.status?.replace('_', ' ')}
          </Badge>
        </div>
        
        <div className="flex flex-col md:flex-row justify-between items-start md:items-center">
          <p className="text-muted-foreground">
            Researcher: <span className="font-medium">{experiment.researcher}</span> | 
            Started: <span className="font-medium">{formatDate(experiment.start_date)}</span>
            {experiment.end_date && ` | Completed: ${formatDate(experiment.end_date)}`}
          </p>
          
          <div className="flex gap-2 mt-4 md:mt-0">
            <Button variant="outline" size="sm" className="flex items-center gap-2">
              <FileDown className="h-4 w-4" />
              Export
            </Button>
            <Button 
              onClick={() => setActiveTab('add-results')} 
              variant="outline" 
              size="sm" 
              className="flex items-center gap-2"
            >
              <Plus className="h-4 w-4" />
              Add Results
            </Button>
            <Button size="sm" className="flex items-center gap-2">
              <Pencil className="h-4 w-4" />
              Edit
            </Button>
          </div>
        </div>
      </div>

      <Tabs value={activeTab} onValueChange={setActiveTab} className="w-full">
        <TabsList className="mb-8">
          <TabsTrigger value="overview">Overview</TabsTrigger>
          <TabsTrigger value="results">Results</TabsTrigger>
          <TabsTrigger value="add-results">Add Results</TabsTrigger>
          <TabsTrigger value="visualization">Visualization</TabsTrigger>
          <TabsTrigger value="protocol">Protocol</TabsTrigger>
          <TabsTrigger value="analysis">Analysis</TabsTrigger>
        </TabsList>
        
        <TabsContent value="overview">
          <div className="grid grid-cols-1 md:grid-cols-3 gap-6 mb-6">
            <div className="bg-card rounded-lg p-6 shadow-sm">
              <h3 className="text-lg font-medium mb-2">Experiment Type</h3>
              <p>{experiment.experiment_type}</p>
            </div>
            <div className="bg-card rounded-lg p-6 shadow-sm">
              <h3 className="text-lg font-medium mb-2">Tissue Type</h3>
              <p>{experiment.tissue_type_id}</p>
            </div>
            <div className="bg-card rounded-lg p-6 shadow-sm">
              <h3 className="text-lg font-medium mb-2">Date Range</h3>
              <p>
                {formatDate(experiment.start_date)} 
                {experiment.end_date && ` - ${formatDate(experiment.end_date)}`}
              </p>
            </div>
          </div>
          
          <div className="bg-card rounded-lg p-6 shadow-sm mb-6">
            <h3 className="text-lg font-medium mb-4">Description</h3>
            <p className="whitespace-pre-line">{experiment.description || 'No description provided.'}</p>
          </div>
          
          {experiment.equipment && experiment.equipment.length > 0 && (
            <div className="bg-card rounded-lg p-6 shadow-sm mb-6">
              <h3 className="text-lg font-medium mb-4">Equipment Used</h3>
              <div className="flex flex-wrap gap-2">
                {experiment.equipment.map((item, index) => (
                  <Badge key={index} variant="outline">{item}</Badge>
                ))}
              </div>
            </div>
          )}
          
          {experiment.tags && experiment.tags.length > 0 && (
            <div className="bg-card rounded-lg p-6 shadow-sm">
              <h3 className="text-lg font-medium mb-4">Tags</h3>
              <div className="flex flex-wrap gap-2">
                {experiment.tags.map((tag, index) => (
                  <Badge key={index}>{tag}</Badge>
                ))}
              </div>
            </div>
          )}
        </TabsContent>
        
        <TabsContent value="results">
          {experiment.results && experiment.results.length > 0 ? (
            <div className="space-y-6">
              <Card>
                <CardHeader>
                  <CardTitle>Experiment Results</CardTitle>
                </CardHeader>
                <CardContent>
                  <ExperimentResultsChart results={experiment.results} />
                </CardContent>
              </Card>
              
              <Card>
                <CardHeader>
                  <CardTitle>Detailed Results</CardTitle>
                </CardHeader>
                <CardContent>
                  <div className="border rounded-md divide-y">
                    {experiment.results.map((result, index) => (
                      <div key={index} className="p-4">
                        <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
                          <div>
                            <p className="text-sm text-muted-foreground">Viability</p>
                            <p className="font-medium">
                              {result.viability_percentage !== undefined 
                                ? `${result.viability_percentage.toFixed(1)}%` 
                                : 'N/A'}
                            </p>
                          </div>
                          <div>
                            <p className="text-sm text-muted-foreground">Recovery Rate</p>
                            <p className="font-medium">
                              {result.recovery_rate !== undefined 
                                ? `${result.recovery_rate.toFixed(1)}%` 
                                : 'N/A'}
                            </p>
                          </div>
                          <div>
                            <p className="text-sm text-muted-foreground">Functionality Score</p>
                            <p className="font-medium">
                              {result.functionality_score !== undefined 
                                ? `${result.functionality_score.toFixed(1)}/10` 
                                : 'N/A'}
                            </p>
                          </div>
                        </div>
                        {result.notes && (
                          <p className="mt-2 text-sm text-muted-foreground">
                            <span className="font-medium">Notes:</span> {result.notes}
                          </p>
                        )}
                      </div>
                    ))}
                  </div>
                </CardContent>
              </Card>
            </div>
          ) : (
            <Card className="text-center py-12">
              <CardContent className="pt-12">
                <p className="mb-4">No results have been entered for this experiment yet.</p>
                <Button 
                  className="flex items-center gap-2"
                  onClick={() => setActiveTab('add-results')}
                >
                  <Plus className="h-4 w-4" />
                  Add Results
                </Button>
              </CardContent>
            </Card>
          )}
        </TabsContent>
        
        <TabsContent value="visualization">
          <Card>
            <CardHeader>
              <CardTitle>Molecular Visualization</CardTitle>
            </CardHeader>
            <CardContent>
              {experiment.results && experiment.results.length > 0 ? (
                <div className="space-y-4">
                  <MolstarViewer 
                    // For testing purposes, use a sample PDB
                    pdbId="1cbn"
                    height={500}
                    style="cartoon"
                    showControls={true}
                    name={`${experiment.name} - Molecular Visualization`}
                    backgroundColor="#f8fafc"
                    experimentalData={experimentalData}
                    experimentalDataLabel="Experimental Results"
                    experimentalDataColorScheme="blueWhiteRed"
                    onEntitySelect={handleEntitySelect}
                    showMeasurementTools={true}
                  />
                  <p className="text-sm text-muted-foreground italic text-center">
                    Note: This is a sample protein visualization. In a production environment, 
                    this would show the actual molecule structures related to your experiment.
                  </p>
                </div>
              ) : (
                <div className="text-center py-12">
                  <p className="mb-4">No molecular data available for visualization.</p>
                  <p className="text-sm text-muted-foreground">
                    Add experiment results with associated molecular data to enable visualization.
                  </p>
                </div>
              )}
            </CardContent>
          </Card>
        </TabsContent>
        
        <TabsContent value="protocol">
          {experiment.protocol ? (
            <Card>
              <CardHeader className="flex flex-row justify-between items-center">
                <CardTitle>Protocol: {experiment.protocol.name}</CardTitle>
                <Button variant="outline" size="sm">View Full Protocol</Button>
              </CardHeader>
              <CardContent className="space-y-6">
                <div>
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
              </CardContent>
            </Card>
          ) : (
            <Card className="text-center py-12">
              <CardContent className="pt-12">
                <p className="mb-4">No protocol is attached to this experiment.</p>
                <Button>Attach Protocol</Button>
              </CardContent>
            </Card>
          )}
        </TabsContent>
        
        <TabsContent value="add-results">
          <ResultsEntryForm 
            experimentId={id as any} 
            onResultsAdded={(newResults) => {
              // Refresh the results after adding new ones
              fetchResults();
              // Switch to the results tab to see the new data
              setActiveTab('results');
            }} 
          />
        </TabsContent>
        
        <TabsContent value="analysis">
          <AnalysisDashboard experimentId={id as any} />
        </TabsContent>
      </Tabs>
    </>
  );
}
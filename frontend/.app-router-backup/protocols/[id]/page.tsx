'use client'

import React, { useState } from 'react';
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs';
import { Button } from '@/components/ui/button';
import { Badge } from '@/components/ui/badge';
import { ProtocolStepsVisualization } from '@/features/protocols/components/protocol-steps-visualization';
import { useProtocol, useProtocolVersions } from '@/features/protocols/hooks';
import { ProtocolBuilder } from '@/features/protocols/components/protocol-builder';
import Link from 'next/link';
import { Loader2, Download, Edit, Share, Clock, FileText, History, AlertTriangle, CheckCircle } from 'lucide-react';
import { ProtocolStep } from '@/features/experiments/services/experiment-service';
import { Dialog, DialogContent, DialogHeader, DialogTitle } from '@/components/ui/dialog';
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card';
import { Separator } from '@/components/ui/separator';
import { AlertCircle } from 'lucide-react';

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
  return (
    <div className="container mx-auto px-4 py-8">
      <ProtocolDetail id={params.id} />
    </div>
  );
}

function ProtocolDetail({ id }: { id: string }) {
  const [isEditing, setIsEditing] = useState(false);
  const [activeTab, setActiveTab] = useState('steps');
  const [validationInProgress, setValidationInProgress] = useState(false);
  const [validationResults, setValidationResults] = useState<any>(null);
  
  const {
    protocol,
    loading: protocolLoading,
    error: protocolError,
    validateProtocol,
    exportProtocol,
    refreshProtocol
  } = useProtocol(id);
  
  const {
    versions,
    loading: versionsLoading,
    error: versionsError,
    compareVersions
  } = useProtocolVersions(id);
  
  const handleValidate = async () => {
    setValidationInProgress(true);
    const results = await validateProtocol();
    setValidationResults(results);
    setValidationInProgress(false);
    setActiveTab('validation');
  };
  
  const handleExport = async (format: 'json' | 'yaml' | 'pdf' | 'human-readable') => {
    await exportProtocol(format);
  };
  
  const handleProtocolSave = () => {
    setIsEditing(false);
    refreshProtocol();
  };
  
  // Calculate the total duration of all steps
  const getTotalDuration = () => {
    if (!protocol || !protocol.steps) return 'Not specified';
    
    const total = protocol.steps.reduce((acc, step) => {
      return acc + (step.duration || 0);
    }, 0);
    
    if (total === 0) return 'Not specified';
    
    // Get the duration unit from the first step, default to 'minutes'
    const unit = protocol.steps[0]?.duration_unit || 'minutes';
    return `${total} ${unit}`;
  };
  
  // Get the temperature range across all steps
  const getTemperatureRange = () => {
    if (!protocol || !protocol.steps) return 'Not specified';
    
    const temperaturesWithValues = protocol.steps
      .filter(step => step.temperature !== undefined)
      .map(step => step.temperature as number);
    
    if (temperaturesWithValues.length === 0) return 'Not specified';
    
    const min = Math.min(...temperaturesWithValues);
    const max = Math.max(...temperaturesWithValues);
    const unit = protocol.steps.find(step => step.temperature !== undefined)?.temperature_unit || 'C';
    
    return min === max ? `${min}°${unit}` : `${min}° to ${max}°${unit}`;
  };
  
  if (protocolLoading) {
    return (
      <div className="flex items-center justify-center py-12">
        <Loader2 className="h-8 w-8 animate-spin mr-2" />
        <span>Loading protocol data...</span>
      </div>
    );
  }

  if (protocolError || !protocol) {
    return (
      <div className="text-center py-12">
        <AlertCircle className="h-12 w-12 text-red-500 mx-auto mb-4" />
        <p className="text-red-500 text-lg mb-4">
          Error loading protocol: {protocolError?.message || 'Protocol not found'}
        </p>
        <Link href="/protocols" passHref>
          <Button variant="outline">Back to Protocols</Button>
        </Link>
      </div>
    );
  }

  if (isEditing) {
    return (
      <div>
        <div className="mb-6">
          <Link href="/protocols" className="text-sm text-muted-foreground hover:text-foreground">
            ← Back to Protocols
          </Link>
          <h1 className="text-2xl font-bold mt-2">Edit Protocol</h1>
        </div>
        
        <ProtocolBuilder 
          protocolId={id as any}
          onSave={handleProtocolSave}
          onCancel={() => setIsEditing(false)}
        />
      </div>
    );
  }

  return (
    <>
      <div className="mb-6">
        <Link href="/protocols" className="text-sm text-muted-foreground hover:text-foreground">
          ← Back to Protocols
        </Link>
      </div>
      
      <div className="flex flex-col md:flex-row justify-between items-start md:items-center mb-8">
        <div>
          <div className="flex flex-wrap items-center gap-2 mb-2">
            <h1 className="text-3xl font-bold">{protocol.name}</h1>
            <Badge variant="secondary">v{protocol.version}</Badge>
            {protocol.is_template && <Badge variant="outline">Template</Badge>}
          </div>
          <p className="text-muted-foreground">
            Created by {protocol.created_by || 'Unknown'} on {formatDate(protocol.created_at)}
            {protocol.created_at !== protocol.updated_at && 
              ` • Updated on ${formatDate(protocol.updated_at)}`}
          </p>
        </div>
        
        <div className="flex flex-wrap gap-2 mt-4 md:mt-0">
          <div className="dropdown">
            <Button variant="outline" onClick={() => handleExport('json')}>
              <Download className="h-4 w-4 mr-2" /> Export
            </Button>
          </div>
          
          <Button variant="outline" onClick={handleValidate}>
            {validationInProgress ? (
              <>
                <Loader2 className="h-4 w-4 mr-2 animate-spin" /> Validating
              </>
            ) : (
              <>
                <AlertTriangle className="h-4 w-4 mr-2" /> Validate
              </>
            )}
          </Button>
          
          {!protocol.is_template && (
            <Button onClick={() => setIsEditing(true)}>
              <Edit className="h-4 w-4 mr-2" /> Edit Protocol
            </Button>
          )}
        </div>
      </div>

      <Tabs value={activeTab} onValueChange={setActiveTab} className="w-full">
        <TabsList className="mb-8">
          <TabsTrigger value="steps" className="flex items-center gap-1">
            <FileText className="h-4 w-4" /> Steps
          </TabsTrigger>
          <TabsTrigger value="overview" className="flex items-center gap-1">
            <Clock className="h-4 w-4" /> Overview
          </TabsTrigger>
          <TabsTrigger value="versions" className="flex items-center gap-1">
            <History className="h-4 w-4" /> Versions
          </TabsTrigger>
          <TabsTrigger value="validation" className="flex items-center gap-1">
            <AlertTriangle className="h-4 w-4" /> Validation
          </TabsTrigger>
        </TabsList>
        
        <TabsContent value="steps">
          <Card>
            <CardHeader className="pb-2">
              <CardTitle>Protocol Steps</CardTitle>
            </CardHeader>
            <CardContent>
              {protocol.steps && protocol.steps.length > 0 ? (
                <ProtocolStepsVisualization protocol={protocol} />
              ) : (
                <div className="text-center py-12 border-2 border-dashed rounded-md">
                  <p>This protocol has no steps defined.</p>
                  <Button className="mt-4" onClick={() => setIsEditing(true)}>
                    <Edit className="h-4 w-4 mr-2" /> Add Steps
                  </Button>
                </div>
              )}
            </CardContent>
          </Card>
        </TabsContent>
        
        <TabsContent value="overview">
          <div className="grid grid-cols-1 md:grid-cols-3 gap-6 mb-6">
            <Card>
              <CardHeader className="pb-2">
                <CardTitle className="text-sm text-muted-foreground">PROTOCOL TYPE</CardTitle>
              </CardHeader>
              <CardContent>
                <p className="text-lg font-medium">{protocol.type || 'Standard Protocol'}</p>
              </CardContent>
            </Card>
            
            <Card>
              <CardHeader className="pb-2">
                <CardTitle className="text-sm text-muted-foreground">TOTAL DURATION</CardTitle>
              </CardHeader>
              <CardContent>
                <p className="text-lg font-medium">{getTotalDuration()}</p>
              </CardContent>
            </Card>
            
            <Card>
              <CardHeader className="pb-2">
                <CardTitle className="text-sm text-muted-foreground">TEMPERATURE RANGE</CardTitle>
              </CardHeader>
              <CardContent>
                <p className="text-lg font-medium">{getTemperatureRange()}</p>
              </CardContent>
            </Card>
          </div>
          
          <Card className="mb-6">
            <CardHeader className="pb-2">
              <CardTitle>Description</CardTitle>
            </CardHeader>
            <CardContent>
              {protocol.description ? (
                <p className="whitespace-pre-line">{protocol.description}</p>
              ) : (
                <p className="text-muted-foreground italic">No description provided</p>
              )}
            </CardContent>
          </Card>
          
          <Card>
            <CardHeader className="pb-2">
              <CardTitle>Tags</CardTitle>
            </CardHeader>
            <CardContent>
              <div className="flex flex-wrap gap-2">
                {protocol.tags && protocol.tags.length > 0 ? (
                  protocol.tags.map((tag, i) => (
                    <Badge key={i} variant="secondary">{tag}</Badge>
                  ))
                ) : (
                  <p className="text-muted-foreground italic">No tags added</p>
                )}
              </div>
            </CardContent>
          </Card>
        </TabsContent>
        
        <TabsContent value="versions">
          <Card>
            <CardHeader className="pb-2">
              <CardTitle>Version History</CardTitle>
            </CardHeader>
            <CardContent>
              {versionsLoading ? (
                <div className="flex items-center justify-center py-12">
                  <Loader2 className="h-6 w-6 animate-spin mr-2" />
                  <span>Loading version history...</span>
                </div>
              ) : versionsError ? (
                <div className="text-center py-8 text-red-500">
                  <p>Error loading version history</p>
                </div>
              ) : versions && versions.length > 0 ? (
                <div className="border rounded-md divide-y">
                  {versions.map((version, index) => (
                    <div key={index} className="p-4 flex justify-between items-center">
                      <div>
                        <div className="flex items-center gap-2">
                          <h4 className="font-medium">Version {version.version}</h4>
                          {version.id === protocol.id && (
                            <Badge variant="secondary">Current</Badge>
                          )}
                        </div>
                        <p className="text-muted-foreground">
                          Updated on {formatDate(version.created_at)} by {version.created_by}
                        </p>
                      </div>
                      <div className="flex gap-2">
                        <Button variant="outline" size="sm">View</Button>
                        {version.id !== protocol.id && (
                          <Button variant="outline" size="sm">Compare</Button>
                        )}
                      </div>
                    </div>
                  ))}
                </div>
              ) : (
                <div className="text-center py-12">
                  <p className="text-muted-foreground">
                    No version history is available for this protocol.
                  </p>
                </div>
              )}
            </CardContent>
          </Card>
        </TabsContent>
        
        <TabsContent value="validation">
          <Card>
            <CardHeader className="pb-2 flex flex-row items-center justify-between">
              <CardTitle>Protocol Validation</CardTitle>
              <Button 
                onClick={handleValidate}
                disabled={validationInProgress}
              >
                {validationInProgress ? (
                  <>
                    <Loader2 className="h-4 w-4 mr-2 animate-spin" /> Validating
                  </>
                ) : (
                  'Validate Now'
                )}
              </Button>
            </CardHeader>
            <CardContent>
              {validationInProgress ? (
                <div className="flex items-center justify-center py-12">
                  <Loader2 className="h-8 w-8 animate-spin mr-2" />
                  <span>Validating protocol...</span>
                </div>
              ) : validationResults ? (
                <div>
                  <div className="flex items-center gap-2 mb-6">
                    <div className={`p-2 rounded-full ${
                      validationResults.valid ? 'bg-green-100' : 'bg-red-100'
                    }`}>
                      {validationResults.valid ? (
                        <CheckCircle className="h-6 w-6 text-green-600" />
                      ) : (
                        <AlertTriangle className="h-6 w-6 text-red-600" />
                      )}
                    </div>
                    <div>
                      <p className="font-medium">
                        {validationResults.valid 
                          ? 'Protocol is valid' 
                          : `Protocol validation failed with ${validationResults.errors.length} errors`
                        }
                      </p>
                      <p className="text-sm text-muted-foreground">
                        Last validated just now
                      </p>
                    </div>
                  </div>
                  
                  {!validationResults.valid && validationResults.errors.length > 0 && (
                    <div className="space-y-4 mt-4">
                      <h3 className="font-medium">Validation Issues</h3>
                      <div className="space-y-2">
                        {validationResults.errors.map((error: any, i: number) => (
                          <div 
                            key={i} 
                            className="p-3 border border-red-200 bg-red-50 rounded-md flex gap-2"
                          >
                            <AlertCircle className="h-5 w-5 text-red-500 flex-shrink-0 mt-0.5" />
                            <div>
                              <p className="font-medium text-red-700">
                                {error.step_id 
                                  ? `Error in step ${protocol.steps?.find(s => s.id === error.step_id)?.order || 'unknown'}`
                                  : error.field 
                                  ? `Error in ${error.field}`
                                  : 'Protocol Error'
                                }
                              </p>
                              <p className="text-red-700">{error.message}</p>
                            </div>
                          </div>
                        ))}
                      </div>
                    </div>
                  )}
                  
                  {validationResults.warnings && validationResults.warnings.length > 0 && (
                    <div className="space-y-4 mt-4">
                      <h3 className="font-medium">Warnings</h3>
                      <div className="space-y-2">
                        {validationResults.warnings.map((warning: any, i: number) => (
                          <div 
                            key={i} 
                            className="p-3 border border-amber-200 bg-amber-50 rounded-md flex gap-2"
                          >
                            <AlertTriangle className="h-5 w-5 text-amber-500 flex-shrink-0 mt-0.5" />
                            <div>
                              <p className="font-medium text-amber-700">
                                {warning.step_id 
                                  ? `Warning in step ${protocol.steps?.find(s => s.id === warning.step_id)?.order || 'unknown'}`
                                  : warning.field 
                                  ? `Warning for ${warning.field}`
                                  : 'Protocol Warning'
                                }
                              </p>
                              <p className="text-amber-700">{warning.message}</p>
                            </div>
                          </div>
                        ))}
                      </div>
                    </div>
                  )}
                  
                  {validationResults.suggestions && validationResults.suggestions.length > 0 && (
                    <div className="space-y-4 mt-4">
                      <h3 className="font-medium">Suggestions</h3>
                      <div className="space-y-2">
                        {validationResults.suggestions.map((suggestion: any, i: number) => (
                          <div 
                            key={i} 
                            className="p-3 border border-blue-200 bg-blue-50 rounded-md flex gap-2"
                          >
                            <AlertCircle className="h-5 w-5 text-blue-500 flex-shrink-0 mt-0.5" />
                            <div>
                              <p className="font-medium text-blue-700">
                                {suggestion.step_id 
                                  ? `Suggestion for step ${protocol.steps?.find(s => s.id === suggestion.step_id)?.order || 'unknown'}`
                                  : suggestion.field 
                                  ? `Suggestion for ${suggestion.field}`
                                  : 'Protocol Suggestion'
                                }
                              </p>
                              <p className="text-blue-700">{suggestion.message}</p>
                              {suggestion.suggested_value !== undefined && (
                                <p className="text-sm text-blue-700 mt-1">
                                  Suggested value: {suggestion.suggested_value}
                                </p>
                              )}
                            </div>
                          </div>
                        ))}
                      </div>
                    </div>
                  )}
                </div>
              ) : (
                <div className="text-center py-8">
                  <p className="mb-4">This protocol has not been validated yet.</p>
                  <p className="text-muted-foreground">
                    Validation checks for scientific accuracy, completeness, and consistency.
                  </p>
                </div>
              )}
            </CardContent>
          </Card>
        </TabsContent>
      </Tabs>
    </>
  );
}
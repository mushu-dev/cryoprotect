/**
 * Protocol Template Comparison Component
 * For comparing different versions of a protocol template
 */

import React, { useState, useEffect } from 'react';
import { useTemplateVersions, useTemplateComparison } from '../hooks/use-convex-protocol-templates';
import { useProtocolTemplates } from '../hooks/use-convex-protocol-templates';
import { Card, CardContent, CardHeader, CardTitle, CardDescription, CardFooter } from '../../ui/card';
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from '../../ui/select';
import { Alert, AlertDescription, AlertTitle } from '../../ui/alert';
import { Button } from '../../ui/button';
import { Badge } from '../../ui/badge';
import { useToast } from '../../ui/use-toast';
import { Id } from "../../../../convex/_generated/dataModel";

export interface TemplateComparisonProps {
  templateId: string;
}

export function TemplateComparison({ templateId }: TemplateComparisonProps) {
  const { toast } = useToast();
  const { getTemplate } = useProtocolTemplates();
  const { compareTemplates } = useTemplateComparison();
  
  const [template, setTemplate] = useState<any>(null);
  const [loading, setLoading] = useState(true);
  const [comparison, setComparison] = useState<any>(null);
  const [compareVersion, setCompareVersion] = useState<string>('');
  
  // Load template and its versions
  useEffect(() => {
    const loadTemplateData = async () => {
      try {
        setLoading(true);
        const templateData = await getTemplate({ protocolId: templateId as Id<"protocols"> });
        
        if (templateData) {
          setTemplate(templateData);
        }
      } catch (error) {
        console.error('Error loading template:', error);
        toast({
          title: 'Error',
          description: 'Failed to load template',
          variant: 'destructive'
        });
      } finally {
        setLoading(false);
      }
    };
    
    loadTemplateData();
  }, [templateId, getTemplate, toast]);
  
  // Get all versions of this template
  const { versions, loading: versionsLoading } = useTemplateVersions(template?.name);
  
  // Handle version selection for comparison
  const handleVersionSelect = async (versionId: string) => {
    if (versionId === templateId) {
      setComparison(null);
      return;
    }
    
    setCompareVersion(versionId);
    
    try {
      const comparisonData = await compareTemplates(templateId as Id<"protocols">, versionId as Id<"protocols">);
      setComparison(comparisonData);
    } catch (error) {
      console.error('Error comparing templates:', error);
      toast({
        title: 'Error',
        description: 'Failed to compare template versions',
        variant: 'destructive'
      });
    }
  };
  
  if (loading || versionsLoading) {
    return (
      <div className="py-8 text-center">
        <div className="animate-spin rounded-full h-12 w-12 border-b-2 border-primary mx-auto"></div>
        <p className="mt-4 text-muted-foreground">Loading template data...</p>
      </div>
    );
  }
  
  if (!template) {
    return (
      <div className="py-8 text-center">
        <p className="text-muted-foreground">Template not found</p>
      </div>
    );
  }
  
  const otherVersions = versions.filter(v => v._id.toString() !== templateId);
  
  return (
    <div className="space-y-6">
      <div className="flex flex-col md:flex-row justify-between gap-4">
        <div>
          <h2 className="text-2xl font-bold mb-2">Template Version Comparison</h2>
          <p className="text-muted-foreground">
            Compare current version with previous versions of this template
          </p>
        </div>
        
        <div className="w-full md:w-[240px]">
          <Select
            value={compareVersion}
            onValueChange={handleVersionSelect}
          >
            <SelectTrigger>
              <SelectValue placeholder="Select version to compare" />
            </SelectTrigger>
            <SelectContent>
              {otherVersions.length === 0 ? (
                <SelectItem value="" disabled>No other versions available</SelectItem>
              ) : (
                <>
                  <SelectItem value="">Select version</SelectItem>
                  {otherVersions.map(version => (
                    <SelectItem key={version._id.toString()} value={version._id.toString()}>
                      Version {version.version} ({new Date(version.createdAt).toLocaleDateString()})
                    </SelectItem>
                  ))}
                </>
              )}
            </SelectContent>
          </Select>
        </div>
      </div>
      
      {compareVersion && comparison ? (
        <div className="space-y-6">
          <Card>
            <CardHeader>
              <CardTitle>Steps Comparison</CardTitle>
              <CardDescription>
                Changes between Version {template.version} and selected version
              </CardDescription>
            </CardHeader>
            <CardContent>
              {comparison.steps.added.length === 0 && 
               comparison.steps.removed.length === 0 && 
               comparison.steps.modified.length === 0 ? (
                <div className="py-4 text-center">
                  <p className="text-muted-foreground">No step changes found</p>
                </div>
              ) : (
                <div className="space-y-4">
                  {comparison.steps.added.length > 0 && (
                    <div>
                      <h3 className="text-md font-semibold mb-2">Added Steps</h3>
                      <ul className="list-disc list-inside space-y-1 ml-4">
                        {comparison.steps.added.map((step: any, index: number) => (
                          <li key={index} className="text-green-600">
                            {step.name}
                          </li>
                        ))}
                      </ul>
                    </div>
                  )}
                  
                  {comparison.steps.removed.length > 0 && (
                    <div>
                      <h3 className="text-md font-semibold mb-2">Removed Steps</h3>
                      <ul className="list-disc list-inside space-y-1 ml-4">
                        {comparison.steps.removed.map((step: any, index: number) => (
                          <li key={index} className="text-red-600">
                            {step.name}
                          </li>
                        ))}
                      </ul>
                    </div>
                  )}
                  
                  {comparison.steps.modified.length > 0 && (
                    <div>
                      <h3 className="text-md font-semibold mb-2">Modified Steps</h3>
                      <div className="space-y-2">
                        {comparison.steps.modified.map((modified: any, index: number) => (
                          <Card key={index} className="bg-muted/30">
                            <CardHeader className="pb-3">
                              <CardTitle className="text-md">{modified.name}</CardTitle>
                            </CardHeader>
                            <CardContent>
                              <ul className="list-disc list-inside space-y-1 ml-4">
                                {modified.changes.map((change: any, cIndex: number) => (
                                  <li key={cIndex} className="text-amber-600">
                                    {change}
                                  </li>
                                ))}
                              </ul>
                            </CardContent>
                          </Card>
                        ))}
                      </div>
                    </div>
                  )}
                </div>
              )}
            </CardContent>
          </Card>
          
          <Card>
            <CardHeader>
              <CardTitle>Parameters Comparison</CardTitle>
              <CardDescription>
                Parameter changes between versions
              </CardDescription>
            </CardHeader>
            <CardContent>
              {comparison.parameters.added.length === 0 && 
               comparison.parameters.removed.length === 0 && 
               comparison.parameters.modified.length === 0 ? (
                <div className="py-4 text-center">
                  <p className="text-muted-foreground">No parameter changes found</p>
                </div>
              ) : (
                <div className="space-y-4">
                  {comparison.parameters.added.length > 0 && (
                    <div>
                      <h3 className="text-md font-semibold mb-2">Added Parameters</h3>
                      <ul className="list-disc list-inside space-y-1 ml-4">
                        {comparison.parameters.added.map((param: string, index: number) => (
                          <li key={index} className="text-green-600">
                            {param}
                          </li>
                        ))}
                      </ul>
                    </div>
                  )}
                  
                  {comparison.parameters.removed.length > 0 && (
                    <div>
                      <h3 className="text-md font-semibold mb-2">Removed Parameters</h3>
                      <ul className="list-disc list-inside space-y-1 ml-4">
                        {comparison.parameters.removed.map((param: string, index: number) => (
                          <li key={index} className="text-red-600">
                            {param}
                          </li>
                        ))}
                      </ul>
                    </div>
                  )}
                  
                  {comparison.parameters.modified.length > 0 && (
                    <div>
                      <h3 className="text-md font-semibold mb-2">Modified Parameters</h3>
                      <div className="space-y-2">
                        {comparison.parameters.modified.map((modified: any, index: number) => (
                          <div key={index} className="border rounded-md p-3">
                            <div className="font-semibold">{modified.key}</div>
                            <div className="grid grid-cols-2 gap-4 mt-2">
                              <div className="space-y-1">
                                <div className="text-xs text-muted-foreground">From:</div>
                                <div className="text-sm p-2 bg-muted rounded-md">
                                  {JSON.stringify(modified.from)}
                                </div>
                              </div>
                              <div className="space-y-1">
                                <div className="text-xs text-muted-foreground">To:</div>
                                <div className="text-sm p-2 bg-muted rounded-md">
                                  {JSON.stringify(modified.to)}
                                </div>
                              </div>
                            </div>
                          </div>
                        ))}
                      </div>
                    </div>
                  )}
                </div>
              )}
            </CardContent>
          </Card>
          
          <div className="flex justify-end">
            <Button variant="outline" onClick={() => setCompareVersion('')}>
              Close Comparison
            </Button>
          </div>
        </div>
      ) : (
        <Card>
          <CardContent className="py-8 text-center">
            <p className="text-muted-foreground">
              {otherVersions.length === 0 
                ? "No other versions available for comparison" 
                : "Select a version to compare with current template"}
            </p>
          </CardContent>
        </Card>
      )}
    </div>
  );
}

export default TemplateComparison;
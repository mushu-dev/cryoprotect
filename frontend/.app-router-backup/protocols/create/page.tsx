import React from 'react';
import { Metadata } from 'next';
import { Button } from '../../../components/ui/button';
import { ProtocolBuilder } from '../../../features/protocols/components/protocol-builder';
import Link from 'next/link';
import { Alert, AlertDescription, AlertTitle } from '../../../components/ui/alert';
import { AlertCircle } from 'lucide-react';
import { Tabs, TabsContent, TabsList, TabsTrigger } from '../../../components/ui/tabs';

export const metadata: Metadata = {
  title: 'Create Protocol | CryoProtect',
  description: 'Create a new standardized protocol for cryopreservation experiments',
};

export default function CreateProtocolPage() {
  return (
    <div className="container mx-auto px-4 py-8">
      <div className="flex flex-col md:flex-row justify-between items-start md:items-center mb-8">
        <div>
          <h1 className="text-3xl font-bold mb-2">Create New Protocol</h1>
          <p className="text-muted-foreground">
            Design a standardized protocol for cryopreservation procedures
          </p>
        </div>
        
        <div className="mt-4 md:mt-0">
          <Link href="/protocols" passHref>
            <Button variant="outline">Cancel</Button>
          </Link>
        </div>
      </div>

      <ProtocolCreationWizard />
    </div>
  );
}

// Client component wrapper
function ProtocolCreationWizard() {
  'use client';
  
  const [activeTab, setActiveTab] = React.useState('create');
  const [error, setError] = React.useState<string | null>(null);
  const [success, setSuccess] = React.useState(false);
  const [protocolId, setProtocolId] = React.useState<string | null>(null);
  
  // In a real implementation, this would handle form submission logic
  const handleSave = (protocolData: any) => {
    try {
      console.log('Saving protocol:', protocolData);
      // Simulating successful protocol creation
      setError(null);
      setSuccess(true);
      setProtocolId('123456'); // Would be a real ID from the API
      
      // In a real implementation, we would redirect to the newly created protocol
      // window.location.href = `/protocols/${newProtocol.id}`;
    } catch (err) {
      setError('Failed to create protocol. Please try again.');
      setSuccess(false);
    }
  };
  
  // Called when user imports a protocol from a file
  const handleImportProtocol = (file: File) => {
    // This would use the protocol service to import a protocol
    try {
      console.log('Importing protocol from file:', file.name);
      setError('File import is not implemented in this demo.');
    } catch (err) {
      setError('Failed to import protocol. Please check file format and try again.');
    }
  };
  
  // Called when user wants to use a template
  const handleUseTemplate = (templateId: string) => {
    // This would fetch the template and populate the protocol builder
    console.log('Using template:', templateId);
    setActiveTab('create');
  };
  
  const handleCancel = () => {
    // Navigate back to protocols list
    window.location.href = '/protocols';
  };

  return (
    <>
      {error && (
        <Alert variant="destructive" className="mb-6">
          <AlertCircle className="h-4 w-4" />
          <AlertTitle>Error</AlertTitle>
          <AlertDescription>{error}</AlertDescription>
        </Alert>
      )}
      
      {success && !error && (
        <Alert className="mb-6 border-green-500 text-green-700">
          <AlertTitle>Success</AlertTitle>
          <AlertDescription>
            Protocol created successfully!{' '}
            <Link href={`/protocols/${protocolId}`} className="font-medium underline">
              View protocol
            </Link>
          </AlertDescription>
        </Alert>
      )}
    
      <Tabs defaultValue="create" value={activeTab} onValueChange={setActiveTab}>
        <TabsList className="mb-6">
          <TabsTrigger value="create">Create From Scratch</TabsTrigger>
          <TabsTrigger value="template">Use Template</TabsTrigger>
          <TabsTrigger value="import">Import</TabsTrigger>
        </TabsList>
        
        <TabsContent value="create">
          <div className="bg-card rounded-lg shadow-sm overflow-hidden">
            <ProtocolBuilder 
              onSave={handleSave}
              onCancel={handleCancel}
            />
          </div>
        </TabsContent>
        
        <TabsContent value="template">
          <div className="bg-card rounded-lg shadow-sm p-6">
            <h2 className="text-xl font-semibold mb-4">Protocol Templates</h2>
            <p className="text-muted-foreground mb-6">
              Select a template to use as a starting point for your protocol. You can modify it to suit your specific needs.
            </p>
            
            <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4">
              {/* Template cards - in a real implementation these would be dynamically loaded */}
              {['Standard Slow Freezing', 'Rapid Vitrification', 'Cell Line Vitrification', 'Tissue Transport'].map((template, index) => (
                <div key={index} className="border rounded-lg p-4 hover:border-primary cursor-pointer" onClick={() => handleUseTemplate(`template-${index}`)}>
                  <h3 className="font-medium">{template}</h3>
                  <p className="text-sm text-muted-foreground">Standard protocol for {template.toLowerCase()}</p>
                  <p className="text-xs mt-2">8-12 steps • 2-3 hours</p>
                </div>
              ))}
            </div>
          </div>
        </TabsContent>
        
        <TabsContent value="import">
          <div className="bg-card rounded-lg shadow-sm p-6">
            <h2 className="text-xl font-semibold mb-4">Import Protocol</h2>
            <p className="text-muted-foreground mb-6">
              Import a protocol from a file. Supported formats: JSON, YAML, or CSV.
            </p>
            
            <div className="border-2 border-dashed border-muted-foreground/20 rounded-lg p-8 text-center">
              <input
                type="file"
                id="protocol-file"
                className="hidden"
                accept=".json,.yaml,.yml,.csv"
                onChange={(e) => e.target.files && handleImportProtocol(e.target.files[0])}
              />
              <label htmlFor="protocol-file" className="cursor-pointer">
                <div className="flex flex-col items-center justify-center gap-2">
                  <div className="rounded-full bg-primary/10 p-2">
                    <svg xmlns="http://www.w3.org/2000/svg" width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" className="text-primary">
                      <path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4"></path>
                      <polyline points="17 8 12 3 7 8"></polyline>
                      <line x1="12" y1="3" x2="12" y2="15"></line>
                    </svg>
                  </div>
                  <p className="font-medium text-primary">Click to upload a file</p>
                  <p className="text-sm text-muted-foreground">or drag and drop</p>
                </div>
              </label>
            </div>
            
            <div className="mt-6">
              <h3 className="text-md font-medium mb-2">Protocol Import Format</h3>
              <p className="text-sm text-muted-foreground mb-2">
                The protocol file should have the following structure:
              </p>
              <pre className="bg-muted p-4 rounded-md text-xs overflow-auto">
{`{
  "name": "Protocol Name",
  "description": "Protocol description",
  "steps": [
    {
      "name": "Step 1",
      "description": "Step description",
      "order": 1,
      "duration": 60,
      "duration_unit": "minutes"
    }
  ]
}`}
              </pre>
            </div>
          </div>
        </TabsContent>
      </Tabs>
    </>
  );
}
import React, { useState, useEffect } from 'react';
import { useRouter } from 'next/router';
import Head from 'next/head';
import Link from 'next/link';
import { ProtocolBuilder } from '../../../features/protocols/components/protocol-builder';
import { 
  useProtocol, 
  useProtocolVersions 
} from '../../../features/protocols/hooks/use-protocols';
import { 
  ArrowLeft, 
  Clock, 
  FileClock, 
  FileText, 
  History, 
  Users 
} from 'lucide-react';
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { Input } from '@/components/ui/input';
import { Label } from '@/components/ui/label';
import { Textarea } from '@/components/ui/textarea';

export default function NewVersionPage() {
  const router = useRouter();
  const { id } = router.query;
  const { protocol, loading: protocolLoading } = useProtocol(id);
  const { versions, loading: versionsLoading, createVersion } = useProtocolVersions(id);
  
  const [versionNotes, setVersionNotes] = useState('');
  const [submitting, setSubmitting] = useState(false);
  const [error, setError] = useState(null);
  
  // Calculate next version number based on current version
  const [nextVersion, setNextVersion] = useState('');
  
  useEffect(() => {
    if (protocol) {
      // Parse the current version (assuming semver format: major.minor.patch)
      const parts = protocol.version.split('.');
      if (parts.length === 3) {
        const major = parseInt(parts[0]);
        const minor = parseInt(parts[1]);
        const patch = parseInt(parts[2]);
        
        // Increment minor version by default
        setNextVersion(`${major}.${minor + 1}.0`);
      } else {
        // Fallback if version format is unexpected
        setNextVersion(`${protocol.version}_new`);
      }
    }
  }, [protocol]);
  
  const handleVersionNotesChange = (e) => {
    setVersionNotes(e.target.value);
  };
  
  const handleNextVersionChange = (e) => {
    setNextVersion(e.target.value);
  };
  
  const handleCreateVersion = async (updatedProtocol) => {
    try {
      setSubmitting(true);
      setError(null);
      
      // Include version notes and specific protocol changes
      const changes = {
        ...updatedProtocol,
        version: nextVersion,
      };
      
      const newVersionProtocol = await createVersion(changes, versionNotes);
      
      if (newVersionProtocol) {
        // Redirect to the new protocol version page
        router.push(`/protocols/${newVersionProtocol.id}`);
      }
    } catch (err) {
      console.error('Error creating new version:', err);
      setError(err.message || 'Failed to create new protocol version');
    } finally {
      setSubmitting(false);
    }
  };
  
  const handleCancel = () => {
    router.back();
  };
  
  if (protocolLoading) {
    return (
      <div className="container mx-auto px-4 py-8">
        <div className="flex justify-center items-center h-64">
          <div className="text-center">
            <div className="animate-spin rounded-full h-12 w-12 border-t-2 border-b-2 border-primary mx-auto mb-4"></div>
            <p className="text-muted-foreground">Loading protocol...</p>
          </div>
        </div>
      </div>
    );
  }
  
  if (!protocol) {
    return (
      <div className="container mx-auto px-4 py-8">
        <div className="text-center py-12">
          <h2 className="text-2xl font-bold mb-2">Protocol Not Found</h2>
          <p className="text-muted-foreground mb-6">The protocol you're looking for doesn't exist or you don't have permission to access it.</p>
          <Button asChild>
            <Link href="/protocols">Back to Protocols</Link>
          </Button>
        </div>
      </div>
    );
  }
  
  return (
    <>
      <Head>
        <title>Create New Version - {protocol.name} | CryoProtect</title>
        <meta name="description" content={`Create a new version of ${protocol.name}`} />
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
            <li>
              <div className="flex items-center">
                <svg className="w-3 h-3 text-muted-foreground mx-1" aria-hidden="true" xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 6 10">
                  <path stroke="currentColor" strokeLinecap="round" strokeLinejoin="round" strokeWidth="2" d="m1 9 4-4-4-4"/>
                </svg>
                <Link href={`/protocols/${id}`}>
                  <a className="ml-1 text-sm font-medium text-muted-foreground hover:text-foreground md:ml-2 truncate max-w-[120px]">
                    {protocol.name}
                  </a>
                </Link>
              </div>
            </li>
            <li aria-current="page">
              <div className="flex items-center">
                <svg className="w-3 h-3 text-muted-foreground mx-1" aria-hidden="true" xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 6 10">
                  <path stroke="currentColor" strokeLinecap="round" strokeLinejoin="round" strokeWidth="2" d="m1 9 4-4-4-4"/>
                </svg>
                <span className="ml-1 text-sm font-medium text-foreground md:ml-2">
                  New Version
                </span>
              </div>
            </li>
          </ol>
        </nav>
        
        {/* Page Header */}
        <div className="mb-6">
          <div className="flex items-center gap-2">
            <Button 
              variant="outline" 
              size="sm"
              onClick={handleCancel}
              className="mb-4"
            >
              <ArrowLeft className="h-4 w-4 mr-1" />
              Back
            </Button>
          </div>
          
          <h1 className="text-3xl font-bold">Create New Version</h1>
          <p className="text-muted-foreground mt-1">
            Creating a new version of <span className="font-medium">{protocol.name}</span>
          </p>
        </div>
        
        {/* Version Information */}
        <Card className="mb-6">
          <CardHeader>
            <CardTitle>Version Information</CardTitle>
          </CardHeader>
          <CardContent className="space-y-4">
            <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
              <div className="space-y-4">
                <div className="flex items-center gap-4">
                  <div className="flex items-center gap-2 px-4 py-2 bg-muted rounded">
                    <History className="h-4 w-4" />
                    <span className="font-medium">Current Version: {protocol.version}</span>
                  </div>
                  
                  <div className="flex items-center gap-2 px-4 py-2 bg-primary/10 text-primary rounded">
                    <History className="h-4 w-4" />
                    <span className="font-medium">New Version: {nextVersion}</span>
                  </div>
                </div>
                
                <div className="space-y-2">
                  <Label htmlFor="next-version">Version Number</Label>
                  <Input
                    id="next-version"
                    value={nextVersion}
                    onChange={handleNextVersionChange}
                    placeholder="1.0.0"
                  />
                  <p className="text-sm text-muted-foreground">
                    Specify the version number for this new version
                  </p>
                </div>
                
                <div className="space-y-2">
                  <Label htmlFor="version-notes">Version Notes</Label>
                  <Textarea
                    id="version-notes"
                    value={versionNotes}
                    onChange={handleVersionNotesChange}
                    placeholder="Describe what changed in this version..."
                    rows={4}
                  />
                  <p className="text-sm text-muted-foreground">
                    Include notes about what has changed in this version, such as improvements, fixes, or new features.
                  </p>
                </div>
              </div>
              
              <div className="space-y-4">
                <div className="p-4 border rounded-md bg-amber-50 border-amber-200">
                  <h3 className="font-medium text-amber-800 mb-2">About Creating New Versions</h3>
                  <ul className="list-disc list-inside space-y-1 text-sm text-amber-700">
                    <li>This will create a new version of the protocol, preserving the original.</li>
                    <li>You can modify any aspect of the protocol in the new version.</li>
                    <li>The original version will remain unchanged and accessible.</li>
                    <li>Version comparison will be available to track changes between versions.</li>
                    <li>Consider following semantic versioning (major.minor.patch) for version numbers.</li>
                  </ul>
                </div>
                
                <div className="p-4 border rounded-md">
                  <h3 className="font-medium mb-2">Protocol Information</h3>
                  <div className="space-y-2 text-sm">
                    <div className="flex items-center gap-2">
                      <FileText className="h-4 w-4 text-muted-foreground" />
                      <span className="text-muted-foreground">Name:</span>
                      <span>{protocol.name}</span>
                    </div>
                    <div className="flex items-center gap-2">
                      <Users className="h-4 w-4 text-muted-foreground" />
                      <span className="text-muted-foreground">Author:</span>
                      <span>{protocol.author?.name || 'Unknown'}</span>
                    </div>
                    <div className="flex items-center gap-2">
                      <Clock className="h-4 w-4 text-muted-foreground" />
                      <span className="text-muted-foreground">Created:</span>
                      <span>{new Date(protocol.created_at).toLocaleDateString()}</span>
                    </div>
                    <div className="flex items-center gap-2">
                      <FileClock className="h-4 w-4 text-muted-foreground" />
                      <span className="text-muted-foreground">Last Updated:</span>
                      <span>{new Date(protocol.updated_at).toLocaleDateString()}</span>
                    </div>
                  </div>
                </div>
              </div>
            </div>
            
            {error && (
              <div className="p-4 border rounded-md bg-red-50 border-red-200 text-red-600">
                {error}
              </div>
            )}
          </CardContent>
        </Card>
        
        {/* Protocol Builder */}
        <div className="mb-6">
          <h2 className="text-xl font-bold mb-4">Modify Protocol</h2>
          <ProtocolBuilder 
            protocolId={id}
            onSave={handleCreateVersion}
            onCancel={handleCancel}
          />
        </div>
      </div>
    </>
  );
}
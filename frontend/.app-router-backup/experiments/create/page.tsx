'use client'

import { useState } from 'react'
import { useRouter, useSearchParams } from 'next/navigation'
import Link from 'next/link'
import { ArrowLeft, ClipboardList, Upload } from 'lucide-react'
import { Button } from '@/components/ui/button'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { Input } from '@/components/ui/input'
import { Tabs, TabsList, TabsTrigger, TabsContent } from '@/components/ui/tabs'
import { useExperimentCreation } from '@/features/experiments/hooks/use-experiments'
import { ExperimentCreationWizard } from '@/features/experiments/components/experiment-creation-wizard'

export default function CreateExperimentPage() {
  const router = useRouter()
  const searchParams = useSearchParams()
  const protocolId = searchParams.get('protocol_id')
  const tissueTypeId = searchParams.get('tissue_type_id')
  
  const [activeTab, setActiveTab] = useState('wizard')
  const { creating, error, importExperiment } = useExperimentCreation()
  const [importFile, setImportFile] = useState<File | null>(null)
  
  const handleImport = async (e: React.FormEvent) => {
    e.preventDefault()
    if (!importFile) return
    
    try {
      const importedExperiment = await importExperiment(importFile)
      if (importedExperiment) {
        router.push(`/experiments/${importedExperiment.id}`)
      }
    } catch (err) {
      console.error('Failed to import experiment:', err)
    }
  }
  
  const handleFileChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    if (e.target.files && e.target.files.length > 0) {
      setImportFile(e.target.files[0])
    }
  }
  
  return (
    <div className="container mx-auto px-4 py-8">
      <div className="mb-8">
        <div className="flex items-center gap-2 mb-2">
          <Link href="/experiments" passHref>
            <Button variant="ghost" size="sm" className="h-8 w-8 p-0">
              <ArrowLeft className="h-4 w-4" />
            </Button>
          </Link>
          <h1 className="text-3xl font-bold">Create New Experiment</h1>
        </div>
        <p className="text-muted-foreground">
          Design and set up a new cryopreservation experiment
        </p>
      </div>
      
      <Tabs defaultValue="wizard" value={activeTab} onValueChange={setActiveTab}>
        <TabsList className="grid w-full grid-cols-2 mb-8">
          <TabsTrigger value="wizard" className="flex items-center gap-2">
            <ClipboardList className="h-4 w-4" />
            <span>Experiment Wizard</span>
          </TabsTrigger>
          <TabsTrigger value="import" className="flex items-center gap-2">
            <Upload className="h-4 w-4" />
            <span>Import Experiment</span>
          </TabsTrigger>
        </TabsList>
        
        <TabsContent value="wizard">
          <ExperimentCreationWizard
            initialProtocolId={protocolId as any}
            initialTissueTypeId={tissueTypeId as any}
          />
        </TabsContent>
        
        <TabsContent value="import">
          <Card>
            <CardHeader>
              <CardTitle>Import Experiment</CardTitle>
            </CardHeader>
            <CardContent>
              <form onSubmit={handleImport} className="space-y-4">
                <div className="space-y-2">
                  <label htmlFor="importFile" className="text-sm font-medium">
                    Experiment File (JSON)
                  </label>
                  <Input 
                    id="importFile" 
                    type="file" 
                    accept=".json" 
                    onChange={handleFileChange}
                    required
                  />
                  <p className="text-sm text-muted-foreground">
                    Upload a previously exported experiment file (.json format)
                  </p>
                </div>
                
                {error && (
                  <div className="bg-destructive/10 text-destructive p-3 rounded-md text-sm">
                    Error: {error.message}
                  </div>
                )}
                
                <div className="flex gap-2 justify-end mt-4">
                  <Button 
                    type="button" 
                    variant="outline" 
                    onClick={() => router.push('/experiments')}
                  >
                    Cancel
                  </Button>
                  <Button 
                    type="submit" 
                    disabled={creating || !importFile}
                  >
                    {creating ? 'Importing...' : 'Import Experiment'}
                  </Button>
                </div>
              </form>
            </CardContent>
          </Card>
        </TabsContent>
      </Tabs>
    </div>
  )
}
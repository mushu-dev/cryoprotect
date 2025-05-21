'use client'

import { useState } from 'react'
import { useForm } from 'react-hook-form'
import { z } from 'zod'
import { zodResolver } from '@hookform/resolvers/zod'
import { useImportFromPubChem } from '@/features/molecules/hooks/use-molecules'
import { Card, CardContent, CardDescription, CardFooter, CardHeader, CardTitle } from '@/components/ui/card'
import { Input } from '@/components/ui/input'
import { Button } from '@/components/ui/button'
import { Label } from '@/components/ui/label'
import { ArrowLeft, Database, Import } from 'lucide-react'
import { useRouter } from 'next/navigation'
import Link from 'next/link'
import { Alert, AlertDescription, AlertTitle } from '@/components/ui/alert'

// Validation schema
const importFormSchema = z.object({
  pubchem_cid: z.string().min(1, 'PubChem CID is required')
    .regex(/^\d+$/, 'PubChem CID must be a number')
})

type ImportFormValues = z.infer<typeof importFormSchema>

export default function ImportMoleculePage() {
  const router = useRouter()
  const [error, setError] = useState<string | null>(null)
  
  // Import mutation
  const { mutate: importMolecule, isPending, isError, isSuccess, data } = useImportFromPubChem()
  
  // Form setup
  const form = useForm<ImportFormValues>({
    resolver: zodResolver(importFormSchema),
    defaultValues: {
      pubchem_cid: ''
    }
  })
  
  // Submit handler
  const onSubmit = (values: ImportFormValues) => {
    setError(null)
    importMolecule(values.pubchem_cid, {
      onSuccess: (data) => {
        // Redirect to the molecule detail page
        router.push(`/molecules/${data.id}`)
      },
      onError: (error) => {
        setError(error instanceof Error 
          ? error.message 
          : 'Failed to import molecule. Please try again.'
        )
      }
    })
  }
  
  return (
    <div className="max-w-2xl mx-auto py-6">
      <div className="flex items-center mb-8">
        <Button 
          variant="outline" 
          size="icon" 
          className="mr-4"
          asChild
        >
          <Link href="/molecules">
            <ArrowLeft className="h-4 w-4" />
          </Link>
        </Button>
        <div>
          <h1 className="text-3xl font-bold tracking-tight">Import Molecule</h1>
          <p className="text-muted-foreground mt-1">
            Import a molecule from PubChem by its Compound ID (CID)
          </p>
        </div>
      </div>
      
      <Card>
        <CardHeader>
          <CardTitle>Import from PubChem</CardTitle>
          <CardDescription>
            Enter a PubChem Compound ID (CID) to import the molecule into the CryoProtect database.
          </CardDescription>
        </CardHeader>
        <CardContent>
          {error && (
            <Alert variant="destructive" className="mb-4">
              <AlertTitle>Import Failed</AlertTitle>
              <AlertDescription>{error}</AlertDescription>
            </Alert>
          )}
          
          <form onSubmit={form.handleSubmit(onSubmit)} className="space-y-4">
            <div className="space-y-2">
              <Label htmlFor="pubchem_cid">PubChem CID</Label>
              <div className="flex">
                <Input
                  id="pubchem_cid"
                  placeholder="e.g. 753"
                  {...form.register('pubchem_cid')}
                  className="flex-grow"
                />
                <Button 
                  type="submit"
                  className="ml-2"
                  disabled={isPending}
                >
                  {isPending ? 'Importing...' : 'Import'}
                </Button>
              </div>
              {form.formState.errors.pubchem_cid && (
                <p className="text-destructive text-sm mt-1">
                  {form.formState.errors.pubchem_cid.message}
                </p>
              )}
            </div>
          </form>
          
          <div className="mt-8 bg-muted/40 p-4 rounded-lg">
            <h3 className="font-medium mb-2">How to find a PubChem CID</h3>
            <ol className="list-decimal list-inside space-y-2 text-sm">
              <li>Visit the <a href="https://pubchem.ncbi.nlm.nih.gov/" target="_blank" rel="noopener noreferrer" className="text-primary hover:underline">PubChem website</a></li>
              <li>Search for a compound by name (e.g., "glycerol")</li>
              <li>Click on the compound in the search results</li>
              <li>Look for the "CID" number in the compound information</li>
              <li>Enter this number in the field above</li>
            </ol>
            <p className="text-sm text-muted-foreground mt-4">
              Example CIDs: 753 (Glycerol), 887 (DMSO), 174 (Ethylene glycol)
            </p>
          </div>
        </CardContent>
        <CardFooter className="flex justify-between">
          <Button variant="outline" asChild>
            <Link href="/molecules">
              Cancel
            </Link>
          </Button>
          <div className="text-sm text-muted-foreground">
            <Link href="https://pubchem.ncbi.nlm.nih.gov/" className="text-primary hover:underline" target="_blank" rel="noopener noreferrer">
              Visit PubChem
              <ExternalLink className="h-3 w-3 inline ml-1" />
            </Link>
          </div>
        </CardFooter>
      </Card>
    </div>
  )
}

function ExternalLink(props: React.SVGProps<SVGSVGElement>) {
  return (
    <svg
      {...props}
      xmlns="http://www.w3.org/2000/svg"
      width="24"
      height="24"
      viewBox="0 0 24 24"
      fill="none"
      stroke="currentColor"
      strokeWidth="2"
      strokeLinecap="round"
      strokeLinejoin="round"
    >
      <path d="M18 13v6a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2V8a2 2 0 0 1 2-2h6" />
      <polyline points="15 3 21 3 21 9" />
      <line x1="10" x2="21" y1="14" y2="3" />
    </svg>
  )
}
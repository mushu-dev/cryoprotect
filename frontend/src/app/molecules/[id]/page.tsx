'use client'

import { useParams } from 'next/navigation'
import { useMolecule, useMoleculeProperties } from '@/features/molecules/hooks/use-molecules'
import { MoleculeViewer3D } from '@/features/molecules/components/molecule-viewer-3d'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { Badge } from '@/components/ui/badge'
import { Button } from '@/components/ui/button'
import { 
  Tabs, 
  TabsContent, 
  TabsList, 
  TabsTrigger 
} from '@/components/ui/tabs'
import { 
  Table, 
  TableBody, 
  TableCell, 
  TableHead, 
  TableHeader, 
  TableRow 
} from '@/components/ui/table'
import Link from 'next/link'
import { ArrowLeft, FileText, ExternalLink, Beaker } from 'lucide-react'

export default function MoleculeDetailPage() {
  const params = useParams()
  const id = params.id as string
  
  const { data: molecule, isLoading, isError } = useMolecule(id)
  const { data: properties, isLoading: isLoadingProps } = useMoleculeProperties(id)
  
  if (isLoading) {
    return (
      <div className="text-center py-16">
        <p>Loading molecule information...</p>
      </div>
    )
  }
  
  if (isError || !molecule) {
    return (
      <div className="text-center py-16">
        <p className="text-destructive">Error loading molecule information</p>
        <Button asChild className="mt-4">
          <Link href="/molecules">Back to Molecules</Link>
        </Button>
      </div>
    )
  }
  
  return (
    <div className="space-y-6">
      <div className="flex flex-col space-y-2 sm:flex-row sm:justify-between sm:space-y-0">
        <div className="flex items-center space-x-2">
          <Button 
            variant="outline" 
            size="icon" 
            asChild
          >
            <Link href="/molecules">
              <ArrowLeft className="h-4 w-4" />
            </Link>
          </Button>
          <h1 className="text-2xl sm:text-3xl font-bold tracking-tight">
            {molecule.name || 'Unnamed Molecule'}
          </h1>
          {molecule.is_cryoprotectant && (
            <Badge variant="outline" className="ml-2 bg-primary/10 border-primary/20 text-primary">
              Cryoprotectant
            </Badge>
          )}
        </div>
        
        {molecule.pubchem_cid && (
          <Button variant="outline" size="sm" asChild>
            <a 
              href={`https://pubchem.ncbi.nlm.nih.gov/compound/${molecule.pubchem_cid}`}
              target="_blank"
              rel="noopener noreferrer"
            >
              <ExternalLink className="h-4 w-4 mr-2" />
              View on PubChem
            </a>
          </Button>
        )}
      </div>
      
      <div className="grid grid-cols-1 md:grid-cols-3 gap-6">
        <div className="md:col-span-2">
          <Card>
            <CardHeader className="pb-2">
              <CardTitle>Molecule Viewer</CardTitle>
            </CardHeader>
            <CardContent className="-mx-6">
              <MoleculeViewer3D 
                smiles={molecule.smiles}
                name={molecule.name}
                height={400}
                spin={true}
                style="stick"
              />
            </CardContent>
          </Card>
        </div>
        
        <div>
          <Card className="h-full">
            <CardHeader className="pb-2">
              <CardTitle>Molecule Information</CardTitle>
            </CardHeader>
            <CardContent>
              <dl className="space-y-4">
                {molecule.formula && (
                  <div>
                    <dt className="text-sm font-medium text-muted-foreground">Formula</dt>
                    <dd className="font-mono">{molecule.formula}</dd>
                  </div>
                )}
                
                {molecule.molecular_weight && (
                  <div>
                    <dt className="text-sm font-medium text-muted-foreground">Molecular Weight</dt>
                    <dd>{molecule.molecular_weight.toFixed(2)} g/mol</dd>
                  </div>
                )}
                
                <div>
                  <dt className="text-sm font-medium text-muted-foreground">SMILES</dt>
                  <dd className="font-mono text-xs break-all">{molecule.smiles}</dd>
                </div>
                
                {molecule.pubchem_cid && (
                  <div>
                    <dt className="text-sm font-medium text-muted-foreground">PubChem CID</dt>
                    <dd>{molecule.pubchem_cid}</dd>
                  </div>
                )}
                
                <div>
                  <dt className="text-sm font-medium text-muted-foreground">Cryoprotectant</dt>
                  <dd>{molecule.is_cryoprotectant ? 'Yes' : 'No'}</dd>
                </div>
              </dl>
            </CardContent>
          </Card>
        </div>
      </div>
      
      <Tabs defaultValue="properties">
        <TabsList>
          <TabsTrigger value="properties">
            <FileText className="h-4 w-4 mr-2" />
            Properties
          </TabsTrigger>
          <TabsTrigger value="mixtures">
            <Beaker className="h-4 w-4 mr-2" />
            Related Mixtures
          </TabsTrigger>
        </TabsList>
        
        <TabsContent value="properties" className="mt-6">
          <Card>
            <CardHeader className="pb-2">
              <CardTitle>Molecular Properties</CardTitle>
            </CardHeader>
            <CardContent>
              {isLoadingProps && <p>Loading properties...</p>}
              
              {!isLoadingProps && (!properties || properties.length === 0) && (
                <p className="text-muted-foreground">No properties available for this molecule.</p>
              )}
              
              {properties && properties.length > 0 && (
                <Table>
                  <TableHeader>
                    <TableRow>
                      <TableHead>Property</TableHead>
                      <TableHead>Value</TableHead>
                      <TableHead>Unit</TableHead>
                      <TableHead>Calculation Method</TableHead>
                    </TableRow>
                  </TableHeader>
                  <TableBody>
                    {properties.map((property) => (
                      <TableRow key={property.id}>
                        <TableCell className="font-medium">{property.property_name}</TableCell>
                        <TableCell>
                          {typeof property.property_value === 'number'
                            ? property.property_value.toFixed(4)
                            : property.property_value}
                        </TableCell>
                        <TableCell>{property.property_unit || '-'}</TableCell>
                        <TableCell>{property.calculation_method || 'Not specified'}</TableCell>
                      </TableRow>
                    ))}
                  </TableBody>
                </Table>
              )}
            </CardContent>
          </Card>
        </TabsContent>
        
        <TabsContent value="mixtures" className="mt-6">
          <Card>
            <CardHeader className="pb-2">
              <CardTitle>Related Mixtures</CardTitle>
            </CardHeader>
            <CardContent>
              <p className="text-muted-foreground">
                Mixtures containing this molecule will appear here.
              </p>
              
              {/* Mixtures will be implemented in a future update */}
              <div className="text-center py-8">
                <p className="text-muted-foreground">
                  Mixture data is being loaded...
                </p>
                <Button 
                  variant="outline" 
                  className="mt-4"
                  asChild
                >
                  <Link href="/mixtures">
                    <Beaker className="h-4 w-4 mr-2" />
                    Browse All Mixtures
                  </Link>
                </Button>
              </div>
            </CardContent>
          </Card>
        </TabsContent>
      </Tabs>
    </div>
  )
}
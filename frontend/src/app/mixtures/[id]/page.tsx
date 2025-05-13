'use client'

import { useParams } from 'next/navigation'
import { useMixture, useCryoprotectionScore } from '@/features/mixtures/hooks/use-mixtures'
import { MixtureCompositionChart } from '@/features/mixtures/components/mixture-composition-chart'
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
import { ArrowLeft, Edit, Trash, Plus, ExternalLink } from 'lucide-react'
import { ChartBar } from '@/components/ui/icons'
import { useState } from 'react'

export default function MixtureDetailPage() {
  const params = useParams()
  const id = params.id as string
  
  const { data: mixture, isLoading, isError } = useMixture(id)
  const { data: scoreData, isLoading: isLoadingScore } = useCryoprotectionScore(id)
  
  const [selectedMoleculeIndex, setSelectedMoleculeIndex] = useState<number | null>(null)

  // When a user clicks on a component in the table, show its 3D structure
  const handleSelectComponent = (index: number) => {
    setSelectedMoleculeIndex(index === selectedMoleculeIndex ? null : index)
  }
  
  if (isLoading) {
    return (
      <div className="text-center py-16">
        <p>Loading mixture information...</p>
      </div>
    )
  }
  
  if (isError || !mixture) {
    return (
      <div className="text-center py-16">
        <p className="text-destructive">Error loading mixture information</p>
        <Button asChild className="mt-4">
          <Link href="/mixtures">Back to Mixtures</Link>
        </Button>
      </div>
    )
  }
  
  // Prepare components data for the chart
  const chartComponents = mixture.components?.map(component => ({
    name: component.molecule?.name || 'Unknown Molecule',
    concentration: component.concentration,
    concentration_unit: component.concentration_unit || '%',
  })) || []
  
  // Get the currently selected molecule
  const selectedMolecule = selectedMoleculeIndex !== null && mixture.components && 
    mixture.components[selectedMoleculeIndex]?.molecule

  return (
    <div className="space-y-6">
      <div className="flex flex-col space-y-2 sm:flex-row sm:justify-between sm:space-y-0">
        <div className="flex items-center space-x-2">
          <Button 
            variant="outline" 
            size="icon" 
            asChild
          >
            <Link href="/mixtures">
              <ArrowLeft className="h-4 w-4" />
            </Link>
          </Button>
          <h1 className="text-2xl sm:text-3xl font-bold tracking-tight">
            {mixture.name || 'Unnamed Mixture'}
          </h1>
          {mixture.is_cryoprotectant_mixture && (
            <Badge variant="outline" className="ml-2 bg-primary/10 border-primary/20 text-primary">
              Cryoprotectant
            </Badge>
          )}
        </div>
        
        <div className="flex items-center space-x-2">
          <Button variant="outline" size="sm" asChild>
            <Link href={`/mixtures/${mixture.id}/edit`}>
              <Edit className="h-4 w-4 mr-2" />
              Edit
            </Link>
          </Button>
          <Button variant="destructive" size="sm">
            <Trash className="h-4 w-4 mr-2" />
            Delete
          </Button>
        </div>
      </div>
      
      {mixture.description && (
        <p className="text-muted-foreground max-w-2xl">
          {mixture.description}
        </p>
      )}
      
      <div className="grid grid-cols-1 md:grid-cols-3 gap-6">
        <div className="md:col-span-2">
          <Card>
            <CardHeader className="pb-2">
              <CardTitle>Mixture Composition</CardTitle>
            </CardHeader>
            <CardContent>
              {chartComponents.length > 0 ? (
                <MixtureCompositionChart 
                  components={chartComponents}
                  height={300}
                  showLegend={true}
                />
              ) : (
                <div className="flex items-center justify-center h-[300px] bg-muted/30 rounded-md">
                  <p className="text-muted-foreground">No components in this mixture</p>
                </div>
              )}
            </CardContent>
          </Card>
        </div>
        
        <div>
          <Card className="h-full">
            <CardHeader className="pb-2">
              <CardTitle>Mixture Information</CardTitle>
            </CardHeader>
            <CardContent>
              <dl className="space-y-4">
                <div>
                  <dt className="text-sm font-medium text-muted-foreground">Components</dt>
                  <dd>{mixture.components?.length || 0}</dd>
                </div>
                
                {mixture.total_concentration !== undefined && (
                  <div>
                    <dt className="text-sm font-medium text-muted-foreground">Total Concentration</dt>
                    <dd>{mixture.total_concentration}%</dd>
                  </div>
                )}
                
                {(mixture.cryoprotection_score !== undefined || scoreData?.score !== undefined) && (
                  <div>
                    <dt className="text-sm font-medium text-muted-foreground">Cryoprotection Score</dt>
                    <dd>{(mixture.cryoprotection_score || scoreData?.score || 0).toFixed(2)}</dd>
                  </div>
                )}
                
                {mixture.created_at && (
                  <div>
                    <dt className="text-sm font-medium text-muted-foreground">Created</dt>
                    <dd>{new Date(mixture.created_at).toLocaleDateString()}</dd>
                  </div>
                )}
                
                {mixture.updated_at && (
                  <div>
                    <dt className="text-sm font-medium text-muted-foreground">Last Updated</dt>
                    <dd>{new Date(mixture.updated_at).toLocaleDateString()}</dd>
                  </div>
                )}
              </dl>
            </CardContent>
          </Card>
        </div>
      </div>
      
      <Tabs defaultValue="components">
        <TabsList>
          <TabsTrigger value="components">
            Components
          </TabsTrigger>
          <TabsTrigger value="visualization">
            Visualization
          </TabsTrigger>
          <TabsTrigger value="analytics">
            <ChartBar className="h-4 w-4 mr-2" />
            Analytics
          </TabsTrigger>
        </TabsList>
        
        <TabsContent value="components" className="mt-6">
          <Card>
            <CardHeader className="pb-2 flex flex-row items-center justify-between">
              <CardTitle>Mixture Components</CardTitle>
              <Button size="sm" asChild>
                <Link href={`/mixtures/${mixture.id}/add-component`}>
                  <Plus className="h-4 w-4 mr-2" />
                  Add Component
                </Link>
              </Button>
            </CardHeader>
            <CardContent>
              {!mixture.components || mixture.components.length === 0 ? (
                <div className="text-center py-8">
                  <p className="text-muted-foreground">This mixture has no components.</p>
                  <Button 
                    variant="outline" 
                    className="mt-4"
                    asChild
                  >
                    <Link href={`/mixtures/${mixture.id}/add-component`}>
                      <Plus className="h-4 w-4 mr-2" />
                      Add Component
                    </Link>
                  </Button>
                </div>
              ) : (
                <Table>
                  <TableHeader>
                    <TableRow>
                      <TableHead>Molecule</TableHead>
                      <TableHead>Concentration</TableHead>
                      <TableHead>Formula</TableHead>
                      <TableHead className="text-right">Actions</TableHead>
                    </TableRow>
                  </TableHeader>
                  <TableBody>
                    {mixture.components.map((component, index) => (
                      <TableRow 
                        key={component.id} 
                        className={selectedMoleculeIndex === index ? "bg-muted" : ""}
                        onClick={() => handleSelectComponent(index)}
                      >
                        <TableCell className="font-medium">
                          {component.molecule?.name || 'Unknown Molecule'}
                        </TableCell>
                        <TableCell>
                          {component.concentration}{component.concentration_unit || '%'}
                        </TableCell>
                        <TableCell className="font-mono text-xs">
                          {component.molecule?.formula || '-'}
                        </TableCell>
                        <TableCell className="text-right">
                          <div className="flex justify-end gap-2">
                            <Button variant="ghost" size="icon" asChild>
                              <Link href={component.molecule ? `/molecules/${component.molecule.id}` : '#'}>
                                <ExternalLink className="h-4 w-4" />
                              </Link>
                            </Button>
                            <Button variant="ghost" size="icon">
                              <Edit className="h-4 w-4" />
                            </Button>
                            <Button variant="ghost" size="icon" className="text-destructive">
                              <Trash className="h-4 w-4" />
                            </Button>
                          </div>
                        </TableCell>
                      </TableRow>
                    ))}
                  </TableBody>
                </Table>
              )}
            </CardContent>
          </Card>
        </TabsContent>
        
        <TabsContent value="visualization" className="mt-6">
          <Card>
            <CardHeader className="pb-2">
              <CardTitle>Molecular Visualization</CardTitle>
            </CardHeader>
            <CardContent>
              {selectedMolecule ? (
                <div>
                  <div className="mb-4">
                    <h3 className="font-medium">{selectedMolecule.name}</h3>
                    <p className="text-sm text-muted-foreground">
                      Click on a component in the Components tab to view its 3D structure.
                    </p>
                  </div>
                  <MoleculeViewer3D 
                    smiles={selectedMolecule.smiles}
                    name={selectedMolecule.name}
                    height={400}
                    spin={true}
                  />
                </div>
              ) : (
                <div className="text-center py-16">
                  <p className="text-muted-foreground">
                    Select a component from the Components tab to view its 3D structure.
                  </p>
                </div>
              )}
            </CardContent>
          </Card>
        </TabsContent>
        
        <TabsContent value="analytics" className="mt-6">
          <Card>
            <CardHeader className="pb-2">
              <CardTitle>Mixture Analytics</CardTitle>
            </CardHeader>
            <CardContent>
              <div className="space-y-6">
                <div>
                  <h3 className="text-lg font-medium mb-2">Cryoprotection Analysis</h3>
                  {isLoadingScore ? (
                    <p>Loading cryoprotection score...</p>
                  ) : (
                    <div className="space-y-4">
                      <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
                        <Card>
                          <CardContent className="pt-6">
                            <div className="text-center">
                              <p className="text-sm text-muted-foreground">Overall Score</p>
                              <p className="text-3xl font-bold">
                                {(scoreData?.score || mixture.cryoprotection_score || 0).toFixed(2)}
                              </p>
                            </div>
                          </CardContent>
                        </Card>
                        
                        <Card>
                          <CardContent className="pt-6">
                            <div className="text-center">
                              <p className="text-sm text-muted-foreground">Components</p>
                              <p className="text-3xl font-bold">{mixture.components?.length || 0}</p>
                            </div>
                          </CardContent>
                        </Card>
                        
                        <Card>
                          <CardContent className="pt-6">
                            <div className="text-center">
                              <p className="text-sm text-muted-foreground">Total Concentration</p>
                              <p className="text-3xl font-bold">{mixture.total_concentration || 0}%</p>
                            </div>
                          </CardContent>
                        </Card>
                      </div>
                      
                      <div>
                        <h4 className="font-medium mb-2">Cryoprotection Factors</h4>
                        <p className="text-sm text-muted-foreground">
                          Detailed analytics about cryoprotection properties will appear here once the API provides this data.
                        </p>
                      </div>
                    </div>
                  )}
                </div>
                
                <div>
                  <h3 className="text-lg font-medium mb-2">Property Analysis</h3>
                  <p className="text-sm text-muted-foreground">
                    Visualization of properties across mixture components will appear here
                    once the API provides this data.
                  </p>
                </div>
              </div>
            </CardContent>
          </Card>
        </TabsContent>
      </Tabs>
    </div>
  )
}
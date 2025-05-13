import { DashboardStats } from '@/features/dashboard/components/dashboard-stats'
import { MoleculeViewer3D } from '@/features/molecules/components/molecule-viewer-3d'
import { MixtureCompositionChart } from '@/features/mixtures/components/mixture-composition-chart'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { Button } from '@/components/ui/button'
import Link from 'next/link'

export default function Home() {
  return (
    <div className="space-y-8">
      <div className="flex flex-col space-y-4">
        <h1 className="text-3xl font-bold tracking-tight">Welcome to CryoProtect</h1>
        <p className="text-muted-foreground">
          A comprehensive platform for analyzing and optimizing cryoprotectant molecules and mixtures
        </p>
      </div>
      
      <DashboardStats />
      
      <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
        <Card>
          <CardHeader>
            <CardTitle>Molecular Viewer</CardTitle>
          </CardHeader>
          <CardContent>
            <MoleculeViewer3D 
              smiles="C(C(CO)O)O" 
              height={300} 
              name="Glycerol"
            />
            <div className="mt-4 text-center">
              <Button asChild>
                <Link href="/molecules">Browse Molecules</Link>
              </Button>
            </div>
          </CardContent>
        </Card>
        
        <Card>
          <CardHeader>
            <CardTitle>Mixture Composition</CardTitle>
          </CardHeader>
          <CardContent>
            <MixtureCompositionChart 
              components={[
                { name: 'Glycerol', concentration: 40 },
                { name: 'DMSO', concentration: 30 },
                { name: 'Ethylene Glycol', concentration: 20 },
                { name: 'Trehalose', concentration: 10 }
              ]}
            />
            <div className="mt-4 text-center">
              <Button asChild>
                <Link href="/mixtures">Explore Mixtures</Link>
              </Button>
            </div>
          </CardContent>
        </Card>
      </div>
      
      <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
        <Card>
          <CardHeader>
            <CardTitle>Import Molecule</CardTitle>
          </CardHeader>
          <CardContent>
            <p className="mb-4 text-muted-foreground">Import a molecule from PubChem by its Compound ID (CID).</p>
            <Button asChild className="w-full">
              <Link href="/molecules/import">Import Molecule</Link>
            </Button>
          </CardContent>
        </Card>
        
        <Card>
          <CardHeader>
            <CardTitle>Create Mixture</CardTitle>
          </CardHeader>
          <CardContent>
            <p className="mb-4 text-muted-foreground">Create a new mixture by combining multiple molecules.</p>
            <Button asChild className="w-full">
              <Link href="/mixtures/new">Create Mixture</Link>
            </Button>
          </CardContent>
        </Card>
        
        <Card>
          <CardHeader>
            <CardTitle>Property Explorer</CardTitle>
          </CardHeader>
          <CardContent>
            <p className="mb-4 text-muted-foreground">Explore molecular properties and find patterns.</p>
            <Button asChild className="w-full">
              <Link href="/properties">Explore Properties</Link>
            </Button>
          </CardContent>
        </Card>
      </div>
    </div>
  )
}
import { MoleculesList } from '@/features/molecules/components/molecules-list'
import { Button } from '@/components/ui/button'
import Link from 'next/link'
import { Plus } from 'lucide-react'

export const metadata = {
  title: 'Molecules | CryoProtect',
  description: 'Browse and search cryoprotectant molecules in the CryoProtect database',
}

export default function MoleculesPage() {
  return (
    <div className="space-y-6">
      <div className="flex flex-col sm:flex-row justify-between gap-4">
        <div>
          <h1 className="text-3xl font-bold tracking-tight">Molecules</h1>
          <p className="text-muted-foreground mt-1">
            Browse and search cryoprotectant molecules in the database
          </p>
        </div>
        <div className="flex gap-2">
          <Button asChild>
            <Link href="/molecules/import">
              <Plus className="h-4 w-4 mr-2" />
              Import Molecule
            </Link>
          </Button>
        </div>
      </div>
      
      <MoleculesList />
    </div>
  )
}
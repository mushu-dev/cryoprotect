import { MixturesList } from '@/features/mixtures/components/mixtures-list'
import { Button } from '@/components/ui/button'
import Link from 'next/link'
import { Plus } from 'lucide-react'

export const metadata = {
  title: 'Mixtures | CryoProtect',
  description: 'Browse and analyze cryoprotectant mixtures in the CryoProtect database',
}

export default function MixturesPage() {
  return (
    <div className="space-y-6">
      <div className="flex flex-col sm:flex-row justify-between gap-4">
        <div>
          <h1 className="text-3xl font-bold tracking-tight">Mixtures</h1>
          <p className="text-muted-foreground mt-1">
            Browse and analyze cryoprotectant mixtures in the database
          </p>
        </div>
        <div className="flex gap-2">
          <Button asChild>
            <Link href="/mixtures/new">
              <Plus className="h-4 w-4 mr-2" />
              Create Mixture
            </Link>
          </Button>
        </div>
      </div>
      
      <MixturesList />
    </div>
  )
}
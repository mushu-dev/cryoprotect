'use client'

import { ConvexMoleculesList } from '@/features/molecules/components/convex-molecules-list'
import { Button } from '@/components/ui/button'
import { Card } from '@/components/ui/card'
import Link from 'next/link'

export default function ConvexTestPage() {
  const isConvexEnabled = process.env.NEXT_PUBLIC_USE_CONVEX === 'true'

  if (!isConvexEnabled) {
    return (
      <div className="container mx-auto p-6">
        <Card className="p-6 max-w-3xl mx-auto">
          <h1 className="text-2xl font-bold mb-4">Convex Integration Test</h1>
          <div className="p-4 bg-amber-100 dark:bg-amber-900/30 rounded-md mb-6">
            <p className="text-amber-800 dark:text-amber-200">
              Convex integration is not enabled. Set <code>NEXT_PUBLIC_USE_CONVEX=true</code> in your environment to enable it.
            </p>
          </div>
          <Button asChild>
            <Link href="/">Back to Home</Link>
          </Button>
        </Card>
      </div>
    )
  }

  return (
    <div className="container mx-auto p-6">
      <Card className="p-6 max-w-3xl mx-auto">
        <h1 className="text-2xl font-bold mb-4">Convex Integration Test</h1>
        <p className="mb-6">
          This page demonstrates integration with the Convex backend. You can create and view molecules stored in the Convex database.
        </p>
        
        <div className="mb-8">
          <h2 className="text-xl font-semibold mb-3">Molecules</h2>
          <ConvexMoleculesList />
        </div>
        
        <div className="flex justify-between">
          <Button asChild variant="outline">
            <Link href="/">Back to Home</Link>
          </Button>
          <Button asChild>
            <Link href="/molecules">Go to Molecules</Link>
          </Button>
        </div>
      </Card>
    </div>
  )
}
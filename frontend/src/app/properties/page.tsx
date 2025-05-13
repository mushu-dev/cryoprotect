import { PropertyExplorer } from '@/features/properties/components/property-explorer'

export const metadata = {
  title: 'Property Explorer | CryoProtect',
  description: 'Explore and analyze molecular properties in the CryoProtect database',
}

export default function PropertiesPage() {
  return (
    <div className="space-y-6">
      <div>
        <h1 className="text-3xl font-bold tracking-tight">Property Explorer</h1>
        <p className="text-muted-foreground mt-1">
          Explore and analyze molecular properties in the database
        </p>
      </div>
      
      <PropertyExplorer />
    </div>
  )
}
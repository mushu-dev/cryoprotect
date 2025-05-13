'use client'

import { Card, CardContent, CardFooter, CardHeader, CardTitle } from '@/components/ui/card'
import { Button } from '@/components/ui/button'
import { MixtureCompositionChart } from './mixture-composition-chart'
import { type Mixture } from '../services/mixture-service'
import { Badge } from '@/components/ui/badge'
import Link from 'next/link'
import { truncateString } from '@/lib/utils'
import { formatConcentration } from '@/lib/utils'

interface MixtureCardProps {
  mixture: Mixture
  showChart?: boolean
  className?: string
}

export function MixtureCard({ 
  mixture, 
  showChart = true,
  className = '' 
}: MixtureCardProps) {
  // Prepare components data for the chart
  const chartComponents = mixture.components?.map(component => ({
    name: component.molecule?.name || 'Unknown Molecule',
    concentration: component.concentration,
    concentration_unit: component.concentration_unit || '%',
    color: undefined, // Let the chart assign colors
  })) || []
  
  return (
    <Card className={`overflow-hidden h-full flex flex-col ${className}`}>
      <CardHeader className="pb-2">
        <div className="flex justify-between items-start">
          <CardTitle className="text-lg font-semibold">
            {truncateString(mixture.name || 'Unnamed Mixture', 30)}
          </CardTitle>
          {mixture.is_cryoprotectant_mixture && (
            <Badge variant="outline" className="bg-primary/10 border-primary/20 text-primary">
              Cryoprotectant
            </Badge>
          )}
        </div>
        {mixture.description && (
          <p className="text-sm text-muted-foreground mt-1">
            {truncateString(mixture.description, 100)}
          </p>
        )}
      </CardHeader>
      
      <CardContent className="flex-grow">
        {showChart && chartComponents.length > 0 && (
          <MixtureCompositionChart
            components={chartComponents}
            height={200}
            title=""
          />
        )}
        
        <div className="mt-4 space-y-2">
          {/* Component counts */}
          <div className="flex justify-between text-sm">
            <span className="text-muted-foreground">Components:</span>
            <span>{mixture.components?.length || 0}</span>
          </div>
          
          {/* Total concentration if available */}
          {mixture.total_concentration !== undefined && (
            <div className="flex justify-between text-sm">
              <span className="text-muted-foreground">Total Concentration:</span>
              <span>{mixture.total_concentration}%</span>
            </div>
          )}
          
          {/* Cryoprotection score if available */}
          {mixture.cryoprotection_score !== undefined && (
            <div className="flex justify-between text-sm">
              <span className="text-muted-foreground">Cryoprotection Score:</span>
              <span>{mixture.cryoprotection_score.toFixed(2)}</span>
            </div>
          )}
        </div>
      </CardContent>
      
      <CardFooter className="pt-2">
        <Button asChild className="w-full" variant="outline">
          <Link href={`/mixtures/${mixture.id}`}>
            View Details
          </Link>
        </Button>
      </CardFooter>
    </Card>
  )
}
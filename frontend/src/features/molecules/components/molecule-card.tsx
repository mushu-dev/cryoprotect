'use client'

import { Card, CardContent, CardFooter, CardHeader, CardTitle } from '@/components/ui/card'
import { Button } from '@/components/ui/button'
import { MoleculeViewer3D } from './molecule-viewer-3d'
import { type Molecule } from '../services/molecule-service'
import { Badge } from '@/components/ui/badge'
import Link from 'next/link'
import { truncateString } from '@/lib/utils'

interface MoleculeCardProps {
  molecule: Molecule
  showViewer?: boolean
  className?: string
}

export function MoleculeCard({ 
  molecule, 
  showViewer = true,
  className = '' 
}: MoleculeCardProps) {
  return (
    <Card className={`overflow-hidden h-full flex flex-col ${className}`}>
      <CardHeader className="pb-2">
        <div className="flex justify-between items-start">
          <CardTitle className="text-lg font-semibold">
            {truncateString(molecule.name || 'Unnamed Molecule', 30)}
          </CardTitle>
          {molecule.is_cryoprotectant && (
            <Badge variant="outline" className="bg-primary/10 border-primary/20 text-primary">
              Cryoprotectant
            </Badge>
          )}
        </div>
        <div className="text-xs text-muted-foreground">
          {molecule.formula && <div className="font-mono">{molecule.formula}</div>}
          {molecule.pubchem_cid && (
            <div className="mt-1">
              PubChem CID: {molecule.pubchem_cid}
            </div>
          )}
        </div>
      </CardHeader>
      
      <CardContent className="flex-grow">
        {showViewer && molecule.smiles && (
          <div className="-mx-6">
            <MoleculeViewer3D
              smiles={molecule.smiles}
              name={molecule.name}
              height={200}
              spin={false}
            />
          </div>
        )}
        
        <div className="mt-4 space-y-2">
          {molecule.molecular_weight && (
            <div className="flex justify-between text-sm">
              <span className="text-muted-foreground">Molecular Weight:</span>
              <span>{molecule.molecular_weight.toFixed(2)} g/mol</span>
            </div>
          )}
          
          <div className="flex justify-between text-sm">
            <span className="text-muted-foreground">SMILES:</span>
            <span className="font-mono text-xs max-w-[180px] truncate">
              {molecule.smiles}
            </span>
          </div>
        </div>
      </CardContent>
      
      <CardFooter className="pt-2">
        <Button asChild className="w-full" variant="outline">
          <Link href={`/molecules/${molecule.id}`}>
            View Details
          </Link>
        </Button>
      </CardFooter>
    </Card>
  )
}
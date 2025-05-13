'use client'

import { useState } from 'react'
import { useMolecules } from '../hooks/use-molecules'
import { MoleculeCard } from './molecule-card'
import { Input } from '@/components/ui/input'
import { Button } from '@/components/ui/button'
import { 
  Select, 
  SelectContent, 
  SelectItem, 
  SelectTrigger, 
  SelectValue 
} from '@/components/ui/select'
import { 
  Pagination, 
  PaginationContent, 
  PaginationEllipsis, 
  PaginationItem, 
  PaginationLink, 
  PaginationNext, 
  PaginationPrevious 
} from '@/components/ui/pagination'
import { Switch } from '@/components/ui/switch'
import { Label } from '@/components/ui/label'
import { Search, Filter } from 'lucide-react'
import { type MoleculeParams } from '../services/molecule-service'

export function MoleculesList() {
  const [params, setParams] = useState<MoleculeParams>({
    page: 1,
    limit: 12,
    sort_by: 'name',
    sort_order: 'asc',
  })
  
  const [search, setSearch] = useState('')
  
  const { data, isLoading, isError, error } = useMolecules(params)
  
  // Handle search input
  const handleSearch = () => {
    setParams(prev => ({ ...prev, search, page: 1 }))
  }
  
  // Handle filter for cryoprotectants only
  const handleCryoprotectantsFilter = (checked: boolean) => {
    setParams(prev => ({ 
      ...prev, 
      is_cryoprotectant: checked || undefined,
      page: 1
    }))
  }
  
  // Handle sort change
  const handleSortChange = (value: string) => {
    const [sort_by, sort_order] = value.split(':')
    setParams(prev => ({ 
      ...prev, 
      sort_by, 
      sort_order: sort_order as 'asc' | 'desc' 
    }))
  }
  
  // Handle pagination
  const handlePageChange = (page: number) => {
    setParams(prev => ({ ...prev, page }))
  }
  
  // Calculate total pages
  const totalPages = data ? Math.ceil(data.total / data.limit) : 0
  
  // Generate pagination items
  const getPaginationItems = () => {
    const items = []
    const currentPage = params.page || 1
    
    // Always show first page
    items.push(
      <PaginationItem key="first">
        <PaginationLink 
          onClick={() => handlePageChange(1)}
          isActive={currentPage === 1}
        >
          1
        </PaginationLink>
      </PaginationItem>
    )
    
    // Show ellipsis if needed
    if (currentPage > 3) {
      items.push(
        <PaginationItem key="ellipsis-1">
          <PaginationEllipsis />
        </PaginationItem>
      )
    }
    
    // Show pages around current page
    for (
      let i = Math.max(2, currentPage - 1); 
      i <= Math.min(totalPages - 1, currentPage + 1); 
      i++
    ) {
      if (i <= 1 || i >= totalPages) continue
      
      items.push(
        <PaginationItem key={i}>
          <PaginationLink 
            onClick={() => handlePageChange(i)}
            isActive={currentPage === i}
          >
            {i}
          </PaginationLink>
        </PaginationItem>
      )
    }
    
    // Show ellipsis if needed
    if (currentPage < totalPages - 2) {
      items.push(
        <PaginationItem key="ellipsis-2">
          <PaginationEllipsis />
        </PaginationItem>
      )
    }
    
    // Always show last page
    if (totalPages > 1) {
      items.push(
        <PaginationItem key="last">
          <PaginationLink 
            onClick={() => handlePageChange(totalPages)}
            isActive={currentPage === totalPages}
          >
            {totalPages}
          </PaginationLink>
        </PaginationItem>
      )
    }
    
    return items
  }
  
  return (
    <div className="space-y-6">
      {/* Search and filter row */}
      <div className="flex flex-col sm:flex-row gap-4">
        <div className="relative flex-grow">
          <Search className="absolute left-2.5 top-2.5 h-4 w-4 text-muted-foreground" />
          <Input 
            placeholder="Search molecules..."
            className="pl-9"
            value={search}
            onChange={(e) => setSearch(e.target.value)}
            onKeyDown={(e) => {
              if (e.key === 'Enter') handleSearch()
            }}
          />
        </div>
        
        <div className="flex gap-2">
          <Button 
            variant="outline" 
            className="w-[90px]"
            onClick={handleSearch}
          >
            Search
          </Button>
          
          <Select defaultValue="name:asc" onValueChange={handleSortChange}>
            <SelectTrigger className="w-[165px]">
              <Filter className="h-4 w-4 mr-2" />
              <SelectValue placeholder="Sort by" />
            </SelectTrigger>
            <SelectContent>
              <SelectItem value="name:asc">Name (A-Z)</SelectItem>
              <SelectItem value="name:desc">Name (Z-A)</SelectItem>
              <SelectItem value="molecular_weight:asc">Weight (Low-High)</SelectItem>
              <SelectItem value="molecular_weight:desc">Weight (High-Low)</SelectItem>
              <SelectItem value="created_at:desc">Newest First</SelectItem>
              <SelectItem value="created_at:asc">Oldest First</SelectItem>
            </SelectContent>
          </Select>
        </div>
      </div>
      
      {/* Cryoprotectant filter */}
      <div className="flex items-center space-x-2">
        <Switch 
          id="cryoprotectant" 
          checked={!!params.is_cryoprotectant} 
          onCheckedChange={handleCryoprotectantsFilter}
        />
        <Label htmlFor="cryoprotectant">Show cryoprotectants only</Label>
      </div>
      
      {/* Loading state */}
      {isLoading && (
        <div className="text-center py-12">
          <p>Loading molecules...</p>
        </div>
      )}
      
      {/* Error state */}
      {isError && (
        <div className="text-center py-12 text-destructive">
          <p>Error loading molecules: {error instanceof Error ? error.message : 'Unknown error'}</p>
          <Button 
            variant="outline" 
            className="mt-4"
            onClick={() => setParams({ ...params })} // Retry
          >
            Retry
          </Button>
        </div>
      )}
      
      {/* Empty state */}
      {data?.molecules?.length === 0 && (
        <div className="text-center py-12">
          <p>No molecules found.</p>
          {(params.search || params.is_cryoprotectant) && (
            <Button 
              variant="outline" 
              className="mt-4"
              onClick={() => {
                setParams({
                  page: 1,
                  limit: 12,
                  sort_by: 'name',
                  sort_order: 'asc',
                })
                setSearch('')
              }}
            >
              Clear Filters
            </Button>
          )}
        </div>
      )}
      
      {/* Molecules grid */}
      {data?.molecules && data.molecules.length > 0 && (
        <div className="grid grid-cols-1 sm:grid-cols-2 md:grid-cols-3 lg:grid-cols-4 gap-6">
          {data.molecules.map((molecule) => (
            <MoleculeCard key={molecule.id} molecule={molecule} />
          ))}
        </div>
      )}
      
      {/* Pagination */}
      {data && data.total > data.limit && (
        <Pagination className="mt-8">
          <PaginationContent>
            <PaginationPrevious 
              onClick={() => handlePageChange(Math.max(1, (params.page || 1) - 1))}
            />
            
            {getPaginationItems()}
            
            <PaginationNext 
              onClick={() => handlePageChange(Math.min(totalPages, (params.page || 1) + 1))}
            />
          </PaginationContent>
        </Pagination>
      )}
    </div>
  )
}
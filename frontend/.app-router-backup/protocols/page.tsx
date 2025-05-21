'use client'

import { useState } from 'react'
import { Metadata } from 'next'
import Link from 'next/link'
import { useProtocols } from '@/features/protocols/hooks'
import { Protocol } from '@/features/experiments/services/experiment-service'
import { ProtocolBuilder } from '@/features/protocols/components/protocol-builder'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { Button } from '@/components/ui/button'
import { Input } from '@/components/ui/input'
import { Plus, Search, Filter, SortAsc, SortDesc, Clock, User, Tag } from 'lucide-react'
import { UUID } from 'crypto'
import {
  Dialog,
  DialogContent,
  DialogHeader,
  DialogTitle,
} from '@/components/ui/dialog'
import {
  DropdownMenu,
  DropdownMenuContent,
  DropdownMenuItem,
  DropdownMenuTrigger,
} from '@/components/ui/dropdown-menu'
import {
  Select,
  SelectContent,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from '@/components/ui/select'
import { Badge } from '@/components/ui/badge'
import { Pagination, PaginationContent, PaginationItem, PaginationLink, PaginationNext, PaginationPrevious } from '@/components/ui/pagination'
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs'
import { ProtocolCard } from '@/features/protocols/components/protocol-card'

// This metadata is moved to layout.tsx since this is now a client component
// export const metadata: Metadata = {
//   title: 'Protocols | CryoProtect',
//   description: 'Browse and manage standardized protocols for cryopreservation experiments',
// };

export default function ProtocolsPage() {
  const [selectedProtocolId, setSelectedProtocolId] = useState<UUID | null>(null)
  const [isCreating, setIsCreating] = useState(false)
  const [searchQuery, setSearchQuery] = useState('')
  const [sortField, setSortField] = useState<string>('updated_at')
  const [sortOrder, setSortOrder] = useState<'asc' | 'desc'>('desc')
  const [activeTab, setActiveTab] = useState('all')
  
  const {
    protocols,
    loading,
    error,
    totalCount,
    page,
    perPage,
    totalPages,
    changePage,
    updateParams,
    refreshProtocols
  } = useProtocols({
    page: 1,
    per_page: 10,
    sort_by: sortField,
    sort_order: sortOrder,
    search: searchQuery,
    ...(activeTab === 'templates' ? { is_template: true } : 
        activeTab === 'my' ? { is_template: false } : {})
  })
  
  const handleSearch = (e: React.FormEvent) => {
    e.preventDefault()
    updateParams({ search: searchQuery })
  }
  
  const handleSortChange = (field: string) => {
    if (field === sortField) {
      // Toggle sort order if clicking on the same field
      const newOrder = sortOrder === 'asc' ? 'desc' : 'asc'
      setSortOrder(newOrder)
      updateParams({ sort_order: newOrder })
    } else {
      // New field, default to descending
      setSortField(field)
      setSortOrder('desc')
      updateParams({ sort_by: field, sort_order: 'desc' })
    }
  }
  
  const handleProtocolSave = (protocol: Protocol) => {
    setSelectedProtocolId(null)
    setIsCreating(false)
    refreshProtocols()
  }
  
  const formatDate = (dateString: string) => {
    return new Date(dateString).toLocaleDateString('en-US', {
      year: 'numeric',
      month: 'short',
      day: 'numeric',
    })
  }
  
  const handleUseTemplate = (templateId: UUID) => {
    // Logic to use template to create a new protocol
    console.log('Using template:', templateId)
    // This would typically open the protocol builder with a template
  }
  
  return (
    <div className="container mx-auto px-4 py-8">
      <div className="flex flex-col md:flex-row justify-between items-start md:items-center mb-8">
        <div>
          <h1 className="text-3xl font-bold mb-2">Protocols</h1>
          <p className="text-muted-foreground">
            Standardized protocols for reproducible cryopreservation procedures
          </p>
        </div>
        
        <div className="mt-4 md:mt-0">
          <Button onClick={() => setIsCreating(true)}>
            <Plus className="h-4 w-4 mr-2" /> Create New Protocol
          </Button>
        </div>
      </div>
      
      <div className="mb-8">
        <form onSubmit={handleSearch} className="max-w-md mx-auto mb-6 flex gap-2">
          <div className="relative flex-grow">
            <Search className="absolute left-3 top-1/2 transform -translate-y-1/2 h-4 w-4 text-muted-foreground" />
            <Input 
              type="search" 
              placeholder="Search protocols..." 
              className="w-full pl-10 pr-4 focus-visible:ring-primary"
              value={searchQuery}
              onChange={(e) => setSearchQuery(e.target.value)}
            />
          </div>
          <Button 
            type="submit" 
            variant="secondary"
            className="shrink-0 hover:bg-secondary/90 transition-colors"
          >
            <span className="sr-only sm:not-sr-only sm:inline-block sm:mr-2">Search</span>
            <Search className="h-4 w-4 sm:hidden" />
          </Button>
        </form>

        <Tabs 
          defaultValue="all" 
          value={activeTab} 
          onValueChange={(value) => {
            setActiveTab(value)
            // Reset page when changing tabs
            changePage(1)
          }}
          className="w-full"
        >
          <TabsList className="mb-6 w-full sm:w-auto grid grid-cols-3 sm:inline-flex">
            <TabsTrigger value="all" className="text-sm flex items-center gap-1 sm:min-w-[120px]">
              <span className="hidden sm:inline">All</span> Protocols
            </TabsTrigger>
            <TabsTrigger value="my" className="text-sm flex items-center gap-1 sm:min-w-[120px]">
              <span className="hidden sm:inline">My</span> Protocols
            </TabsTrigger>
            <TabsTrigger value="templates" className="text-sm flex items-center gap-1 sm:min-w-[120px]">
              Templates
            </TabsTrigger>
          </TabsList>
          
          <TabsContent value="all">
            {renderProtocolGrid(protocols, loading, error)}
          </TabsContent>
          
          <TabsContent value="my">
            {protocols.length > 0 
              ? renderProtocolGrid(protocols, loading, error)
              : (
                <div className="text-center py-12 bg-card rounded-lg">
                  <p className="mb-4">You haven't created any protocols yet.</p>
                  <Button onClick={() => setIsCreating(true)}>
                    Create Your First Protocol
                  </Button>
                </div>
              )
            }
          </TabsContent>
          
          <TabsContent value="templates">
            {renderProtocolGrid(protocols, loading, error)}
          </TabsContent>
        </Tabs>
        
        {totalPages > 1 && (
          <div className="mt-8 flex flex-col items-center">
            <div className="text-sm text-muted-foreground mb-4">
              Showing <span className="font-medium">{protocols.length}</span> of <span className="font-medium">{totalCount}</span> protocols
            </div>
            <Pagination className="mt-2">
              <PaginationContent>
                <PaginationItem>
                  <PaginationPrevious 
                    href="#" 
                    onClick={(e) => {
                      e.preventDefault();
                      if (page > 1) changePage(page - 1);
                    }} 
                    className={`transition-colors ${page === 1 ? 'pointer-events-none opacity-50' : 'hover:bg-secondary/20'}`}
                  />
                </PaginationItem>
                
                {/* For small screens, simplified pagination */}
                <div className="sm:hidden flex items-center px-2">
                  <span className="text-sm font-medium">Page {page} of {totalPages}</span>
                </div>
                
                {/* For larger screens, full pagination */}
                <div className="hidden sm:flex">
                  {Array.from({ length: Math.min(5, totalPages) }).map((_, i) => {
                    // Show 5 pages max with current page in middle when possible
                    let pageNum = i + 1;
                    if (totalPages > 5) {
                      if (page > 3) {
                        pageNum = page - 3 + i;
                      }
                      if (page > totalPages - 2) {
                        pageNum = totalPages - 4 + i;
                      }
                    }
                    
                    return (
                      <PaginationItem key={i}>
                        <PaginationLink 
                          href="#" 
                          isActive={page === pageNum}
                          onClick={(e) => {
                            e.preventDefault();
                            changePage(pageNum);
                          }}
                          className="transition-colors hover:bg-secondary/20"
                        >
                          {pageNum}
                        </PaginationLink>
                      </PaginationItem>
                    );
                  })}
                </div>
                
                <PaginationItem>
                  <PaginationNext 
                    href="#" 
                    onClick={(e) => {
                      e.preventDefault();
                      if (page < totalPages) changePage(page + 1);
                    }}
                    className={`transition-colors ${page === totalPages ? 'pointer-events-none opacity-50' : 'hover:bg-secondary/20'}`}
                  />
                </PaginationItem>
              </PaginationContent>
            </Pagination>
          </div>
        )}
      </div>
      
      {/* Create Protocol Dialog */}
      <Dialog open={isCreating} onOpenChange={setIsCreating}>
        <DialogContent className="max-w-5xl max-h-[90vh] overflow-y-auto">
          <DialogHeader>
            <DialogTitle>Create New Protocol</DialogTitle>
          </DialogHeader>
          <ProtocolBuilder 
            onSave={handleProtocolSave}
            onCancel={() => setIsCreating(false)}
          />
        </DialogContent>
      </Dialog>
      
      {/* Edit Protocol Dialog */}
      <Dialog open={!!selectedProtocolId} onOpenChange={(open) => !open && setSelectedProtocolId(null)}>
        <DialogContent className="max-w-5xl max-h-[90vh] overflow-y-auto">
          <DialogHeader>
            <DialogTitle>Edit Protocol</DialogTitle>
          </DialogHeader>
          {selectedProtocolId && (
            <ProtocolBuilder 
              protocolId={selectedProtocolId}
              onSave={handleProtocolSave}
              onCancel={() => setSelectedProtocolId(null)}
            />
          )}
        </DialogContent>
      </Dialog>
    </div>
  )
  
  function renderProtocolGrid(protocols: Protocol[], loading: boolean, error: Error | null) {
    if (loading) {
      return <div className="text-center py-8">Loading protocols...</div>;
    }
    
    if (error) {
      return <div className="text-center py-8 text-red-500">Error loading protocols: {error.message}</div>;
    }
    
    if (protocols.length === 0) {
      return (
        <div className="text-center py-12 bg-card rounded-lg">
          <p className="mb-4">No protocols found matching your search criteria.</p>
          <Button onClick={() => {
            setSearchQuery('');
            updateParams({ search: '' });
          }}>Clear Search</Button>
        </div>
      );
    }
    
    return (
      <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 xl:grid-cols-4 gap-4 md:gap-6">
        {protocols.map(protocol => (
          <div 
            key={protocol.id.toString()} 
            onClick={() => setSelectedProtocolId(protocol.id)}
            className="cursor-pointer transition-transform hover:scale-[1.01] focus-within:scale-[1.01]"
          >
            <ProtocolCard 
              protocol={{
                ...protocol,
                step_count: protocol.steps?.length || 0,
                duration: protocol.steps?.reduce((total, step) => 
                  total + (step.duration || 0), 0)?.toString() || '',
                author: { name: protocol.created_by, id: protocol.created_by }
              }}
              onUseTemplate={protocol.is_template ? () => handleUseTemplate(protocol.id) : undefined}
            />
          </div>
        ))}
      </div>
    );
  }
}
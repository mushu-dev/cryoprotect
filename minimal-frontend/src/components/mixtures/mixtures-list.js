/**
 * MixturesList Component
 * Display a list of mixtures with search, filtering, and pagination
 */
import React, { useState } from 'react';
import { useMixtures } from '../../hooks/use-mixtures';
import { MixtureCard } from './mixture-card';
import { Input } from '../ui/input';
import { Button } from '../ui/button';
import { 
  Select, 
  SelectContent, 
  SelectItem, 
  SelectTrigger, 
  SelectValue 
} from '../ui/select';
import { 
  Pagination, 
  PaginationContent, 
  PaginationEllipsis, 
  PaginationItem, 
  PaginationLink, 
  PaginationNext, 
  PaginationPrevious 
} from '../ui/pagination';
import { Switch } from '../ui/switch';
import { Label } from '../ui/label';

// Icons
const SearchIcon = () => (
  <svg xmlns="http://www.w3.org/2000/svg" width="24" height="24" viewBox="0 0 24 24" fill="none" 
    stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" 
    className="h-4 w-4">
    <circle cx="11" cy="11" r="8"></circle>
    <path d="m21 21-4.3-4.3"></path>
  </svg>
);

const FilterIcon = () => (
  <svg xmlns="http://www.w3.org/2000/svg" width="24" height="24" viewBox="0 0 24 24" fill="none" 
    stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" 
    className="h-4 w-4 mr-2">
    <polygon points="22 3 2 3 10 12.46 10 19 14 21 14 12.46 22 3"></polygon>
  </svg>
);

export function MixturesList() {
  const [params, setParams] = useState({
    page: 1,
    limit: 12,
    sort_by: 'name',
    sort_order: 'asc',
  });
  
  const [search, setSearch] = useState('');
  
  const { data, isLoading, isError, error } = useMixtures(params);
  
  // Handle search input
  const handleSearch = () => {
    setParams(prev => ({ ...prev, search, page: 1 }));
  };
  
  // Handle filter for cryoprotectants only
  const handleCryoprotectantsFilter = (checked) => {
    setParams(prev => ({ 
      ...prev, 
      is_cryoprotectant_mixture: checked || undefined,
      page: 1
    }));
  };
  
  // Handle sort change
  const handleSortChange = (value) => {
    const [sort_by, sort_order] = value.split(':');
    setParams(prev => ({ 
      ...prev, 
      sort_by, 
      sort_order
    }));
  };
  
  // Handle pagination
  const handlePageChange = (page) => {
    setParams(prev => ({ ...prev, page }));
  };
  
  // Calculate total pages
  const totalPages = data ? Math.ceil(data.total / data.limit) : 0;
  
  // Generate pagination items
  const getPaginationItems = () => {
    const items = [];
    const currentPage = params.page || 1;
    
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
    );
    
    // Show ellipsis if needed
    if (currentPage > 3) {
      items.push(
        <PaginationItem key="ellipsis-1">
          <PaginationEllipsis />
        </PaginationItem>
      );
    }
    
    // Show pages around current page
    for (
      let i = Math.max(2, currentPage - 1); 
      i <= Math.min(totalPages - 1, currentPage + 1); 
      i++
    ) {
      if (i <= 1 || i >= totalPages) continue;
      
      items.push(
        <PaginationItem key={i}>
          <PaginationLink 
            onClick={() => handlePageChange(i)}
            isActive={currentPage === i}
          >
            {i}
          </PaginationLink>
        </PaginationItem>
      );
    }
    
    // Show ellipsis if needed
    if (currentPage < totalPages - 2) {
      items.push(
        <PaginationItem key="ellipsis-2">
          <PaginationEllipsis />
        </PaginationItem>
      );
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
      );
    }
    
    return items;
  };
  
  return (
    <div className="space-y-6">
      {/* Search and filter row */}
      <div className="flex flex-col sm:flex-row gap-4">
        <div className="relative flex-grow">
          <div className="absolute left-2.5 top-2.5 h-4 w-4 text-muted-foreground">
            <SearchIcon />
          </div>
          <Input 
            placeholder="Search mixtures..."
            className="pl-9"
            value={search}
            onChange={(e) => setSearch(e.target.value)}
            onKeyDown={(e) => {
              if (e.key === 'Enter') handleSearch();
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
              <FilterIcon />
              <SelectValue placeholder="Sort by" />
            </SelectTrigger>
            <SelectContent>
              <SelectItem value="name:asc">Name (A-Z)</SelectItem>
              <SelectItem value="name:desc">Name (Z-A)</SelectItem>
              <SelectItem value="cryoprotection_score:desc">Score (High-Low)</SelectItem>
              <SelectItem value="cryoprotection_score:asc">Score (Low-High)</SelectItem>
              <SelectItem value="created_at:desc">Newest First</SelectItem>
              <SelectItem value="created_at:asc">Oldest First</SelectItem>
            </SelectContent>
          </Select>
        </div>
      </div>
      
      {/* Cryoprotectant filter */}
      <div className="flex items-center space-x-2">
        <Switch 
          id="cryoprotectant-mixture" 
          checked={!!params.is_cryoprotectant_mixture} 
          onCheckedChange={handleCryoprotectantsFilter}
        />
        <Label htmlFor="cryoprotectant-mixture">Show cryoprotectant mixtures only</Label>
      </div>
      
      {/* Loading state */}
      {isLoading && (
        <div className="text-center py-12">
          <p>Loading mixtures...</p>
        </div>
      )}
      
      {/* Error state */}
      {isError && (
        <div className="text-center py-12 text-destructive">
          <p>Error loading mixtures: {error instanceof Error ? error.message : 'Unknown error'}</p>
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
      {data?.mixtures?.length === 0 && (
        <div className="text-center py-12">
          <p>No mixtures found.</p>
          {(params.search || params.is_cryoprotectant_mixture) && (
            <Button 
              variant="outline" 
              className="mt-4"
              onClick={() => {
                setParams({
                  page: 1,
                  limit: 12,
                  sort_by: 'name',
                  sort_order: 'asc',
                });
                setSearch('');
              }}
            >
              Clear Filters
            </Button>
          )}
        </div>
      )}
      
      {/* Mixtures grid */}
      {data?.mixtures && data.mixtures.length > 0 && (
        <div className="grid grid-cols-1 sm:grid-cols-2 md:grid-cols-3 gap-6">
          {data.mixtures.map((mixture) => (
            <MixtureCard key={mixture.id} mixture={mixture} />
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
  );
}
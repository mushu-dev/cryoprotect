/**
 * Enhanced Molecules List Component
 * Display molecules with real-time updates, advanced filtering, and performance monitoring
 */
import React, { useState, useEffect } from 'react';
import { useSimpleEnhancedMolecules, useSimpleClientMetrics, useSimpleRefreshMolecules } from '../../hooks/use-simple-enhanced-molecules';
import { EnhancedMoleculeCard, EnhancedMoleculeCardSkeleton } from './enhanced-molecule-card';
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
import { Card } from '../ui/card';
import { Badge } from '../ui/badge';

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

const RefreshIcon = () => (
  <svg xmlns="http://www.w3.org/2000/svg" width="24" height="24" viewBox="0 0 24 24" fill="none" 
    stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" 
    className="h-4 w-4">
    <path d="M3 12a9 9 0 0 1 9-9 9.75 9.75 0 0 1 6.74 2.74L21 8"/>
    <path d="M21 3v5h-5"/>
    <path d="M21 12a9 9 0 0 1-9 9 9.75 9.75 0 0 1-6.74-2.74L3 16"/>
    <path d="M3 21v-5h5"/>
  </svg>
);

const MetricsIcon = () => (
  <svg xmlns="http://www.w3.org/2000/svg" width="24" height="24" viewBox="0 0 24 24" fill="none" 
    stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" 
    className="h-4 w-4">
    <path d="M3 3v18h18"/>
    <path d="M7 12l4-4 4 4 4-4"/>
  </svg>
);

export function EnhancedMoleculesList() {
  const [params, setParams] = useState({
    page: 1,
    limit: 12,
    sort_by: 'name',
    sort_order: 'asc',
  });
  
  const [search, setSearch] = useState('');
  const [showMetrics, setShowMetrics] = useState(false);
  const [showDetails, setShowDetails] = useState(false);
  
  const { data, isLoading, error, refetch } = useSimpleEnhancedMolecules(params);
  const { data: clientMetrics } = useSimpleClientMetrics();
  const refreshMolecules = useSimpleRefreshMolecules();
  
  const isError = !!error;
  const isFetching = isLoading;
  
  // Auto-refresh every 2 minutes for real-time feel
  useEffect(() => {
    const interval = setInterval(() => {
      if (!isFetching) {
        refreshMolecules();
        refetch();
      }
    }, 2 * 60 * 1000); // 2 minutes
    
    return () => clearInterval(interval);
  }, [refreshMolecules, refetch, isFetching]);
  
  // Handle search input
  const handleSearch = () => {
    setParams(prev => ({ ...prev, search, page: 1 }));
  };
  
  // Handle filter for cryoprotectants only
  const handleCryoprotectantsFilter = (checked) => {
    setParams(prev => ({ 
      ...prev, 
      is_cryoprotectant: checked || undefined,
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
  
  // Handle molecule click
  const handleMoleculeClick = (molecule) => {
    console.log('Selected molecule:', molecule);
    // Could navigate to detail page or open modal
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
      {/* Performance Metrics (optional) */}
      {showMetrics && clientMetrics && (
        <Card className="p-4 bg-blue-50 border-blue-200">
          <div className="flex items-center justify-between mb-2">
            <h4 className="font-semibold flex items-center">
              <MetricsIcon className="mr-2" />
              Performance Metrics
            </h4>
            <Button
              variant="ghost"
              size="sm"
              onClick={() => setShowMetrics(false)}
            >
              ×
            </Button>
          </div>
          <div className="grid grid-cols-2 md:grid-cols-4 gap-4 text-sm">
            <div>
              <span className="text-gray-600">Cache Hit Rate:</span>
              <span className="ml-1 font-semibold text-green-600">
                {clientMetrics.hitRate}%
              </span>
            </div>
            <div>
              <span className="text-gray-600">Avg Response:</span>
              <span className="ml-1 font-semibold">
                {Math.round(clientMetrics.avgResponseTime)}ms
              </span>
            </div>
            <div>
              <span className="text-gray-600">Requests:</span>
              <span className="ml-1 font-semibold">
                {clientMetrics.requests}
              </span>
            </div>
            <div>
              <span className="text-gray-600">Cache Size:</span>
              <span className="ml-1 font-semibold">
                {clientMetrics.cacheSize}
              </span>
            </div>
          </div>
        </Card>
      )}

      {/* Controls */}
      <div className="flex flex-col space-y-4">
        {/* Search and filter row */}
        <div className="flex flex-col sm:flex-row gap-4">
          <div className="relative flex-grow">
            <div className="absolute left-2.5 top-2.5 h-4 w-4 text-muted-foreground">
              <SearchIcon />
            </div>
            <Input 
              placeholder="Search molecules..."
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
              disabled={isFetching}
            >
              {isFetching ? '...' : 'Search'}
            </Button>
            
            <Select defaultValue="name:asc" onValueChange={handleSortChange}>
              <SelectTrigger className="w-[165px]">
                <FilterIcon />
                <SelectValue placeholder="Sort by" />
              </SelectTrigger>
              <SelectContent>
                <SelectItem value="name:asc">Name (A-Z)</SelectItem>
                <SelectItem value="name:desc">Name (Z-A)</SelectItem>
                <SelectItem value="molecular_weight:asc">Weight (Low-High)</SelectItem>
                <SelectItem value="molecular_weight:desc">Weight (High-Low)</SelectItem>
                <SelectItem value="cryoprotectantScore:desc">Score (High-Low)</SelectItem>
                <SelectItem value="cryoprotectantScore:asc">Score (Low-High)</SelectItem>
                <SelectItem value="created_at:desc">Newest First</SelectItem>
                <SelectItem value="created_at:asc">Oldest First</SelectItem>
              </SelectContent>
            </Select>
            
            <Button
              variant="outline"
              onClick={refreshMolecules}
              disabled={isFetching}
              className="px-3"
            >
              <RefreshIcon className={isFetching ? 'animate-spin' : ''} />
            </Button>
          </div>
        </div>
        
        {/* Filter switches */}
        <div className="flex flex-wrap items-center gap-6">
          <div className="flex items-center space-x-2">
            <Switch 
              id="cryoprotectant" 
              checked={!!params.is_cryoprotectant} 
              onCheckedChange={handleCryoprotectantsFilter}
            />
            <Label htmlFor="cryoprotectant">Show cryoprotectants only</Label>
          </div>
          
          <div className="flex items-center space-x-2">
            <Switch 
              id="details" 
              checked={showDetails} 
              onCheckedChange={setShowDetails}
            />
            <Label htmlFor="details">Show detailed properties</Label>
          </div>
          
          <div className="flex items-center space-x-2">
            <Switch 
              id="metrics" 
              checked={showMetrics} 
              onCheckedChange={setShowMetrics}
            />
            <Label htmlFor="metrics">Show performance metrics</Label>
          </div>
        </div>
      </div>
      
      {/* Status indicators */}
      {data && (
        <div className="flex items-center gap-4 text-sm text-gray-600">
          <span>
            Showing {data.molecules?.length || 0} of {data.total || 0} molecules
          </span>
          {clientMetrics && clientMetrics.hitRate >= 80 && (
            <Badge variant="outline" className="bg-green-100 text-green-800 border-green-200">
              🚀 Fast Cache
            </Badge>
          )}
          {isFetching && (
            <Badge variant="outline">
              🔄 Updating...
            </Badge>
          )}
        </div>
      )}
      
      {/* Loading state */}
      {isLoading && (
        <div className="grid grid-cols-1 sm:grid-cols-2 md:grid-cols-3 lg:grid-cols-4 gap-6">
          {[...Array(params.limit)].map((_, i) => (
            <EnhancedMoleculeCardSkeleton key={i} />
          ))}
        </div>
      )}
      
      {/* Error state */}
      {isError && (
        <div className="text-center py-12">
          <Card className="p-6 border-red-200 bg-red-50">
            <p className="text-red-600 mb-4">
              Error loading enhanced molecules: {error instanceof Error ? error.message : 'Unknown error'}
            </p>
            <div className="flex gap-2 justify-center">
              <Button 
                variant="outline" 
                onClick={refetch} // Retry
              >
                Retry
              </Button>
              <Button 
                variant="outline"
                onClick={refreshMolecules}
              >
                <RefreshIcon className="mr-2" />
                Refresh Cache
              </Button>
            </div>
          </Card>
        </div>
      )}
      
      {/* Empty state */}
      {data?.molecules?.length === 0 && !isLoading && (
        <div className="text-center py-12">
          <Card className="p-6">
            <p className="text-gray-600 mb-4">No molecules found.</p>
            {(params.search || params.is_cryoprotectant) && (
              <Button 
                variant="outline"
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
          </Card>
        </div>
      )}
      
      {/* Molecules grid */}
      {data?.molecules && data.molecules.length > 0 && (
        <div className="grid grid-cols-1 sm:grid-cols-2 md:grid-cols-3 lg:grid-cols-4 gap-6">
          {data.molecules.map((molecule) => (
            <EnhancedMoleculeCard 
              key={molecule.id} 
              molecule={molecule}
              onClick={handleMoleculeClick}
              showDetails={showDetails}
            />
          ))}
        </div>
      )}
      
      {/* Pagination */}
      {data && data.total > data.limit && (
        <Pagination className="mt-8">
          <PaginationContent>
            <PaginationPrevious 
              onClick={() => handlePageChange(Math.max(1, (params.page || 1) - 1))}
              disabled={params.page <= 1}
            />
            
            {getPaginationItems()}
            
            <PaginationNext 
              onClick={() => handlePageChange(Math.min(totalPages, (params.page || 1) + 1))}
              disabled={params.page >= totalPages}
            />
          </PaginationContent>
        </Pagination>
      )}
    </div>
  );
}
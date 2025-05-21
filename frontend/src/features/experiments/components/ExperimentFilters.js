import React, { useState } from 'react';

/**
 * Component for filtering experiments by various criteria
 */
export default function ExperimentFilters({ filters, setFilters, experiments }) {
  const [filtersVisible, setFiltersVisible] = useState(true);

  // Extract unique values for filter dropdowns
  const getUniqueValues = (fieldName) => {
    if (!experiments || !Array.isArray(experiments)) return [];
    
    const values = new Set();
    
    experiments.forEach(exp => {
      if (fieldName === 'cryoprotectants') {
        // Special handling for cryoprotectants array
        if (exp.cryoprotectants && Array.isArray(exp.cryoprotectants)) {
          exp.cryoprotectants.forEach(cp => {
            values.add(cp.name);
          });
        }
      } else if (exp[fieldName]) {
        values.add(exp[fieldName]);
      }
    });
    
    return Array.from(values);
  };

  const uniqueStatuses = getUniqueValues('status');
  const uniqueCellTypes = getUniqueValues('cellType');
  const uniqueCryoprotectants = getUniqueValues('cryoprotectants');
  
  // Get unique protocol names
  const uniqueProtocols = experiments
    ? [...new Set(experiments.map(exp => exp.protocol?.name).filter(Boolean))]
    : [];

  // Handle filter changes
  const handleFilterChange = (field, value) => {
    setFilters(prev => ({
      ...prev,
      [field]: value
    }));
  };

  // Handle threshold changes (for numeric filters)
  const handleThresholdChange = (field, value) => {
    const numValue = value === '' ? null : Number(value);
    setFilters(prev => ({
      ...prev,
      [field]: numValue
    }));
  };

  // Clear all filters
  const clearFilters = () => {
    setFilters({});
  };

  return (
    <div className="mb-6 border rounded-lg overflow-hidden" data-testid="filter-controls">
      <div 
        className="p-4 bg-gray-50 border-b font-medium flex items-center justify-between cursor-pointer"
        onClick={() => setFiltersVisible(!filtersVisible)}
      >
        <div className="flex items-center">
          <svg width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" className="mr-2">
            <polygon points="22 3 2 3 10 12.46 10 19 14 21 14 12.46 22 3"/>
          </svg>
          <span>Filter Experiments</span>
        </div>
        <button className="text-gray-500 hover:text-gray-700">
          {filtersVisible ? (
            <svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
              <polyline points="18 15 12 9 6 15"/>
            </svg>
          ) : (
            <svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
              <polyline points="6 9 12 15 18 9"/>
            </svg>
          )}
        </button>
      </div>
      
      {filtersVisible && (
        <div className="p-4 bg-white">
          <div className="grid grid-cols-1 md:grid-cols-3 gap-4 mb-4">
            {/* Status Filter */}
            <div>
              <label htmlFor="statusFilter" className="block text-sm font-medium mb-1">
                Status
              </label>
              <select
                id="statusFilter"
                className="w-full px-3 py-2 border rounded-md"
                value={filters.status || ''}
                onChange={(e) => handleFilterChange('status', e.target.value || null)}
              >
                <option value="">Any Status</option>
                {uniqueStatuses.map(status => (
                  <option key={status} value={status}>{status}</option>
                ))}
              </select>
            </div>
            
            {/* Cell Type Filter */}
            <div>
              <label htmlFor="cellTypeFilter" className="block text-sm font-medium mb-1">
                Cell Type
              </label>
              <select
                id="cellTypeFilter"
                className="w-full px-3 py-2 border rounded-md"
                value={filters.cellType || ''}
                onChange={(e) => handleFilterChange('cellType', e.target.value || null)}
              >
                <option value="">Any Cell Type</option>
                {uniqueCellTypes.map(cellType => (
                  <option key={cellType} value={cellType}>{cellType}</option>
                ))}
              </select>
            </div>
            
            {/* Cryoprotectant Filter */}
            <div>
              <label htmlFor="cryoprotectantFilter" className="block text-sm font-medium mb-1">
                Cryoprotectant
              </label>
              <select
                id="cryoprotectantFilter"
                className="w-full px-3 py-2 border rounded-md"
                value={filters.cryoprotectant || ''}
                onChange={(e) => handleFilterChange('cryoprotectant', e.target.value || null)}
              >
                <option value="">Any Cryoprotectant</option>
                {uniqueCryoprotectants.map(cp => (
                  <option key={cp} value={cp}>{cp}</option>
                ))}
              </select>
            </div>
          </div>
          
          <div className="grid grid-cols-1 md:grid-cols-3 gap-4 mb-4">
            {/* Date Range Filter */}
            <div>
              <label htmlFor="dateRangeFilter" className="block text-sm font-medium mb-1">
                Date Range
              </label>
              <input
                id="dateRangeFilter"
                type="date"
                className="w-full px-3 py-2 border rounded-md"
                value={filters.dateFrom || ''}
                onChange={(e) => handleFilterChange('dateFrom', e.target.value || null)}
              />
            </div>
            
            {/* Protocol Filter */}
            <div>
              <label htmlFor="protocolFilter" className="block text-sm font-medium mb-1">
                Protocol
              </label>
              <select
                id="protocolFilter"
                className="w-full px-3 py-2 border rounded-md"
                value={filters.protocol || ''}
                onChange={(e) => handleFilterChange('protocol', e.target.value || null)}
              >
                <option value="">Any Protocol</option>
                {uniqueProtocols.map(protocol => (
                  <option key={protocol} value={protocol}>{protocol}</option>
                ))}
              </select>
            </div>
            
            {/* Sort Order */}
            <div className="md:col-span-1" data-testid="sort-controls">
              <label htmlFor="sortOrder" className="block text-sm font-medium mb-1">
                Sort By
              </label>
              <select
                id="sortOrder"
                className="w-full px-3 py-2 border rounded-md"
                value={filters.sortBy || 'date'}
                onChange={(e) => handleFilterChange('sortBy', e.target.value)}
              >
                <option value="date">Date (newest first)</option>
                <option value="dateAsc">Date (oldest first)</option>
                <option value="title">Title (A-Z)</option>
                <option value="viability">Viability (highest first)</option>
                <option value="recovery">Recovery (highest first)</option>
              </select>
            </div>
          </div>
          
          <div className="grid grid-cols-1 md:grid-cols-2 gap-4 mb-4">
            {/* Viability Threshold */}
            <div>
              <label htmlFor="viabilityFilter" className="block text-sm font-medium mb-1">
                Minimum Viability (%)
              </label>
              <input
                id="viabilityFilter"
                type="number"
                min="0"
                max="100"
                className="w-full px-3 py-2 border rounded-md"
                value={filters.minViability || ''}
                onChange={(e) => handleThresholdChange('minViability', e.target.value)}
                placeholder="Any viability"
              />
            </div>
            
            {/* Recovery Threshold */}
            <div>
              <label htmlFor="recoveryFilter" className="block text-sm font-medium mb-1">
                Minimum Recovery (%)
              </label>
              <input
                id="recoveryFilter"
                type="number"
                min="0"
                max="100"
                className="w-full px-3 py-2 border rounded-md"
                value={filters.minRecovery || ''}
                onChange={(e) => handleThresholdChange('minRecovery', e.target.value)}
                placeholder="Any recovery rate"
              />
            </div>
          </div>
          
          <div className="flex justify-end">
            <button
              className="px-4 py-2 text-sm font-medium text-blue-600 hover:text-blue-800"
              onClick={clearFilters}
            >
              Clear Filters
            </button>
          </div>
        </div>
      )}
    </div>
  );
}
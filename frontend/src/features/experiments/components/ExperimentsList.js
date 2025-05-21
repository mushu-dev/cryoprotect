import React, { useState, useEffect } from 'react';
import Link from 'next/link';
import useExperimentData from '../hooks/useExperimentData';
import useConvexExperimentData from '../hooks/useConvexExperimentData';
import ExperimentCard from './ExperimentCard';
import ExperimentFilters from './ExperimentFilters';
import ExperimentComparison from './ExperimentComparison';
import { analyzeCryoprotectantEffectiveness } from '../utils/data-transformation';

/**
 * Main experiments list page with filtering, comparison, and data analysis
 * Supports both traditional API and Convex for real-time collaboration
 */
export default function ExperimentsList() {
  // Environment configuration for determining data source
  const useConvex = process.env.NEXT_PUBLIC_USE_CONVEX === 'true';
  
  // Standard API data hook
  const standardData = useExperimentData();
  
  // Convex real-time data hook
  const convexData = useConvexExperimentData(
    {}, // initial filters
    {   // options
      includeResults: true,
      includeProtocol: true,
      includeMixture: false
    }
  );

  // Use the appropriate data source based on configuration
  const { 
    experiments, 
    filteredExperiments, 
    loading, 
    error, 
    filters, 
    setFilters,
    selectedExperiments,
    toggleExperimentSelection,
    selectedExperimentData,
    viewMode,
    setViewMode
  } = useConvex ? convexData : standardData;

  // UI state
  const [showComparison, setShowComparison] = useState(false);
  const [showEffectiveness, setShowEffectiveness] = useState(false);

  // Handle removing experiments from comparison
  const handleRemoveFromComparison = (experimentId) => {
    if (Array.isArray(experimentId)) {
      // Clear all
      if (useConvex) {
        convexData.setSelectedExperiments([]);
      } else {
        standardData.setSelectedExperiments([]);
      }
    } else {
      // Remove specific experiment
      toggleExperimentSelection(experimentId);
    }
  };

  // Generate cryoprotectant effectiveness data
  const effectivenessData = analyzeCryoprotectantEffectiveness(experiments);

  if (loading) {
    return (
      <div className="py-8 text-center">
        <div className="animate-spin rounded-full h-12 w-12 border-b-2 border-primary mx-auto"></div>
        <p className="mt-4 text-gray-500">Loading experiments...</p>
      </div>
    );
  }

  if (error) {
    return (
      <div className="py-8 text-center">
        <div className="text-red-500 mb-2">Error loading experiments</div>
        <p className="text-gray-500">{error}</p>
      </div>
    );
  }

  return (
    <div className="pb-12">
      <div className="flex flex-col md:flex-row justify-between items-start md:items-center mb-8">
        <div>
          <h1 className="text-3xl font-bold mb-2">Experiments</h1>
          <p className="text-muted-foreground">
            Design, track, and analyze cryopreservation experiments
          </p>
        </div>
        
        <div className="mt-4 md:mt-0 flex flex-wrap gap-2">
          <div className="flex rounded-md overflow-hidden">
            <button
              onClick={() => setViewMode('grid')}
              className={`px-3 py-2 text-sm font-medium ${
                viewMode === 'grid' 
                  ? 'bg-primary text-white' 
                  : 'bg-white text-gray-700 hover:bg-gray-100'
              } border border-r-0`}
            >
              <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" className="mr-1 inline-block">
                <rect x="3" y="3" width="7" height="7"></rect>
                <rect x="14" y="3" width="7" height="7"></rect>
                <rect x="14" y="14" width="7" height="7"></rect>
                <rect x="3" y="14" width="7" height="7"></rect>
              </svg>
              View as Grid
            </button>
            <button
              onClick={() => setViewMode('list')}
              className={`px-3 py-2 text-sm font-medium ${
                viewMode === 'list' 
                  ? 'bg-primary text-white' 
                  : 'bg-white text-gray-700 hover:bg-gray-100'
              } border`}
            >
              <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" className="mr-1 inline-block">
                <line x1="8" y1="6" x2="21" y2="6"></line>
                <line x1="8" y1="12" x2="21" y2="12"></line>
                <line x1="8" y1="18" x2="21" y2="18"></line>
                <line x1="3" y1="6" x2="3.01" y2="6"></line>
                <line x1="3" y1="12" x2="3.01" y2="12"></line>
                <line x1="3" y1="18" x2="3.01" y2="18"></line>
              </svg>
              View as List
            </button>
          </div>
          
          <button
            onClick={() => setShowComparison(!showComparison)}
            className="px-3 py-2 text-sm font-medium rounded-md bg-blue-50 text-blue-700 border border-blue-200 hover:bg-blue-100"
          >
            {showComparison ? 'Hide Comparison' : 'Show Comparison'}
            {selectedExperiments.length > 0 && (
              <span className="ml-1 px-1.5 py-0.5 text-xs rounded-full bg-blue-700 text-white">
                {selectedExperiments.length}
              </span>
            )}
          </button>
          
          <button
            onClick={() => setShowEffectiveness(!showEffectiveness)}
            className="px-3 py-2 text-sm font-medium rounded-md bg-purple-50 text-purple-700 border border-purple-200 hover:bg-purple-100"
          >
            {showEffectiveness ? 'Hide Analysis' : 'Cryoprotectant Analysis'}
          </button>
          
          <Link href="/experiments/create">
            <a className="inline-flex items-center justify-center rounded-md bg-primary px-4 py-2 text-sm font-medium text-primary-foreground shadow transition-colors hover:bg-primary/90">
              Create New Experiment
            </a>
          </Link>
        </div>
      </div>
      
      {/* Filters */}
      <ExperimentFilters 
        filters={filters} 
        setFilters={setFilters} 
        experiments={experiments}
      />
      
      {/* Comparison panel */}
      {showComparison && (
        <div className="mb-8">
          <ExperimentComparison 
            experiments={selectedExperimentData}
            onRemove={handleRemoveFromComparison}
          />
        </div>
      )}
      
      {/* Effectiveness analysis panel */}
      {showEffectiveness && (
        <div className="mb-8 bg-white rounded-lg border shadow-sm overflow-hidden">
          <div className="p-4 bg-gray-50 border-b font-medium flex justify-between items-center">
            <h3>Cryoprotectant Effectiveness Analysis</h3>
            <button
              onClick={() => setShowEffectiveness(false)}
              className="text-xs text-gray-500 hover:text-gray-700"
            >
              Close
            </button>
          </div>
          
          <div className="p-4">
            <div className="mb-4">
              <p className="text-sm text-gray-600">
                This analysis shows the average performance metrics across all experiments using each cryoprotectant.
              </p>
            </div>
            
            <div className="overflow-x-auto">
              <table className="w-full text-sm">
                <thead>
                  <tr className="bg-gray-50">
                    <th className="px-4 py-2 text-left font-medium text-gray-500">Cryoprotectant</th>
                    <th className="px-4 py-2 text-left font-medium text-gray-500">Avg. Viability</th>
                    <th className="px-4 py-2 text-left font-medium text-gray-500">Avg. Recovery</th>
                    <th className="px-4 py-2 text-left font-medium text-gray-500">Avg. Functionality</th>
                    <th className="px-4 py-2 text-left font-medium text-gray-500">Experiments</th>
                  </tr>
                </thead>
                <tbody className="divide-y divide-gray-200">
                  {effectivenessData.map((cp, index) => (
                    <tr key={index} className="hover:bg-gray-50">
                      <td className="px-4 py-2 font-medium">{cp.name}</td>
                      <td className={`px-4 py-2 ${Number(cp.averageViability) > 90 ? 'text-green-600 font-medium' : ''}`}>
                        {cp.averageViability}%
                      </td>
                      <td className={`px-4 py-2 ${Number(cp.averageRecovery) > 90 ? 'text-green-600 font-medium' : ''}`}>
                        {cp.averageRecovery}%
                      </td>
                      <td className={`px-4 py-2 ${Number(cp.averageFunctionality) > 90 ? 'text-green-600 font-medium' : ''}`}>
                        {cp.averageFunctionality}%
                      </td>
                      <td className="px-4 py-2">{cp.experimentCount}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </div>
        </div>
      )}
      
      {/* Experiments list */}
      {filteredExperiments.length === 0 ? (
        <div className="py-12 text-center border rounded-lg bg-gray-50">
          <p className="text-gray-500 mb-2">No experiments match your filters</p>
          <button
            onClick={() => setFilters({})}
            className="text-blue-600 hover:text-blue-800 hover:underline text-sm"
          >
            Clear all filters
          </button>
        </div>
      ) : (
        <>
          <div className="text-sm text-gray-500 mb-4">
            Showing {filteredExperiments.length} experiments
            {Object.keys(filters).length > 0 && ' (filtered)'}
          </div>
          
          {viewMode === 'grid' ? (
            <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
              {filteredExperiments.map(experiment => (
                <div key={experiment.id} className="experiment-card">
                  <ExperimentCard
                    experiment={experiment}
                    isSelected={selectedExperiments.includes(experiment.id)}
                    onSelect={toggleExperimentSelection}
                    showSelectOption={true}
                  />
                </div>
              ))}
            </div>
          ) : (
            <div className="bg-white rounded-lg border shadow-sm overflow-hidden">
              <table className="w-full text-sm">
                <thead>
                  <tr className="border-b">
                    <th className="px-4 py-3 text-left font-medium text-gray-500">Title</th>
                    <th className="px-4 py-3 text-left font-medium text-gray-500">Status</th>
                    <th className="px-4 py-3 text-left font-medium text-gray-500">Cell Type</th>
                    <th className="px-4 py-3 text-left font-medium text-gray-500">Cryoprotectants</th>
                    <th className="px-4 py-3 text-left font-medium text-gray-500">Viability</th>
                    <th className="px-4 py-3 text-left font-medium text-gray-500">Date</th>
                    <th className="px-4 py-3 text-center font-medium text-gray-500">Compare</th>
                    <th className="px-4 py-3 text-right font-medium text-gray-500">Actions</th>
                  </tr>
                </thead>
                <tbody>
                  {filteredExperiments.map(experiment => (
                    <tr key={experiment.id} className="border-b last:border-b-0 hover:bg-gray-50 experiment-list-item">
                      <td className="px-4 py-3 font-medium">
                        <Link href={`/experiments/${experiment.id}`}>
                          <a className="hover:text-blue-600">{experiment.title}</a>
                        </Link>
                      </td>
                      <td className="px-4 py-3">
                        <span className={`inline-flex items-center rounded-full px-2.5 py-0.5 text-xs font-medium ${
                          experiment.status === 'Completed' 
                            ? 'bg-green-100 text-green-800' 
                            : 'bg-blue-100 text-blue-800'
                        }`}>
                          {experiment.status}
                        </span>
                      </td>
                      <td className="px-4 py-3 text-gray-600">{experiment.cellType}</td>
                      <td className="px-4 py-3 text-gray-600">
                        {experiment.cryoprotectants.map(cp => cp.name).join(', ')}
                      </td>
                      <td className="px-4 py-3 font-medium">{experiment.results.viability}</td>
                      <td className="px-4 py-3 text-gray-600">{experiment.date}</td>
                      <td className="px-4 py-3 text-center">
                        <button
                          onClick={() => toggleExperimentSelection(experiment.id)}
                          className={`w-5 h-5 rounded-md border ${
                            selectedExperiments.includes(experiment.id)
                              ? 'bg-blue-500 border-blue-500 text-white'
                              : 'border-gray-300 hover:border-blue-400'
                          }`}
                        >
                          {selectedExperiments.includes(experiment.id) && (
                            <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="currentColor" width="16" height="16" className="mx-auto">
                              <path fillRule="evenodd" d="M19.916 4.626a.75.75 0 01.208 1.04l-9 13.5a.75.75 0 01-1.154.114l-6-6a.75.75 0 011.06-1.06l5.353 5.353 8.493-12.739a.75.75 0 011.04-.208z" clipRule="evenodd" />
                            </svg>
                          )}
                        </button>
                      </td>
                      <td className="px-4 py-3 text-right">
                        <Link href={`/experiments/${experiment.id}`}>
                          <a className="text-blue-600 hover:text-blue-800 hover:underline text-sm mr-3">
                            View
                          </a>
                        </Link>
                        <Link href={`/experiments/${experiment.id}/edit`}>
                          <a className="text-gray-600 hover:text-gray-800 hover:underline text-sm">
                            Edit
                          </a>
                        </Link>
                      </td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          )}
        </>
      )}
    </div>
  );
}
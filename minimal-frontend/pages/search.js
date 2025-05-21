import React, { useState, useEffect } from 'react';
import { useRouter } from 'next/router';
import Link from 'next/link';
import Layout from '../components/Layout';
import SearchBar from '../components/SearchBar';
import MoleculeCard from '../components/MoleculeCard';
import MixtureCard from '../components/MixtureCard';
import Loading from '../components/Loading';
import ErrorMessage from '../components/ErrorMessage';
import { getMolecules, getMixtures } from '../utils/api';

const searchFilters = [
  { label: 'All', value: 'all' },
  { label: 'Molecules', value: 'molecules' },
  { label: 'Mixtures', value: 'mixtures' },
  { label: 'Cryoprotectants Only', value: 'cryoprotectants' }
];

export default function SearchPage() {
  const router = useRouter();
  const { q, filter } = router.query;
  
  const [searchResults, setSearchResults] = useState({
    molecules: [],
    mixtures: [],
    loading: false,
    error: null
  });
  
  const [searchTerm, setSearchTerm] = useState(q || '');
  const [activeFilter, setActiveFilter] = useState(filter || 'all');
  
  // Initial data fetch on mount and when query parameters change
  useEffect(() => {
    if (q) {
      setSearchTerm(q);
      handleSearch(q, filter || 'all');
    }
  }, [q, filter]);
  
  const handleSearch = async (term, filterValue) => {
    if (!term) return;
    
    try {
      setSearchResults(prev => ({ ...prev, loading: true, error: null }));
      setActiveFilter(filterValue);
      
      // Update the URL with search parameters
      router.push({
        pathname: '/search',
        query: { q: term, filter: filterValue }
      }, undefined, { shallow: true });
      
      // Fetch molecules if the filter allows it
      let molecules = [];
      if (['all', 'molecules', 'cryoprotectants'].includes(filterValue)) {
        molecules = await getMolecules();
        
        // Filter by search term
        molecules = molecules.filter(molecule => {
          const matches = [
            molecule.name,
            molecule.formula,
            molecule.pubchem_cid,
            molecule.description
          ].filter(Boolean).some(field => 
            field.toLowerCase().includes(term.toLowerCase())
          );
          
          // Additional filtering for cryoprotectants
          if (filterValue === 'cryoprotectants') {
            return matches && molecule.is_cryoprotectant;
          }
          
          return matches;
        });
      }
      
      // Fetch mixtures if the filter allows it
      let mixtures = [];
      if (['all', 'mixtures'].includes(filterValue)) {
        mixtures = await getMixtures();
        
        // Filter by search term
        mixtures = mixtures.filter(mixture => {
          const hasMatchingComponents = mixture.components && mixture.components.some(comp => 
            comp.molecule && comp.molecule.name && 
            comp.molecule.name.toLowerCase().includes(term.toLowerCase())
          );
          
          return [
            mixture.name,
            mixture.description
          ].filter(Boolean).some(field => 
            field.toLowerCase().includes(term.toLowerCase())
          ) || hasMatchingComponents;
        });
      }
      
      setSearchResults({
        molecules,
        mixtures,
        loading: false,
        error: null
      });
    } catch (error) {
      console.error('Search failed:', error);
      setSearchResults(prev => ({
        ...prev,
        loading: false,
        error: 'Failed to fetch search results. Please try again.'
      }));
    }
  };
  
  const handleSearchSubmit = (term, filterValue) => {
    setSearchTerm(term);
    handleSearch(term, filterValue);
  };
  
  // Calculate the total number of results
  const totalResults = searchResults.molecules.length + searchResults.mixtures.length;

  return (
    <Layout title="Search">
      <div className="page-header">
        <h1 className="title">Search</h1>
        <p className="description">
          Search for molecules, mixtures, and cryoprotectants in our database.
        </p>
      </div>
      
      <div className="search-container">
        <SearchBar 
          onSearch={handleSearchSubmit} 
          placeholder="Search by name, formula, or ID..."
          initialValue={searchTerm}
          filters={searchFilters}
        />
      </div>
      
      {/* Search Results */}
      <div className="results-container">
        {searchResults.loading ? (
          <Loading />
        ) : searchResults.error ? (
          <ErrorMessage message={searchResults.error} onRetry={() => handleSearch(searchTerm, activeFilter)} />
        ) : searchTerm ? (
          <>
            <div className="results-header">
              <h2 className="results-title">Search Results</h2>
              <p className="results-count">
                {totalResults} {totalResults === 1 ? 'result' : 'results'} for "{searchTerm}"
              </p>
            </div>
            
            {/* Molecule Results */}
            {(activeFilter === 'all' || activeFilter === 'molecules' || activeFilter === 'cryoprotectants') && 
              searchResults.molecules.length > 0 && (
                <div className="results-section">
                  <h3 className="section-title">Molecules</h3>
                  <div className="molecules-grid">
                    {searchResults.molecules.map(molecule => (
                      <MoleculeCard key={molecule.id} molecule={molecule} />
                    ))}
                  </div>
                </div>
              )
            }
            
            {/* Mixture Results */}
            {(activeFilter === 'all' || activeFilter === 'mixtures') && 
              searchResults.mixtures.length > 0 && (
                <div className="results-section">
                  <h3 className="section-title">Mixtures</h3>
                  <div className="mixtures-grid">
                    {searchResults.mixtures.map(mixture => (
                      <MixtureCard key={mixture.id} mixture={mixture} />
                    ))}
                  </div>
                </div>
              )
            }
            
            {/* No Results Message */}
            {totalResults === 0 && (
              <div className="no-results">
                <p>No results found for "{searchTerm}".</p>
                <p>Try using different keywords or filters.</p>
                <div className="search-suggestions">
                  <p>Popular searches:</p>
                  <div className="suggestion-buttons">
                    <button 
                      className="suggestion-button" 
                      onClick={() => handleSearchSubmit('glycerol', activeFilter)}
                    >
                      Glycerol
                    </button>
                    <button 
                      className="suggestion-button" 
                      onClick={() => handleSearchSubmit('DMSO', activeFilter)}
                    >
                      DMSO
                    </button>
                    <button 
                      className="suggestion-button" 
                      onClick={() => handleSearchSubmit('cryoprotectant', activeFilter)}
                    >
                      Cryoprotectant
                    </button>
                  </div>
                </div>
              </div>
            )}
          </>
        ) : (
          <div className="search-prompt">
            <p>Enter a search term to find molecules and mixtures.</p>
            <p>You can search by name, formula, or description.</p>
          </div>
        )}
      </div>
      
      <style jsx>{`
        .search-container {
          width: 100%;
          max-width: 800px;
          margin: 0 auto 20px;
        }
        
        .results-container {
          width: 100%;
          max-width: 1200px;
        }
        
        .results-header {
          margin-bottom: 20px;
          padding-bottom: 10px;
          border-bottom: 1px solid #e2e8f0;
        }
        
        .results-title {
          font-size: 1.5rem;
          color: #2d3748;
          margin-bottom: 5px;
        }
        
        .results-count {
          color: #718096;
          font-size: 0.9rem;
        }
        
        .results-section {
          margin-bottom: 30px;
        }
        
        .section-title {
          font-size: 1.2rem;
          color: #4a5568;
          margin-bottom: 15px;
          padding-left: 10px;
          border-left: 3px solid #3182ce;
        }
        
        .no-results {
          text-align: center;
          padding: 40px 20px;
          color: #4a5568;
        }
        
        .search-suggestions {
          margin-top: 20px;
        }
        
        .suggestion-buttons {
          display: flex;
          flex-wrap: wrap;
          justify-content: center;
          gap: 10px;
          margin-top: 10px;
        }
        
        .suggestion-button {
          padding: 8px 16px;
          background-color: #f7fafc;
          border: 1px solid #e2e8f0;
          border-radius: 4px;
          color: #4a5568;
          font-size: 0.9rem;
          cursor: pointer;
          transition: all 0.2s;
        }
        
        .suggestion-button:hover {
          background-color: #edf2f7;
          border-color: #cbd5e0;
        }
        
        .search-prompt {
          text-align: center;
          padding: 40px 20px;
          color: #4a5568;
        }
      `}</style>
    </Layout>
  );
}
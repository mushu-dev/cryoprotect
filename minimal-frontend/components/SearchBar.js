import React, { useState } from 'react';

/**
 * A reusable search bar component
 * @param {Object} props - Component props
 * @param {Function} props.onSearch - Function to call when search is submitted
 * @param {string} props.placeholder - Placeholder text for search input
 * @param {string} props.initialValue - Initial value for search input
 * @param {Array} props.filters - Array of filter options
 */
export default function SearchBar({ onSearch, placeholder = 'Search...', initialValue = '', filters = [] }) {
  const [searchTerm, setSearchTerm] = useState(initialValue);
  const [selectedFilter, setSelectedFilter] = useState(filters.length > 0 ? filters[0].value : '');
  const [isExpanded, setIsExpanded] = useState(false);

  const handleSubmit = (e) => {
    e.preventDefault();
    onSearch(searchTerm, selectedFilter);
  };

  const handleChange = (e) => {
    setSearchTerm(e.target.value);
  };

  const handleFilterChange = (e) => {
    setSelectedFilter(e.target.value);
  };

  const handleKeyDown = (e) => {
    if (e.key === 'Enter') {
      handleSubmit(e);
    }
  };

  return (
    <div className={`search-bar ${isExpanded ? 'expanded' : ''}`}>
      <div className="search-container">
        <form onSubmit={handleSubmit} className="search-form">
          <div className="search-input-container">
            <input
              type="text"
              className="search-input"
              value={searchTerm}
              onChange={handleChange}
              onKeyDown={handleKeyDown}
              placeholder={placeholder}
              onFocus={() => setIsExpanded(true)}
              onBlur={() => setTimeout(() => setIsExpanded(false), 200)}
            />
            <button type="submit" className="search-button">
              <svg className="search-icon" fill="none" stroke="currentColor" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth="2" d="M21 21l-6-6m2-5a7 7 0 11-14 0 7 7 0 0114 0z"></path>
              </svg>
            </button>
          </div>
          
          {filters.length > 0 && (
            <div className="filter-container">
              <select 
                className="filter-select"
                value={selectedFilter}
                onChange={handleFilterChange}
                aria-label="Search filter"
              >
                {filters.map((filter) => (
                  <option key={filter.value} value={filter.value}>
                    {filter.label}
                  </option>
                ))}
              </select>
            </div>
          )}
        </form>
      </div>
      
      {isExpanded && (
        <div className="search-suggestions">
          <div className="suggestions-header">Try searching for:</div>
          <div className="suggestion-items">
            <div className="suggestion-item" onClick={() => setSearchTerm('Glycerol')}>Glycerol</div>
            <div className="suggestion-item" onClick={() => setSearchTerm('DMSO')}>DMSO</div>
            <div className="suggestion-item" onClick={() => setSearchTerm('C3H8O3')}>C3H8O3</div>
            <div className="suggestion-item" onClick={() => setSearchTerm('cryoprotectant')}>cryoprotectant</div>
          </div>
        </div>
      )}
      
      <style jsx>{`
        .search-bar {
          position: relative;
          width: 100%;
          max-width: 600px;
          margin: 0 auto 20px;
        }
        
        .search-container {
          width: 100%;
        }
        
        .search-form {
          display: flex;
          gap: 10px;
        }
        
        .search-input-container {
          flex: 1;
          position: relative;
        }
        
        .search-input {
          width: 100%;
          padding: 12px 45px 12px 15px;
          border: 1px solid #e2e8f0;
          border-radius: 8px;
          font-size: 16px;
          transition: all 0.3s;
          box-shadow: 0 2px 4px rgba(0, 0, 0, 0.05);
        }
        
        .search-input:focus {
          outline: none;
          border-color: #3182ce;
          box-shadow: 0 2px 8px rgba(49, 130, 206, 0.2);
        }
        
        .search-button {
          position: absolute;
          right: 10px;
          top: 50%;
          transform: translateY(-50%);
          background: none;
          border: none;
          color: #4a5568;
          cursor: pointer;
        }
        
        .search-icon {
          width: 20px;
          height: 20px;
        }
        
        .filter-container {
          width: 150px;
        }
        
        .filter-select {
          width: 100%;
          padding: 12px 10px;
          border: 1px solid #e2e8f0;
          border-radius: 8px;
          background-color: white;
          font-size: 16px;
          cursor: pointer;
          box-shadow: 0 2px 4px rgba(0, 0, 0, 0.05);
        }
        
        .filter-select:focus {
          outline: none;
          border-color: #3182ce;
        }
        
        .search-suggestions {
          position: absolute;
          top: calc(100% + 5px);
          left: 0;
          width: 100%;
          background-color: white;
          border-radius: 8px;
          border: 1px solid #e2e8f0;
          box-shadow: 0 4px 12px rgba(0, 0, 0, 0.1);
          padding: 10px;
          z-index: 1000;
        }
        
        .suggestions-header {
          font-size: 14px;
          color: #718096;
          margin-bottom: 8px;
          padding: 0 5px;
        }
        
        .suggestion-items {
          display: flex;
          flex-wrap: wrap;
          gap: 8px;
        }
        
        .suggestion-item {
          background-color: #f7fafc;
          border: 1px solid #e2e8f0;
          border-radius: 4px;
          padding: 5px 10px;
          font-size: 14px;
          cursor: pointer;
          transition: background-color 0.2s;
        }
        
        .suggestion-item:hover {
          background-color: #edf2f7;
        }
        
        /* Mobile styles */
        @media (max-width: 576px) {
          .search-form {
            flex-direction: column;
            gap: 8px;
          }
          
          .filter-container {
            width: 100%;
          }
          
          .search-suggestions {
            width: 100%;
            position: static;
            margin-top: 10px;
            box-shadow: none;
            border-radius: 6px;
          }
          
          .suggestion-items {
            justify-content: center;
          }
          
          .search-input {
            padding: 10px 36px 10px 12px;
            border-radius: 6px;
          }
          
          .search-button {
            right: 8px;
          }
          
          /* Improve touch target sizes */
          .suggestion-item {
            padding: 8px 12px;
            margin-bottom: 5px;
          }
        }
      `}</style>
    </div>
  );
}
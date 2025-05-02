/**
 * CryoProtect Analyzer Web Interface
 * Pagination Component
 * 
 * Provides pagination controls and manages API requests for paginated data.
 */

class PaginationComponent {
  /**
   * Initialize pagination component
   * @param {string} containerId - ID of container element
   * @param {string} apiEndpoint - API endpoint for fetching data
   * @param {Function} renderCallback - Function to render items
   * @param {Object} options - Additional options
   */
  constructor(containerId, apiEndpoint, renderCallback, options = {}) {
    this.container = document.getElementById(containerId);
    if (!this.container) {
      console.error(`Container element #${containerId} not found`);
      return;
    }
    
    this.apiEndpoint = apiEndpoint;
    this.renderCallback = renderCallback;
    
    // Default options
    this.options = {
      limit: 20,
      initialPage: 1,
      searchEnabled: true,
      sortEnabled: true,
      defaultSort: 'name',
      defaultOrder: 'asc',
      ...options
    };
    
    // State
    this.currentPage = this.options.initialPage;
    this.totalPages = 1;
    this.isLoading = false;
    this.searchTerm = '';
    this.sortField = this.options.defaultSort;
    this.sortOrder = this.options.defaultOrder;
    
    // Initialize component
    this.initialize();
  }
  
  /**
   * Initialize pagination component
   */
  initialize() {
    // Create component elements
    this.createElements();
    
    // Add event listeners
    this.setupEventListeners();
    
    // Load initial data
    this.loadPage(this.currentPage);
  }
  
  /**
   * Create pagination UI elements
   */
  createElements() {
    // Create container structure
    this.container.innerHTML = `
      <div class="pagination-controls">
        <div class="search-container ${!this.options.searchEnabled ? 'hidden' : ''}">
          <input type="text" id="search-input" placeholder="Search..." class="search-input">
          <button id="search-button" class="search-button">Search</button>
        </div>
        
        <div class="sort-container ${!this.options.sortEnabled ? 'hidden' : ''}">
          <select id="sort-select" class="sort-select">
            <option value="name">Name</option>
            <option value="molecular_weight">Molecular Weight</option>
            <option value="formula">Formula</option>
            <option value="created_at">Date Added</option>
          </select>
          <button id="sort-order" class="sort-order" data-order="${this.sortOrder}">
            ${this.sortOrder === 'asc' ? '↑' : '↓'}
          </button>
        </div>
      </div>
      
      <div id="content-container" class="content-container">
        <div id="loading-indicator" class="loading-indicator hidden">
          <div class="spinner"></div>
          <p>Loading...</p>
        </div>
        <div id="content-list" class="content-list"></div>
      </div>
      
      <div class="pagination-nav">
        <button id="prev-page" class="pagination-button" disabled>&lt; Previous</button>
        <span id="page-indicator" class="page-indicator">Page 1 of 1</span>
        <button id="next-page" class="pagination-button" disabled>Next &gt;</button>
      </div>
    `;
    
    // Get references to elements
    this.searchInput = document.getElementById('search-input');
    this.searchButton = document.getElementById('search-button');
    this.sortSelect = document.getElementById('sort-select');
    this.sortOrderButton = document.getElementById('sort-order');
    this.contentList = document.getElementById('content-list');
    this.loadingIndicator = document.getElementById('loading-indicator');
    this.prevButton = document.getElementById('prev-page');
    this.nextButton = document.getElementById('next-page');
    this.pageIndicator = document.getElementById('page-indicator');
    
    // Set initial sort field
    if (this.sortSelect) {
      this.sortSelect.value = this.sortField;
    }
  }
  
  /**
   * Set up event listeners for pagination controls
   */
  setupEventListeners() {
    // Search functionality
    if (this.searchButton) {
      this.searchButton.addEventListener('click', () => {
        this.searchTerm = this.searchInput.value;
        this.currentPage = 1;
        this.loadPage(this.currentPage);
      });
    }
    
    if (this.searchInput) {
      this.searchInput.addEventListener('keypress', (e) => {
        if (e.key === 'Enter') {
          this.searchTerm = this.searchInput.value;
          this.currentPage = 1;
          this.loadPage(this.currentPage);
        }
      });
    }
    
    // Sort functionality
    if (this.sortSelect) {
      this.sortSelect.addEventListener('change', () => {
        this.sortField = this.sortSelect.value;
        this.loadPage(this.currentPage);
      });
    }
    
    if (this.sortOrderButton) {
      this.sortOrderButton.addEventListener('click', () => {
        this.sortOrder = this.sortOrder === 'asc' ? 'desc' : 'asc';
        this.sortOrderButton.setAttribute('data-order', this.sortOrder);
        this.sortOrderButton.textContent = this.sortOrder === 'asc' ? '↑' : '↓';
        this.loadPage(this.currentPage);
      });
    }
    
    // Pagination navigation
    if (this.prevButton) {
      this.prevButton.addEventListener('click', () => {
        if (this.currentPage > 1) {
          this.loadPage(this.currentPage - 1);
        }
      });
    }
    
    if (this.nextButton) {
      this.nextButton.addEventListener('click', () => {
        if (this.currentPage < this.totalPages) {
          this.loadPage(this.currentPage + 1);
        }
      });
    }
  }
  
  /**
   * Load data for specified page
   * @param {number} page - Page number to load
   */
  loadPage(page) {
    // Prevent loading while already in progress
    if (this.isLoading) {
      return;
    }
    
    this.isLoading = true;
    this.showLoading(true);
    
    // Calculate offset
    const offset = (page - 1) * this.options.limit;
    
    // Build API URL with parameters
    let url = `${this.apiEndpoint}?limit=${this.options.limit}&offset=${offset}`;
    
    // Add search parameter if provided
    if (this.searchTerm) {
      url += `&search=${encodeURIComponent(this.searchTerm)}`;
    }
    
    // Add sort parameters
    url += `&sort=${this.sortField}&order=${this.sortOrder}`;
    
    // Fetch data from API
    fetch(url, {
      method: 'GET',
      headers: {
        'Content-Type': 'application/json'
      }
    })
    .then(response => {
      if (!response.ok) {
        throw new Error('Network response was not ok');
      }
      return response.json();
    })
    .then(data => {
      // Update state
      this.currentPage = page;
      this.totalPages = data.pagination.total_pages;
      
      // Update UI
      this.updatePageIndicator();
      this.updateNavigationButtons();
      
      // Render data
      if (typeof this.renderCallback === 'function') {
        this.renderCallback(data.data, this.contentList);
      } else {
        console.error('No render callback provided');
      }
    })
    .catch(error => {
      console.error('Error fetching data:', error);
      this.contentList.innerHTML = `
        <div class="error-message">
          <p>Failed to load data. Please try again.</p>
        </div>
      `;
    })
    .finally(() => {
      this.isLoading = false;
      this.showLoading(false);
    });
  }
  
  /**
   * Show or hide loading indicator
   * @param {boolean} show - Whether to show loading indicator
   */
  showLoading(show) {
    if (this.loadingIndicator) {
      if (show) {
        this.loadingIndicator.classList.remove('hidden');
      } else {
        this.loadingIndicator.classList.add('hidden');
      }
    }
  }
  
  /**
   * Update page indicator text
   */
  updatePageIndicator() {
    if (this.pageIndicator) {
      this.pageIndicator.textContent = `Page ${this.currentPage} of ${this.totalPages}`;
    }
  }
  
  /**
   * Update navigation button states
   */
  updateNavigationButtons() {
    if (this.prevButton) {
      this.prevButton.disabled = this.currentPage <= 1;
    }
    
    if (this.nextButton) {
      this.nextButton.disabled = this.currentPage >= this.totalPages;
    }
  }
}

// Export for use in other modules
if (typeof module !== 'undefined' && module.exports) {
  module.exports = { PaginationComponent };
}
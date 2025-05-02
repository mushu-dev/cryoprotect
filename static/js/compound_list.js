/**
 * CryoProtect Analyzer Web Interface
 * Compound List Module
 * 
 * This module handles the display and pagination of compound lists.
 * It integrates with the pagination component to provide a seamless
 * user experience for browsing large lists of compounds.
 */

// Initialize the compound list when the DOM is fully loaded
document.addEventListener('DOMContentLoaded', function() {
  initializeCompoundList();
});

/**
 * Initialize the compound list with pagination
 */
function initializeCompoundList() {
  const container = document.getElementById('molecules-container');
  if (!container) return;
  
  // Create a container for the pagination component
  const paginationContainer = document.createElement('div');
  paginationContainer.id = 'pagination-container';
  container.parentNode.insertBefore(paginationContainer, container);
  
  // Create a container for the molecules list
  const listContainer = document.createElement('div');
  listContainer.id = 'molecules-list-container';
  container.parentNode.insertBefore(listContainer, container);
  container.parentNode.removeChild(container);
  
  // Initialize the pagination component
  new PaginationComponent(
    'pagination-container',
    '/api/v1/molecules',
    renderMoleculesList,
    {
      limit: 12,
      searchEnabled: true,
      sortEnabled: true,
      defaultSort: 'name',
      defaultOrder: 'asc'
    }
  );
}

/**
 * Render a list of molecules
 * @param {Array} molecules - Array of molecule objects
 * @param {HTMLElement} container - Container element to render into
 */
function renderMoleculesList(molecules, container) {
  if (!container) return;
  
  // Clear the container
  container.innerHTML = '';
  
  // Create a row for the molecules grid
  const row = document.createElement('div');
  row.className = 'row';
  container.appendChild(row);
  
  if (!molecules || molecules.length === 0) {
    // Display a message if no molecules are found
    const emptyMessage = document.createElement('div');
    emptyMessage.className = 'col-12 text-center my-5';
    emptyMessage.innerHTML = `
      <div class="alert alert-info">
        <h4>No molecules found</h4>
        <p>Try adjusting your search criteria or import new molecules.</p>
        <button class="btn btn-primary" onclick="showImportModal()">
          <i class="bi bi-cloud-download"></i> Import Molecule
        </button>
      </div>
    `;
    row.appendChild(emptyMessage);
    return;
  }
  
  // Render each molecule as a card
  molecules.forEach(molecule => {
    const card = document.createElement('div');
    card.className = 'col-md-4 mb-4';
    card.innerHTML = `
      <div class="card molecule-card h-100">
        <div class="card-header">${molecule.name || `CID: ${molecule.cid}`}</div>
        <div class="card-body">
          <p><strong>Formula:</strong> ${molecule.molecular_formula || 'N/A'}</p>
          <p><strong>SMILES:</strong> ${molecule.smiles || 'N/A'}</p>
          ${molecule.properties ? renderMoleculeProperties(molecule.properties) : ''}
        </div>
        <div class="card-footer">
          <a href="/molecules/${molecule.id}" class="btn btn-primary">View Details</a>
        </div>
      </div>
    `;
    row.appendChild(card);
  });
}

/**
 * Render molecule properties
 * @param {Object} properties - Molecule properties
 * @returns {string} HTML string of rendered properties
 */
function renderMoleculeProperties(properties) {
  if (!properties) return '';
  
  // Select key properties to display in the card
  const keyProperties = [
    { key: 'molecular_weight', label: 'Molecular Weight', unit: 'g/mol' },
    { key: 'logp', label: 'LogP', unit: '' },
    { key: 'hydrogen_bond_donors', label: 'H-Bond Donors', unit: '' },
    { key: 'hydrogen_bond_acceptors', label: 'H-Bond Acceptors', unit: '' }
  ];
  
  let html = '<hr><h6>Properties:</h6><ul class="property-list">';
  
  keyProperties.forEach(prop => {
    if (properties[prop.key] !== undefined) {
      const value = typeof properties[prop.key] === 'number' 
        ? properties[prop.key].toFixed(2) 
        : properties[prop.key];
      html += `<li><span>${prop.label}:</span> ${value}${prop.unit ? ' ' + prop.unit : ''}</li>`;
    }
  });
  
  html += '</ul>';
  return html;
}

/**
 * Show the import molecule modal
 */
function showImportModal() {
  // Get the modal element
  const modal = document.getElementById('import-molecule-modal');
  if (!modal) {
    console.error('Import molecule modal not found');
    return;
  }
  
  // Show the modal
  const modalInstance = new bootstrap.Modal(modal);
  modalInstance.show();
}

/**
 * Import a molecule from PubChem
 * @param {Event} event - Form submit event
 */
function importMolecule(event) {
  event.preventDefault();
  
  const cidInput = document.getElementById('pubchem-cid');
  if (!cidInput) return;
  
  const cid = cidInput.value.trim();
  if (!cid) {
    alert('Please enter a valid PubChem CID');
    return;
  }
  
  // Show loading state
  const submitButton = event.target.querySelector('button[type="submit"]');
  const originalText = submitButton.innerHTML;
  submitButton.disabled = true;
  submitButton.innerHTML = '<span class="spinner-border spinner-border-sm" role="status" aria-hidden="true"></span> Importing...';
  
  // Call the API to import the molecule
  API.importMolecule(cid)
    .then(response => {
      // Hide the modal
      const modal = bootstrap.Modal.getInstance(document.getElementById('import-molecule-modal'));
      modal.hide();
      
      // Show success message
      showToast('Molecule Imported', `Successfully imported molecule: ${response.name || 'CID: ' + response.cid}`, 'success');
      
      // Refresh the molecules list
      initializeCompoundList();
    })
    .catch(error => {
      // Show error message
      ErrorHandler.handleApiError(error, {
        toastTitle: 'Import Error',
        showModal: false
      });
    })
    .finally(() => {
      // Reset the form
      submitButton.disabled = false;
      submitButton.innerHTML = originalText;
      cidInput.value = '';
    });
}

// Add CSS styles for the pagination component
function addPaginationStyles() {
  const style = document.createElement('style');
  style.textContent = `
    .pagination-controls {
      display: flex;
      justify-content: space-between;
      margin-bottom: 1rem;
      flex-wrap: wrap;
    }
    
    .search-container {
      display: flex;
      margin-bottom: 0.5rem;
    }
    
    .search-input {
      padding: 0.375rem 0.75rem;
      border: 1px solid #ced4da;
      border-radius: 0.25rem 0 0 0.25rem;
      flex-grow: 1;
    }
    
    .search-button {
      padding: 0.375rem 0.75rem;
      background-color: #007bff;
      color: white;
      border: 1px solid #007bff;
      border-radius: 0 0.25rem 0.25rem 0;
      cursor: pointer;
    }
    
    .sort-container {
      display: flex;
      margin-bottom: 0.5rem;
    }
    
    .sort-select {
      padding: 0.375rem 0.75rem;
      border: 1px solid #ced4da;
      border-radius: 0.25rem 0 0 0.25rem;
    }
    
    .sort-order {
      padding: 0.375rem 0.75rem;
      background-color: #6c757d;
      color: white;
      border: 1px solid #6c757d;
      border-radius: 0 0.25rem 0.25rem 0;
      cursor: pointer;
    }
    
    .content-container {
      position: relative;
      min-height: 200px;
    }
    
    .loading-indicator {
      position: absolute;
      top: 0;
      left: 0;
      right: 0;
      bottom: 0;
      background-color: rgba(255, 255, 255, 0.7);
      display: flex;
      flex-direction: column;
      justify-content: center;
      align-items: center;
      z-index: 10;
    }
    
    .spinner {
      width: 3rem;
      height: 3rem;
      border: 0.25rem solid #f3f3f3;
      border-top: 0.25rem solid #007bff;
      border-radius: 50%;
      animation: spin 1s linear infinite;
    }
    
    @keyframes spin {
      0% { transform: rotate(0deg); }
      100% { transform: rotate(360deg); }
    }
    
    .hidden {
      display: none !important;
    }
    
    .pagination-nav {
      display: flex;
      justify-content: center;
      align-items: center;
      margin-top: 1rem;
    }
    
    .pagination-button {
      padding: 0.375rem 0.75rem;
      background-color: #007bff;
      color: white;
      border: 1px solid #007bff;
      border-radius: 0.25rem;
      cursor: pointer;
      margin: 0 0.5rem;
    }
    
    .pagination-button:disabled {
      background-color: #6c757d;
      border-color: #6c757d;
      cursor: not-allowed;
    }
    
    .page-indicator {
      font-weight: bold;
    }
    
    .property-list {
      list-style: none;
      padding-left: 0;
      margin-bottom: 0;
    }
    
    .property-list li {
      margin-bottom: 0.25rem;
    }
    
    .property-list li span {
      font-weight: bold;
    }
  `;
  document.head.appendChild(style);
}

// Add the pagination styles when the script loads
addPaginationStyles();
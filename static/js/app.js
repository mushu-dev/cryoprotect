/**
 * CryoProtect Analyzer Web Interface
 * Main JavaScript file
 */

// Initialize the application when the DOM is fully loaded
document.addEventListener('DOMContentLoaded', function() {
  // Initialize Bootstrap tooltips and popovers
  initializeBootstrapComponents();
  
  // Initialize authentication
  initializeAuth();
  
  // Set up navigation event handlers
  setupNavigation();
  
  // Initialize page-specific functionality
  initializeCurrentPage();
});

/**
 * Initialize Bootstrap components
 */
function initializeBootstrapComponents() {
  // Initialize tooltips
  const tooltipTriggerList = [].slice.call(document.querySelectorAll('[data-bs-toggle="tooltip"]'));
  tooltipTriggerList.map(function(tooltipTriggerEl) {
    return new bootstrap.Tooltip(tooltipTriggerEl);
  });
  
  // Initialize popovers
  const popoverTriggerList = [].slice.call(document.querySelectorAll('[data-bs-toggle="popover"]'));
  popoverTriggerList.map(function(popoverTriggerEl) {
    return new bootstrap.Popover(popoverTriggerEl);
  });
}

/**
 * Initialize authentication
 */
function initializeAuth() {
  // Check if user is logged in
  Auth.checkSession()
    .then(user => {
      if (user) {
        // User is logged in
        updateUIForAuthenticatedUser(user);
      } else {
        // User is not logged in
        updateUIForUnauthenticatedUser();
      }
    })
    .catch(error => {
      console.error('Authentication error:', error);
      showToast('Authentication error', error.message, 'error');
    });
}

/**
 * Update UI for authenticated user
 */
function updateUIForAuthenticatedUser(user) {
  // Update navigation
  document.querySelectorAll('.nav-authenticated').forEach(el => {
    el.classList.remove('d-none');
  });
  
  document.querySelectorAll('.nav-unauthenticated').forEach(el => {
    el.classList.add('d-none');
  });
  
  // Update user info
  const userNameElements = document.querySelectorAll('.user-name');
  userNameElements.forEach(el => {
    el.textContent = user.email;
  });
  
  // Enable authenticated actions
  document.querySelectorAll('.requires-auth').forEach(el => {
    el.removeAttribute('disabled');
    el.classList.remove('disabled');
  });
}

/**
 * Update UI for unauthenticated user
 */
function updateUIForUnauthenticatedUser() {
  // Update navigation
  document.querySelectorAll('.nav-authenticated').forEach(el => {
    el.classList.add('d-none');
  });
  
  document.querySelectorAll('.nav-unauthenticated').forEach(el => {
    el.classList.remove('d-none');
  });
  
  // Disable authenticated actions
  document.querySelectorAll('.requires-auth').forEach(el => {
    el.setAttribute('disabled', 'disabled');
    el.classList.add('disabled');
  });
  
  // Redirect to login if on a protected page
  const currentPage = window.location.pathname;
  const protectedPages = [
    '/mixtures/new',
    '/predictions/new',
    '/experiments/new'
  ];
  
  if (protectedPages.some(page => currentPage.includes(page))) {
    window.location.href = '/login';
  }
}

/**
 * Set up navigation event handlers
 */
function setupNavigation() {
  // Handle navigation link clicks
  document.querySelectorAll('.nav-link').forEach(link => {
    link.addEventListener('click', function(event) {
      // Only handle links without href (dynamic navigation)
      if (!this.getAttribute('href')) {
        event.preventDefault();
        const target = this.getAttribute('data-target');
        if (target) {
          navigateTo(target);
        }
      }
    });
  });
  
  // Handle logout
  const logoutButton = document.getElementById('logout-button');
  if (logoutButton) {
    logoutButton.addEventListener('click', function(event) {
      event.preventDefault();
      Auth.signOut()
        .then(() => {
          showToast('Logged out', 'You have been successfully logged out', 'success');
          window.location.href = '/login';
        })
        .catch(error => {
          console.error('Logout error:', error);
          showToast('Logout error', error.message, 'error');
        });
    });
  }
}

/**
 * Navigate to a specific page
 */
function navigateTo(target) {
  window.location.href = target;
}

/**
 * Initialize page-specific functionality based on current page
 */
function initializeCurrentPage() {
  const currentPage = window.location.pathname;
  
  // Initialize page-specific functionality
  if (currentPage === '/' || currentPage === '/index.html') {
    initializeDashboard();
  } else if (currentPage.includes('/molecules')) {
    initializeMoleculesPage();
  } else if (currentPage.includes('/mixtures')) {
    initializeMixturesPage();
  } else if (currentPage.includes('/predictions')) {
    initializePredictionsPage();
  } else if (currentPage.includes('/experiments')) {
    initializeExperimentsPage();
  } else if (currentPage.includes('/comparisons')) {
    initializeComparisonsPage();
  } else if (currentPage.includes('/login')) {
    initializeLoginPage();
  }
}

/**
 * Initialize dashboard page
 */
function initializeDashboard() {
  // Load summary data
  Promise.all([
    API.getMolecules(),
    API.getMixtures()
  ])
    .then(([molecules, mixtures]) => {
      // Update dashboard stats
      updateDashboardStats(molecules, mixtures);
      
      // Initialize dashboard charts
      initializeDashboardCharts(molecules, mixtures);
    })
    .catch(error => {
      console.error('Error loading dashboard data:', error);
      showToast('Error', 'Failed to load dashboard data', 'error');
    });
}

/**
 * Initialize molecules page
 */
function initializeMoleculesPage() {
  // Load molecules data
  API.getMolecules()
    .then(molecules => {
      // Render molecules list or detail view
      const moleculeId = getIdFromUrl();
      if (moleculeId) {
        // Detail view for a specific molecule
        renderMoleculeDetail(moleculeId);
      } else {
        // List view for all molecules
        renderMoleculesList(molecules);
      }
    })
    .catch(error => {
      console.error('Error loading molecules:', error);
      showToast('Error', 'Failed to load molecules data', 'error');
    });
}

/**
 * Initialize mixtures page
 */
function initializeMixturesPage() {
  // Load mixtures data
  API.getMixtures()
    .then(mixtures => {
      // Render mixtures list or detail view
      const mixtureId = getIdFromUrl();
      if (mixtureId) {
        // Detail view for a specific mixture
        renderMixtureDetail(mixtureId);
      } else {
        // List view for all mixtures
        renderMixturesList(mixtures);
      }
    })
    .catch(error => {
      console.error('Error loading mixtures:', error);
      showToast('Error', 'Failed to load mixtures data', 'error');
    });
}

/**
 * Initialize predictions page
 */
function initializePredictionsPage() {
  const mixtureId = getIdFromUrl();
  if (!mixtureId) {
    showToast('Error', 'Mixture ID is required', 'error');
    return;
  }
  
  // Load mixture data and predictions
  Promise.all([
    API.getMixture(mixtureId),
    API.getPredictions(mixtureId)
  ])
    .then(([mixture, predictions]) => {
      // Update mixture info
      updateMixtureInfo(mixture);
      
      // Render predictions
      renderPredictions(predictions);
    })
    .catch(error => {
      console.error('Error loading predictions:', error);
      showToast('Error', 'Failed to load predictions data', 'error');
    });
}

/**
 * Initialize experiments page
 */
function initializeExperimentsPage() {
  const mixtureId = getIdFromUrl();
  if (!mixtureId) {
    showToast('Error', 'Mixture ID is required', 'error');
    return;
  }
  
  // Load mixture data and experiments
  Promise.all([
    API.getMixture(mixtureId),
    API.getExperiments(mixtureId)
  ])
    .then(([mixture, experiments]) => {
      // Update mixture info
      updateMixtureInfo(mixture);
      
      // Render experiments
      renderExperiments(experiments);
    })
    .catch(error => {
      console.error('Error loading experiments:', error);
      showToast('Error', 'Failed to load experiments data', 'error');
    });
}

/**
 * Initialize comparisons page
 */
function initializeComparisonsPage() {
  const mixtureId = getIdFromUrl();
  if (!mixtureId) {
    showToast('Error', 'Mixture ID is required', 'error');
    return;
  }
  
  // Load mixture data
  API.getMixture(mixtureId)
    .then(mixture => {
      // Update mixture info
      updateMixtureInfo(mixture);
      
      // Load property types for the dropdown
      return loadPropertyTypes('comparison-property');
    })
    .then(() => {
      // Set up comparison form
      setupComparisonForm(mixtureId);
    })
    .catch(error => {
      console.error('Error initializing comparisons page:', error);
      showToast('Error', 'Failed to initialize comparisons page', 'error');
    });
}

/**
 * Initialize login page
 */
function initializeLoginPage() {
  const loginForm = document.getElementById('login-form');
  if (loginForm) {
    loginForm.addEventListener('submit', function(event) {
      event.preventDefault();
      
      const email = document.getElementById('email').value;
      const password = document.getElementById('password').value;
      
      showSpinner();
      Auth.signIn(email, password)
        .then(user => {
          hideSpinner();
          showToast('Success', 'Logged in successfully', 'success');
          
          // Redirect to the dashboard or the page the user was trying to access
          const redirectUrl = new URLSearchParams(window.location.search).get('redirect') || '/';
          window.location.href = redirectUrl;
        })
        .catch(error => {
          hideSpinner();
          console.error('Login error:', error);
          showToast('Login failed', error.message, 'error');
        });
    });
  }
}

/**
 * Show a toast notification
 */
function showToast(title, message, type = 'info') {
  // Create toast container if it doesn't exist
  let toastContainer = document.querySelector('.toast-container');
  if (!toastContainer) {
    toastContainer = document.createElement('div');
    toastContainer.className = 'toast-container';
    document.body.appendChild(toastContainer);
  }
  
  // Create toast element
  const toastId = 'toast-' + Date.now();
  const toast = document.createElement('div');
  toast.className = `toast ${type}`;
  toast.id = toastId;
  toast.setAttribute('role', 'alert');
  toast.setAttribute('aria-live', 'assertive');
  toast.setAttribute('aria-atomic', 'true');
  
  // Set toast content
  toast.innerHTML = `
    <div class="toast-header">
      <strong class="me-auto">${title}</strong>
      <button type="button" class="btn-close toast-close-button" data-bs-dismiss="toast" aria-label="Close"></button>
    </div>
    <div class="toast-body">
      ${message}
    </div>
  `;
  
  // Add toast to container
  toastContainer.appendChild(toast);
  
  // Initialize and show the toast
  const bsToast = new bootstrap.Toast(toast, {
    autohide: true,
    delay: 5000
  });
  bsToast.show();
}

/**
 * Show loading spinner
 */
function showSpinner() {
  // Create spinner overlay if it doesn't exist
  let spinnerOverlay = document.querySelector('.spinner-overlay');
  if (!spinnerOverlay) {
    spinnerOverlay = document.createElement('div');
    spinnerOverlay.className = 'spinner-overlay';
    spinnerOverlay.innerHTML = `
      <div class="spinner-border" role="status">
        <span class="visually-hidden">Loading...</span>
      </div>
    `;
    document.body.appendChild(spinnerOverlay);
  }
  
  // Show spinner
  spinnerOverlay.style.display = 'flex';
}

/**
 * Hide loading spinner
 */
function hideSpinner() {
  const spinnerOverlay = document.querySelector('.spinner-overlay');
  if (spinnerOverlay) {
    spinnerOverlay.style.display = 'none';
  }
}

/**
 * Get ID from URL
 */
function getIdFromUrl() {
  const pathParts = window.location.pathname.split('/');
  for (let i = 0; i < pathParts.length; i++) {
    // Look for UUIDs in the URL
    if (/^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$/i.test(pathParts[i])) {
      return pathParts[i];
    }
  }
  return null;
}

// Helper functions for rendering UI components
function updateDashboardStats(molecules, mixtures) {
  document.getElementById('molecule-count')?.textContent = molecules.length;
  document.getElementById('mixture-count')?.textContent = mixtures.length;
}

function initializeDashboardCharts(molecules, mixtures) {
  Charts.createMoleculePropertiesChart('molecule-properties-chart', molecules);
  Charts.createMixtureCompositionChart('mixture-composition-chart', mixtures);
}

function renderMoleculeDetail(moleculeId) {
  API.getMolecule(moleculeId)
    .then(molecule => {
      document.getElementById('molecule-name')?.textContent = molecule.name || `CID: ${molecule.cid}`;
      document.getElementById('molecule-formula')?.textContent = molecule.molecular_formula || 'N/A';
      document.getElementById('molecule-smiles')?.textContent = molecule.smiles || 'N/A';
      
      // Render properties
      renderMoleculeProperties(molecule.properties);
    })
    .catch(error => {
      console.error('Error loading molecule detail:', error);
      showToast('Error', 'Failed to load molecule detail', 'error');
    });
}

function renderMoleculesList(molecules) {
  const container = document.getElementById('molecules-container');
  if (!container) return;
  
  container.innerHTML = '';
  
  molecules.forEach(molecule => {
    const card = document.createElement('div');
    card.className = 'col-md-4 mb-4';
    card.innerHTML = `
      <div class="card molecule-card h-100">
        <div class="card-header">${molecule.name || `CID: ${molecule.cid}`}</div>
        <div class="card-body">
          <p><strong>Formula:</strong> ${molecule.molecular_formula || 'N/A'}</p>
          <p><strong>SMILES:</strong> ${molecule.smiles || 'N/A'}</p>
        </div>
        <div class="card-footer">
          <a href="/molecules/${molecule.id}" class="btn btn-primary">View Details</a>
        </div>
      </div>
    `;
    container.appendChild(card);
  });
}

function renderMixtureDetail(mixtureId) {
  API.getMixture(mixtureId)
    .then(mixture => {
      document.getElementById('mixture-name')?.textContent = mixture.name;
      document.getElementById('mixture-description')?.textContent = mixture.description || 'No description';
      
      // Render components
      renderMixtureComponents(mixture.components);
    })
    .catch(error => {
      console.error('Error loading mixture detail:', error);
      showToast('Error', 'Failed to load mixture detail', 'error');
    });
}

function renderMixturesList(mixtures) {
  const container = document.getElementById('mixtures-container');
  if (!container) return;
  
  container.innerHTML = '';
  
  mixtures.forEach(mixture => {
    const card = document.createElement('div');
    card.className = 'col-md-4 mb-4';
    card.innerHTML = `
      <div class="card mixture-card h-100">
        <div class="card-header">${mixture.name}</div>
        <div class="card-body">
          <p>${mixture.description || 'No description'}</p>
          <p><strong>Components:</strong> ${mixture.components?.length || 0}</p>
        </div>
        <div class="card-footer">
          <a href="/mixtures/${mixture.id}" class="btn btn-primary">View Details</a>
        </div>
      </div>
    `;
    container.appendChild(card);
  });
}

function updateMixtureInfo(mixture) {
  document.getElementById('mixture-name-display')?.textContent = mixture.name;
  document.getElementById('mixture-description-display')?.textContent = mixture.description || 'No description';
  
  const componentsElement = document.getElementById('mixture-components-display');
  if (componentsElement && mixture.components) {
    componentsElement.innerHTML = '';
    mixture.components.forEach(component => {
      const item = document.createElement('li');
      item.className = 'list-group-item';
      item.textContent = `${component.name}: ${component.concentration} ${component.concentration_unit}`;
      componentsElement.appendChild(item);
    });
  }
}

function renderPredictions(predictions) {
  const container = document.getElementById('predictions-container');
  if (!container) return;
  
  container.innerHTML = '';
  
  if (predictions.length === 0) {
    container.innerHTML = '<div class="alert alert-info">No predictions found for this mixture.</div>';
    return;
  }
  
  predictions.forEach(prediction => {
    const card = document.createElement('div');
    card.className = 'card mb-3';
    card.innerHTML = `
      <div class="card-header">${prediction.property_name}</div>
      <div class="card-body">
        <p><strong>Value:</strong> ${getPredictionValue(prediction)}</p>
        <p><strong>Confidence:</strong> ${(prediction.confidence * 100).toFixed(1)}%</p>
        <p><strong>Method:</strong> ${prediction.calculation_method}</p>
      </div>
    `;
    container.appendChild(card);
  });
}

function renderExperiments(experiments) {
  const container = document.getElementById('experiments-container');
  if (!container) return;
  
  container.innerHTML = '';
  
  if (experiments.length === 0) {
    container.innerHTML = '<div class="alert alert-info">No experiments found for this mixture.</div>';
    return;
  }
  
  experiments.forEach(experiment => {
    const card = document.createElement('div');
    card.className = 'card mb-3';
    card.innerHTML = `
      <div class="card-header">${experiment.property_name}</div>
      <div class="card-body">
        <p><strong>Value:</strong> ${getExperimentValue(experiment)}</p>
        <p><strong>Conditions:</strong> ${experiment.experimental_conditions || 'Not specified'}</p>
        <p><strong>Date:</strong> ${experiment.date_performed}</p>
      </div>
    `;
    container.appendChild(card);
  });
}

function getPredictionValue(prediction) {
  if (prediction.numeric_value !== null) return prediction.numeric_value;
  if (prediction.text_value !== null) return prediction.text_value;
  if (prediction.boolean_value !== null) return prediction.boolean_value ? 'Yes' : 'No';
  return 'N/A';
}

function getExperimentValue(experiment) {
  if (experiment.numeric_value !== null) return experiment.numeric_value;
  if (experiment.text_value !== null) return experiment.text_value;
  if (experiment.boolean_value !== null) return experiment.boolean_value ? 'Yes' : 'No';
  return 'N/A';
}

function loadPropertyTypes(selectId) {
  const select = document.getElementById(selectId);
  if (!select) return Promise.reject(new Error(`Select element with ID ${selectId} not found`));
  
  return API.getPropertyTypes()
    .then(propertyTypes => {
      select.innerHTML = '<option value="">Select a property</option>';
      propertyTypes.forEach(property => {
        const option = document.createElement('option');
        option.value = property.name;
        option.textContent = property.name + (property.units ? ` (${property.units})` : '');
        select.appendChild(option);
      });
      return propertyTypes;
    });
}

function setupComparisonForm(mixtureId) {
  const form = document.getElementById('comparison-form');
  if (!form) return;
  
  form.addEventListener('submit', function(event) {
    event.preventDefault();
    
    const propertyName = document.getElementById('comparison-property').value;
    
    showSpinner();
    API.getComparison(mixtureId, propertyName)
      .then(comparison => {
        hideSpinner();
        renderComparisonResult(comparison);
      })
      .catch(error => {
        hideSpinner();
        console.error('Error getting comparison:', error);
        showToast('Error', `Failed to get comparison: ${error.message}`, 'error');
      });
  });
}

function renderComparisonResult(comparison) {
  const container = document.getElementById('comparison-result');
  if (!container) return;
  
  container.innerHTML = '';
  
  if (!comparison.prediction || !comparison.experiment) {
    container.innerHTML = '<div class="alert alert-warning">Cannot compare: missing prediction or experiment data.</div>';
    return;
  }
  
  const percentError = comparison.percent_error !== null ? comparison.percent_error.toFixed(2) + '%' : 'N/A';
  const difference = comparison.difference !== null ? comparison.difference.toFixed(2) : 'N/A';
  
  let resultClass = 'good';
  if (comparison.percent_error > 10) {
    resultClass = 'bad';
  } else if (comparison.percent_error > 5) {
    resultClass = 'warning';
  }
  
  container.innerHTML = `
    <div class="comparison-result ${resultClass}">
      <h4>Comparison Results</h4>
      <div class="row">
        <div class="col-md-6">
          <h5>Prediction</h5>
          <p><strong>Value:</strong> ${comparison.prediction.numeric_value}</p>
          <p><strong>Confidence:</strong> ${(comparison.prediction.confidence * 100).toFixed(1)}%</p>
          <p><strong>Method:</strong> ${comparison.prediction.method}</p>
        </div>
        <div class="col-md-6">
          <h5>Experiment</h5>
          <p><strong>Value:</strong> ${comparison.experiment.numeric_value}</p>
          <p><strong>Conditions:</strong> ${comparison.experiment.conditions || 'Not specified'}</p>
          <p><strong>Date:</strong> ${comparison.experiment.date}</p>
        </div>
      </div>
      <hr>
      <div class="row">
        <div class="col-md-6">
          <p class="difference"><strong>Difference:</strong> ${difference}</p>
        </div>
        <div class="col-md-6">
          <p class="percent-error"><strong>Percent Error:</strong> ${percentError}</p>
        </div>
      </div>
    </div>
  `;
}

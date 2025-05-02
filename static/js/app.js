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
      ErrorHandler.handleApiError(error, {
        toastTitle: 'Authentication Error',
        redirectOnAuthError: false
      });
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
          ErrorHandler.handleApiError(error, {
            toastTitle: 'Logout Error'
          });
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
      ErrorHandler.handleApiError(error, {
        toastTitle: 'Dashboard Error',
        modalTitle: 'Error Loading Dashboard',
        showModal: true
      });
      
      // Show empty state for dashboard
      document.getElementById('dashboard-content')?.classList.add('dashboard-error');
      document.getElementById('dashboard-error-message')?.classList.remove('d-none');
    });
}

/**
 * Initialize molecules page
 */
function initializeMoleculesPage() {
  // Check if we're on a detail page or list page
  const moleculeId = getIdFromUrl();
  if (moleculeId) {
    // Detail view for a specific molecule
    renderMoleculeDetail(moleculeId);
  } else {
    // List view is handled by compound_list.js
    // The initialization happens in the DOMContentLoaded event in that file
    // Just ensure we have the container ready
    const container = document.getElementById('molecules-container');
    if (!container) {
      const mainContent = document.querySelector('.main-content');
      if (mainContent) {
        const newContainer = document.createElement('div');
        newContainer.id = 'molecules-container';
        newContainer.className = 'container-fluid';
        mainContent.appendChild(newContainer);
      }
    }
  }
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
      ErrorHandler.handleApiError(error, {
        toastTitle: 'Mixtures Error',
        showModal: false
      });
      
      // Show empty state for mixtures
      const container = document.getElementById('mixtures-container');
      if (container) {
        container.innerHTML = `
          <div class="alert alert-danger" role="alert">
            <h4 class="alert-heading"><i class="bi bi-exclamation-triangle-fill me-2"></i>Error Loading Mixtures</h4>
            <p>We encountered a problem while loading the mixtures data. Please try again later.</p>
            <hr>
            <p class="mb-0">
              <button type="button" class="btn btn-outline-danger" onclick="initializeMixturesPage()">
                <i class="bi bi-arrow-clockwise me-1"></i> Retry
              </button>
            </p>
          </div>
        `;
      }
    });
}

/**
 * Initialize predictions page
 */
function initializePredictionsPage() {
  const mixtureId = getIdFromUrl();
  if (!mixtureId) {
    showErrorModal('Missing Information', 'Mixture ID is required to view predictions. Please select a mixture first.');
    // Redirect to mixtures page after a delay
    setTimeout(() => {
      window.location.href = '/mixtures';
    }, 3000);
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
      ErrorHandler.handleApiError(error, {
        toastTitle: 'Predictions Error',
        showModal: false
      });
      
      // Show empty state for predictions
      const container = document.getElementById('predictions-container');
      if (container) {
        container.innerHTML = `
          <div class="alert alert-danger" role="alert">
            <h4 class="alert-heading"><i class="bi bi-exclamation-triangle-fill me-2"></i>Error Loading Predictions</h4>
            <p>We encountered a problem while loading the predictions data. Please try again later.</p>
            <hr>
            <p class="mb-0">
              <button type="button" class="btn btn-outline-danger" onclick="initializePredictionsPage()">
                <i class="bi bi-arrow-clockwise me-1"></i> Retry
              </button>
            </p>
          </div>
        `;
      }
    });
}

/**
 * Initialize experiments page
 */
function initializeExperimentsPage() {
  const mixtureId = getIdFromUrl();
  if (!mixtureId) {
    showErrorModal('Missing Information', 'Mixture ID is required to view experiments. Please select a mixture first.');
    // Redirect to mixtures page after a delay
    setTimeout(() => {
      window.location.href = '/mixtures';
    }, 3000);
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
      ErrorHandler.handleApiError(error, {
        toastTitle: 'Experiments Error',
        showModal: false
      });
      
      // Show empty state for experiments
      const container = document.getElementById('experiments-container');
      if (container) {
        container.innerHTML = `
          <div class="alert alert-danger" role="alert">
            <h4 class="alert-heading"><i class="bi bi-exclamation-triangle-fill me-2"></i>Error Loading Experiments</h4>
            <p>We encountered a problem while loading the experiments data. Please try again later.</p>
            <hr>
            <p class="mb-0">
              <button type="button" class="btn btn-outline-danger" onclick="initializeExperimentsPage()">
                <i class="bi bi-arrow-clockwise me-1"></i> Retry
              </button>
            </p>
          </div>
        `;
      }
    });
}

/**
 * Initialize comparisons page
 */
function initializeComparisonsPage() {
  const mixtureId = getIdFromUrl();
  if (!mixtureId) {
    showErrorModal('Missing Information', 'Mixture ID is required to view comparisons. Please select a mixture first.');
    // Redirect to mixtures page after a delay
    setTimeout(() => {
      window.location.href = '/mixtures';
    }, 3000);
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
      ErrorHandler.handleApiError(error, {
        toastTitle: 'Comparisons Error',
        showModal: true,
        modalTitle: 'Error Loading Comparisons'
      });
      
      // Show empty state for comparisons
      const container = document.getElementById('comparisons-container');
      if (container) {
        container.innerHTML = `
          <div class="alert alert-danger" role="alert">
            <h4 class="alert-heading"><i class="bi bi-exclamation-triangle-fill me-2"></i>Error Loading Comparisons</h4>
            <p>We encountered a problem while loading the comparisons data. Please try again later.</p>
            <hr>
            <p class="mb-0">
              <button type="button" class="btn btn-outline-danger" onclick="initializeComparisonsPage()">
                <i class="bi bi-arrow-clockwise me-1"></i> Retry
              </button>
            </p>
          </div>
        `;
      }
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
      
      // Show loading state
      ErrorHandler.setFormLoading(loginForm, true);
      
      Auth.signIn(email, password)
        .then(user => {
          // Hide loading state
          ErrorHandler.setFormLoading(loginForm, false);
          
          // Show success message
          showSuccessToast('Login Successful', 'You have been logged in successfully');
          
          // Announce to screen readers
          announceToScreenReader('Login successful. Redirecting to dashboard.');
          
          // Redirect to the dashboard or the page the user was trying to access
          const redirectUrl = new URLSearchParams(window.location.search).get('redirect') || '/';
          window.location.href = redirectUrl;
        })
        .catch(error => {
          // Hide loading state
          ErrorHandler.setFormLoading(loginForm, false);
          
          // Handle the error
          ErrorHandler.handleFormError(loginForm, error);
          
          // Focus back on the email field for retry
          document.getElementById('email')?.focus();
        });
    });
  }
}

/**
 * Show a toast notification
 * @param {string} title - The toast title
 * @param {string} message - The toast message
 * @param {string} type - The toast type (success, error, warning, info)
 * @param {Object} options - Additional options
 */
function showToast(title, message, type = 'info', options = {}) {
  // Get toast container
  const toastContainer = document.getElementById('toast-container');
  if (!toastContainer) {
    console.error('Toast container not found');
    return;
  }
  
  // Create toast element
  const toastId = 'toast-' + Date.now();
  const toast = document.createElement('div');
  toast.className = 'toast';
  toast.id = toastId;
  toast.setAttribute('role', 'status');
  toast.setAttribute('aria-live', type === 'error' ? 'assertive' : 'polite');
  toast.setAttribute('aria-atomic', 'true');
  
  // Set toast content with appropriate icon and color
  let iconClass = '';
  let bgClass = '';
  
  switch (type) {
    case 'success':
      iconClass = 'bi-check-circle-fill text-success';
      bgClass = 'border-success';
      break;
    case 'error':
      iconClass = 'bi-exclamation-triangle-fill text-danger';
      bgClass = 'border-danger';
      break;
    case 'warning':
      iconClass = 'bi-exclamation-circle-fill text-warning';
      bgClass = 'border-warning';
      break;
    case 'info':
    default:
      iconClass = 'bi-info-circle-fill text-info';
      bgClass = 'border-info';
      break;
  }
  
  toast.innerHTML = `
    <div class="toast-header ${bgClass}">
      <i class="bi ${iconClass} me-2"></i>
      <strong class="me-auto">${title}</strong>
      <button type="button" class="btn-close" data-bs-dismiss="toast" aria-label="Close"></button>
    </div>
    <div class="toast-body">
      ${message}
    </div>
  `;
  
  // Add toast to container
  toastContainer.appendChild(toast);
  
  // Initialize and show the toast
  const bsToast = new bootstrap.Toast(toast, {
    autohide: options.autohide !== false,
    delay: options.delay || 5000
  });
  
  bsToast.show();
  
  // Announce to screen readers
  announceToScreenReader(message, type);
  
  // Return the toast instance for potential further manipulation
  return { element: toast, instance: bsToast };
}

/**
 * Show a success toast
 * @param {string} title - The toast title
 * @param {string} message - The toast message
 * @param {Object} options - Additional options
 */
function showSuccessToast(title, message, options = {}) {
  return showToast(title, message, 'success', options);
}

/**
 * Show an error toast
 * @param {string} title - The toast title
 * @param {string} message - The toast message
 * @param {Object} options - Additional options
 */
function showErrorToast(title, message, options = {}) {
  return showToast(title, message, 'error', options);
}

/**
 * Show a warning toast
 * @param {string} title - The toast title
 * @param {string} message - The toast message
 * @param {Object} options - Additional options
 */
function showWarningToast(title, message, options = {}) {
  return showToast(title, message, 'warning', options);
}

/**
 * Show an info toast
 * @param {string} title - The toast title
 * @param {string} message - The toast message
 * @param {Object} options - Additional options
 */
function showInfoToast(title, message, options = {}) {
  return showToast(title, message, 'info', options);
}

/**
 * Show a modal dialog
 * @param {string} title - The modal title
 * @param {string} message - The modal message
 * @param {Object} options - Additional options
 */
function showModal(title, message, options = {}) {
  // Get modal elements
  const modal = document.getElementById('feedback-modal');
  const modalTitle = document.getElementById('feedback-modal-title');
  const modalBody = document.getElementById('feedback-modal-body');
  const modalAction = document.getElementById('feedback-modal-action');
  
  if (!modal || !modalTitle || !modalBody || !modalAction) {
    console.error('Modal elements not found');
    return;
  }
  
  // Set modal content
  modalTitle.textContent = title;
  modalBody.innerHTML = message;
  
  // Configure action button
  if (options.actionText) {
    modalAction.textContent = options.actionText;
    modalAction.style.display = 'block';
    
    // Remove previous event listeners
    const newModalAction = modalAction.cloneNode(true);
    modalAction.parentNode.replaceChild(newModalAction, modalAction);
    
    // Add new event listener
    if (options.onAction && typeof options.onAction === 'function') {
      newModalAction.addEventListener('click', function() {
        options.onAction();
        bootstrap.Modal.getInstance(modal).hide();
      });
    }
  } else {
    modalAction.style.display = 'none';
  }
  
  // Set modal type/style
  const modalDialog = modal.querySelector('.modal-dialog');
  modalDialog.className = 'modal-dialog modal-dialog-centered';
  
  if (options.size) {
    modalDialog.classList.add(`modal-${options.size}`);
  }
  
  if (options.type) {
    const headerClass = `bg-${options.type} text-white`;
    modal.querySelector('.modal-header').className = `modal-header ${headerClass}`;
  } else {
    modal.querySelector('.modal-header').className = 'modal-header';
  }
  
  // Show the modal
  const modalInstance = new bootstrap.Modal(modal);
  modalInstance.show();
  
  // Set focus to the modal when it's shown
  modal.addEventListener('shown.bs.modal', function() {
    modalTitle.focus();
  }, { once: true });
  
  // Announce to screen readers
  announceToScreenReader(message, 'modal');
  
  // Return the modal instance for potential further manipulation
  return { element: modal, instance: modalInstance };
}

/**
 * Show a confirmation modal
 * @param {string} title - The modal title
 * @param {string} message - The modal message
 * @param {Function} onConfirm - Function to call when confirmed
 * @param {Object} options - Additional options
 */
function showConfirmModal(title, message, onConfirm, options = {}) {
  return showModal(title, message, {
    type: options.type || 'primary',
    actionText: options.actionText || 'Confirm',
    onAction: onConfirm,
    size: options.size
  });
}

/**
 * Show an error modal
 * @param {string} title - The modal title
 * @param {string} message - The modal message
 * @param {Object} options - Additional options
 */
function showErrorModal(title, message, options = {}) {
  return showModal(title, message, {
    type: 'danger',
    actionText: options.actionText || 'OK',
    onAction: options.onAction,
    size: options.size
  });
}

/**
 * Announce a message to screen readers
 * @param {string} message - The message to announce
 * @param {string} type - The type of message
 */
function announceToScreenReader(message, type = 'info') {
  const srAnnouncements = document.getElementById('sr-announcements');
  if (!srAnnouncements) return;
  
  // Clear previous announcements
  srAnnouncements.textContent = '';
  
  // Add prefix based on type
  let prefix = '';
  switch (type) {
    case 'error':
      prefix = 'Error: ';
      break;
    case 'warning':
      prefix = 'Warning: ';
      break;
    case 'success':
      prefix = 'Success: ';
      break;
    case 'modal':
      prefix = 'Dialog opened: ';
      break;
    default:
      prefix = '';
  }
  
  // Set the announcement text
  srAnnouncements.textContent = prefix + message;
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
      
      // Set PubChem link
      const pubchemLink = document.getElementById('pubchem-link');
      if (pubchemLink && molecule.cid) {
        pubchemLink.href = `https://pubchem.ncbi.nlm.nih.gov/compound/${molecule.cid}`;
      }
      
      // Set integrated viewer link
      const integratedViewerLink = document.getElementById('integrated-viewer-link');
      if (integratedViewerLink) {
        integratedViewerLink.href = `/molecules/integrated?id=${moleculeId}&smiles=${encodeURIComponent(molecule.smiles)}`;
      }
      
      // Render properties
      renderMoleculeProperties(molecule.properties);
    })
    .catch(error => {
      console.error('Error loading molecule detail:', error);
      showToast('Error', 'Failed to load molecule detail', 'error');
    });
}

// Note: renderMoleculesList function has been moved to compound_list.js
// This comment is kept here to maintain the file structure and prevent confusion

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
        <div class="lab-verification-section mt-3" data-experiment-id="${experiment.id}"></div>
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

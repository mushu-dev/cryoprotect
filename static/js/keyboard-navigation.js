/**
 * CryoProtect Analyzer - Keyboard Navigation Accessibility
 * This file contains functions to improve keyboard navigation throughout the application
 */

document.addEventListener('DOMContentLoaded', function() {
  // Initialize focus management for all interactive components
  initFocusManagement();
});

/**
 * Initialize focus management for interactive components
 */
function initFocusManagement() {
  // Initialize focus management for modals
  initModalFocusManagement();
  
  // Initialize focus management for collapsible sections
  initCollapsibleFocusManagement();
  
  // Initialize focus management for tabs
  initTabFocusManagement();
  
  // Initialize focus management for view toggling
  initViewToggleFocusManagement();
}

/**
 * Initialize focus management for Bootstrap modals
 */
function initModalFocusManagement() {
  // Get all modals
  const modals = document.querySelectorAll('.modal');
  
  // Store the element that had focus before the modal was opened
  let previouslyFocusedElement = null;
  
  modals.forEach(modal => {
    // Store the trigger elements for each modal
    const modalTriggers = document.querySelectorAll(`[data-bs-toggle="modal"][data-bs-target="#${modal.id}"], [data-bs-toggle="modal"][href="#${modal.id}"]`);
    
    // When modal is about to be shown, store the currently focused element
    modal.addEventListener('show.bs.modal', function() {
      previouslyFocusedElement = document.activeElement;
    });
    
    // When modal is shown, focus the first focusable element
    modal.addEventListener('shown.bs.modal', function() {
      const focusableElements = getFocusableElements(modal);
      if (focusableElements.length > 0) {
        focusableElements[0].focus();
      }
      
      // Add a click event listener to the document to prevent clicks outside the modal
      document.addEventListener('click', handleOutsideClick);
    });
    
    // Trap focus within the modal when it's open
    modal.addEventListener('keydown', function(event) {
      // Handle Escape key to close the modal if not handled by Bootstrap
      if (event.key === 'Escape') {
        const bsModal = bootstrap.Modal.getInstance(modal);
        if (bsModal) {
          bsModal.hide();
        }
      }
      
      // Handle Tab key navigation
      trapFocus(event, modal);
    });
    
    // When modal is hidden, restore focus to the previously focused element
    modal.addEventListener('hidden.bs.modal', function() {
      if (previouslyFocusedElement) {
        previouslyFocusedElement.focus();
      } else if (modalTriggers.length > 0) {
        // If we don't have a previously focused element, focus the trigger
        modalTriggers[0].focus();
      }
      
      // Remove the click event listener from the document
      document.removeEventListener('click', handleOutsideClick);
    });
    
    // Function to handle clicks outside the modal
    function handleOutsideClick(event) {
      // If the modal is open and the click is outside the modal content
      if (modal.classList.contains('show') && !modal.querySelector('.modal-content').contains(event.target)) {
        // Refocus the first focusable element in the modal
        const focusableElements = getFocusableElements(modal);
        if (focusableElements.length > 0) {
          event.preventDefault();
          focusableElements[0].focus();
        }
      }
    }
  });
}

/**
 * Initialize focus management for collapsible sections
 */
function initCollapsibleFocusManagement() {
  // Get all collapsible triggers
  const collapseTriggers = document.querySelectorAll('[data-bs-toggle="collapse"]');
  
  collapseTriggers.forEach(trigger => {
    const targetId = trigger.getAttribute('data-bs-target') || trigger.getAttribute('href');
    if (!targetId) return;
    
    const target = document.querySelector(targetId);
    if (!target) return;
    
    // When collapse is shown, focus the first focusable element
    target.addEventListener('shown.bs.collapse', function() {
      const focusableElements = getFocusableElements(target);
      if (focusableElements.length > 0) {
        focusableElements[0].focus();
      }
    });
    
    // When collapse is hidden, return focus to the trigger
    target.addEventListener('hidden.bs.collapse', function() {
      trigger.focus();
    });
  });
}

/**
 * Initialize focus management for Bootstrap tabs
 */
function initTabFocusManagement() {
  // Get all tab buttons
  const tabButtons = document.querySelectorAll('[data-bs-toggle="tab"]');
  
  tabButtons.forEach(button => {
    const targetId = button.getAttribute('data-bs-target') || button.getAttribute('href');
    if (!targetId) return;
    
    // When tab is shown, focus the first focusable element in the tab panel
    button.addEventListener('shown.bs.tab', function() {
      const tabPanel = document.querySelector(targetId);
      if (!tabPanel) return;
      
      const focusableElements = getFocusableElements(tabPanel);
      if (focusableElements.length > 0) {
        focusableElements[0].focus();
      } else {
        // If no focusable elements, focus the tab panel itself
        tabPanel.setAttribute('tabindex', '-1');
        tabPanel.focus();
      }
    });
  });
}

/**
 * Initialize focus management for view toggling
 */
function initViewToggleFocusManagement() {
  // Handle view toggling for mixtures page
  const addExperimentBtn = document.getElementById('record-experiment-btn');
  const backToExperimentsBtn = document.getElementById('back-to-experiments-btn');
  
  if (addExperimentBtn) {
    addExperimentBtn.addEventListener('click', function(event) {
      event.preventDefault();
      const experimentsListView = document.getElementById('experiments-list-view');
      const addExperimentView = document.getElementById('add-experiment-view');
      
      if (experimentsListView && addExperimentView) {
        experimentsListView.style.display = 'none';
        addExperimentView.style.display = 'block';
        
        // Focus the first focusable element in the add experiment view
        const focusableElements = getFocusableElements(addExperimentView);
        if (focusableElements.length > 0) {
          focusableElements[0].focus();
        }
      }
    });
  }
  
  if (backToExperimentsBtn) {
    backToExperimentsBtn.addEventListener('click', function(event) {
      event.preventDefault();
      const experimentsListView = document.getElementById('experiments-list-view');
      const addExperimentView = document.getElementById('add-experiment-view');
      
      if (experimentsListView && addExperimentView) {
        experimentsListView.style.display = 'block';
        addExperimentView.style.display = 'none';
        
        // Focus the "Record Experiment" button
        if (addExperimentBtn) {
          addExperimentBtn.focus();
        }
      }
    });
  }
  
  // Handle view toggling for predictions page
  const addPredictionBtn = document.getElementById('add-prediction-btn');
  const backToPredictionsBtn = document.getElementById('back-to-predictions-btn');
  
  if (addPredictionBtn) {
    addPredictionBtn.addEventListener('click', function(event) {
      event.preventDefault();
      const predictionsListView = document.getElementById('predictions-list-view');
      const addPredictionView = document.getElementById('add-prediction-view');
      
      if (predictionsListView && addPredictionView) {
        predictionsListView.style.display = 'none';
        addPredictionView.style.display = 'block';
        
        // Focus the first focusable element in the add prediction view
        const focusableElements = getFocusableElements(addPredictionView);
        if (focusableElements.length > 0) {
          focusableElements[0].focus();
        }
      }
    });
  }
  
  if (backToPredictionsBtn) {
    backToPredictionsBtn.addEventListener('click', function(event) {
      event.preventDefault();
      const predictionsListView = document.getElementById('predictions-list-view');
      const addPredictionView = document.getElementById('add-prediction-view');
      
      if (predictionsListView && addPredictionView) {
        predictionsListView.style.display = 'block';
        addPredictionView.style.display = 'none';
        
        // Focus the "Add Prediction" button
        if (addPredictionBtn) {
          addPredictionBtn.focus();
        }
      }
    });
  }
}

/**
 * Get all focusable elements within a container
 * @param {HTMLElement} container - The container element
 * @returns {Array} - Array of focusable elements
 */
function getFocusableElements(container) {
  return Array.from(container.querySelectorAll(
    'button, [href], input, select, textarea, [tabindex]:not([tabindex="-1"])'
  )).filter(el => !el.hasAttribute('disabled') && !el.getAttribute('aria-hidden'));
}

/**
 * Trap focus within a container (for modals, dialogs, etc.)
 * @param {Event} event - The keydown event
 * @param {HTMLElement} container - The container to trap focus within
 */
function trapFocus(event, container) {
  // Only process if Tab key is pressed
  if (event.key !== 'Tab') return;
  
  const focusableElements = getFocusableElements(container);
  if (focusableElements.length === 0) return;
  
  const firstElement = focusableElements[0];
  const lastElement = focusableElements[focusableElements.length - 1];
  
  // If Shift+Tab on first element, move to last element
  if (event.shiftKey && document.activeElement === firstElement) {
    event.preventDefault();
    lastElement.focus();
  } 
  // If Tab on last element, move to first element
  else if (!event.shiftKey && document.activeElement === lastElement) {
    event.preventDefault();
    firstElement.focus();
  }
}

/**
 * Add keyboard support for custom interactive components
 * @param {HTMLElement} element - The interactive element
 * @param {Function} clickHandler - The click handler function
 */
function addKeyboardSupport(element, clickHandler) {
  if (!element || typeof clickHandler !== 'function') return;
  
  element.setAttribute('tabindex', '0');
  element.setAttribute('role', 'button');
  
  element.addEventListener('keydown', function(event) {
    // Trigger click on Enter or Space
    if (event.key === 'Enter' || event.key === ' ') {
      event.preventDefault();
      clickHandler.call(this, event);
    }
  });
}
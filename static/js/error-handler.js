/**
 * CryoProtect Analyzer - Error Handler
 * This file contains functions to handle errors and provide consistent feedback
 */

const ErrorHandler = (function() {
  /**
   * Handle API errors and display appropriate feedback
   * @param {Error} error - The error object
   * @param {Object} options - Additional options
   */
  function handleApiError(error, options = {}) {
    console.error('API Error:', error);
    
    // Default options
    const defaultOptions = {
      showToast: true,
      showModal: false,
      toastTitle: 'Error',
      modalTitle: 'Error',
      redirectOnAuthError: true,
      retryCallback: null
    };
    
    // Merge options
    const settings = { ...defaultOptions, ...options };
    
    // Check if it's an API error
    if (error.isApiError) {
      switch (error.errorType) {
        case 'AUTHENTICATION_ERROR':
          // Handle authentication errors
          if (settings.showModal) {
            showErrorModal(settings.modalTitle, error.userMessage, {
              actionText: 'Go to Login',
              onAction: () => {
                window.location.href = '/login?redirect=' + encodeURIComponent(window.location.pathname);
              }
            });
          } else if (settings.showToast) {
            showErrorToast(settings.toastTitle, error.userMessage);
          }
          
          // Redirect to login page if needed
          if (settings.redirectOnAuthError) {
            setTimeout(() => {
              window.location.href = '/login?redirect=' + encodeURIComponent(window.location.pathname);
            }, 2000);
          }
          break;
          
        case 'AUTHORIZATION_ERROR':
          // Handle authorization errors
          if (settings.showModal) {
            showErrorModal(settings.modalTitle, error.userMessage);
          } else if (settings.showToast) {
            showErrorToast(settings.toastTitle, error.userMessage);
          }
          break;
          
        case 'VALIDATION_ERROR':
          // Handle validation errors
          if (error.validationErrors) {
            // If we have a form, highlight the fields with errors
            if (options.form) {
              highlightFormErrors(options.form, error.validationErrors);
            }
            
            // Show detailed validation errors
            if (settings.showModal) {
              showErrorModal(settings.modalTitle, error.userMessage);
            } else if (settings.showToast) {
              showErrorToast(settings.toastTitle, error.userMessage);
            }
          } else {
            // Generic validation error
            if (settings.showModal) {
              showErrorModal(settings.modalTitle, error.userMessage);
            } else if (settings.showToast) {
              showErrorToast(settings.toastTitle, error.userMessage);
            }
          }
          break;
          
        case 'NETWORK_ERROR':
          // Handle network errors
          if (settings.showModal) {
            showErrorModal(settings.modalTitle, error.userMessage, {
              actionText: settings.retryCallback ? 'Retry' : null,
              onAction: settings.retryCallback
            });
          } else if (settings.showToast) {
            showErrorToast(settings.toastTitle, error.userMessage, {
              autohide: false
            });
          }
          break;
          
        case 'SERVER_ERROR':
          // Handle server errors
          if (settings.showModal) {
            showErrorModal(settings.modalTitle, error.userMessage);
          } else if (settings.showToast) {
            showErrorToast(settings.toastTitle, error.userMessage);
          }
          break;
          
        default:
          // Handle other API errors
          if (settings.showModal) {
            showErrorModal(settings.modalTitle, error.userMessage || error.message);
          } else if (settings.showToast) {
            showErrorToast(settings.toastTitle, error.userMessage || error.message);
          }
          break;
      }
    } else {
      // Handle generic errors
      const message = error.message || 'An unexpected error occurred';
      
      if (settings.showModal) {
        showErrorModal(settings.modalTitle, message);
      } else if (settings.showToast) {
        showErrorToast(settings.toastTitle, message);
      }
    }
  }
  
  /**
   * Highlight form fields with validation errors
   * @param {HTMLFormElement} form - The form element
   * @param {Object} validationErrors - Validation errors object
   */
  function highlightFormErrors(form, validationErrors) {
    // Clear previous errors
    clearFormErrors(form);
    
    // Add errors to fields
    for (const field in validationErrors) {
      if (Object.prototype.hasOwnProperty.call(validationErrors, field)) {
        const fieldErrors = validationErrors[field];
        const fieldElement = form.querySelector(`[name="${field}"]`);
        
        if (fieldElement) {
          // Mark field as invalid
          fieldElement.classList.add('is-invalid');
          fieldElement.setAttribute('aria-invalid', 'true');
          
          // Add error message
          const errorContainer = document.getElementById(`${fieldElement.id}-error`);
          if (errorContainer) {
            errorContainer.textContent = Array.isArray(fieldErrors) ? fieldErrors.join(', ') : fieldErrors;
          }
        }
      }
    }
    
    // Focus the first invalid field
    const firstInvalidField = form.querySelector('.is-invalid');
    if (firstInvalidField) {
      firstInvalidField.focus();
    }
    
    // Announce validation errors to screen readers
    announceValidationErrors(form);
  }
  
  /**
   * Handle form submission errors
   * @param {HTMLFormElement} form - The form element
   * @param {Error} error - The error object
   */
  function handleFormError(form, error) {
    // If it's an API error with validation errors
    if (error.isApiError && error.validationErrors) {
      highlightFormErrors(form, error.validationErrors);
    } else {
      // Generic form error
      showErrorToast('Form Error', error.message || 'There was a problem submitting the form');
    }
  }
  
  /**
   * Show a loading indicator for a form
   * @param {HTMLFormElement} form - The form element
   * @param {boolean} isLoading - Whether the form is loading
   */
  function setFormLoading(form, isLoading) {
    const submitButton = form.querySelector('button[type="submit"]');
    
    if (isLoading) {
      // Disable form
      form.classList.add('is-loading');
      
      // Disable submit button and show spinner
      if (submitButton) {
        submitButton.disabled = true;
        
        // Save original text
        submitButton.dataset.originalText = submitButton.innerHTML;
        
        // Show spinner
        submitButton.innerHTML = `
          <span class="spinner-border spinner-border-sm" role="status" aria-hidden="true"></span>
          <span>Processing...</span>
        `;
      }
      
      // Announce to screen readers
      announceToScreenReader('Form is processing, please wait...');
    } else {
      // Enable form
      form.classList.remove('is-loading');
      
      // Enable submit button and restore text
      if (submitButton) {
        submitButton.disabled = false;
        
        // Restore original text
        if (submitButton.dataset.originalText) {
          submitButton.innerHTML = submitButton.dataset.originalText;
        }
      }
    }
  }
  
  /**
   * Show a loading indicator for a button
   * @param {HTMLButtonElement} button - The button element
   * @param {boolean} isLoading - Whether the button is loading
   */
  function setButtonLoading(button, isLoading) {
    if (isLoading) {
      // Disable button
      button.disabled = true;
      
      // Save original text
      button.dataset.originalText = button.innerHTML;
      
      // Show spinner
      button.innerHTML = `
        <span class="spinner-border spinner-border-sm" role="status" aria-hidden="true"></span>
        <span>Processing...</span>
      `;
      
      // Announce to screen readers
      announceToScreenReader('Processing, please wait...');
    } else {
      // Enable button
      button.disabled = false;
      
      // Restore original text
      if (button.dataset.originalText) {
        button.innerHTML = button.dataset.originalText;
      }
    }
  }
  
  /**
   * Show a loading indicator for a section
   * @param {HTMLElement} section - The section element
   * @param {boolean} isLoading - Whether the section is loading
   */
  function setSectionLoading(section, isLoading) {
    if (isLoading) {
      // Add loading class
      section.classList.add('is-loading');
      
      // Create loading overlay if it doesn't exist
      let loadingOverlay = section.querySelector('.section-loading-overlay');
      if (!loadingOverlay) {
        loadingOverlay = document.createElement('div');
        loadingOverlay.className = 'section-loading-overlay';
        loadingOverlay.innerHTML = `
          <div class="spinner-border" role="status">
            <span class="visually-hidden">Loading...</span>
          </div>
        `;
        section.appendChild(loadingOverlay);
      }
      
      // Show loading overlay
      loadingOverlay.style.display = 'flex';
      
      // Announce to screen readers
      announceToScreenReader('Loading content, please wait...');
    } else {
      // Remove loading class
      section.classList.remove('is-loading');
      
      // Hide loading overlay
      const loadingOverlay = section.querySelector('.section-loading-overlay');
      if (loadingOverlay) {
        loadingOverlay.style.display = 'none';
      }
    }
  }
  
  // Public API
  return {
    handleApiError,
    handleFormError,
    setFormLoading,
    setButtonLoading,
    setSectionLoading
  };
})();
/**
 * CryoProtect Analyzer - Accessible Form Validation
 * This file contains functions to improve form accessibility and validation
 */

document.addEventListener('DOMContentLoaded', function() {
  // Initialize accessible form validation for all forms
  initAccessibleForms();
});

/**
 * Initialize accessible form validation
 */
function initAccessibleForms() {
  // Get all forms
  const forms = document.querySelectorAll('form');
  
  forms.forEach(form => {
    // Add accessible validation to the form
    addAccessibleValidation(form);
    
    // Add accessible error handling
    addAccessibleErrorHandling(form);
    
    // Add validation requirement indicators
    addValidationRequirements(form);
    
    // Add real-time validation
    addRealTimeValidation(form);
  });
  
  // Initialize password strength meters
  initPasswordStrengthMeters();
}

/**
 * Add accessible validation to a form
 * @param {HTMLFormElement} form - The form element
 */
function addAccessibleValidation(form) {
  // Add validation on form submission
  form.addEventListener('submit', function(event) {
    // Clear previous errors
    clearFormErrors(form);
    
    // Check for validity
    const isValid = validateForm(form);
    
    // If not valid, prevent submission
    if (!isValid) {
      event.preventDefault();
      
      // Focus the first invalid field
      const firstInvalidField = form.querySelector(':invalid');
      if (firstInvalidField) {
        firstInvalidField.focus();
      }
      
      // Announce validation errors to screen readers
      announceValidationErrors(form);
    }
  });
  
  // Add validation on field blur
  const formFields = form.querySelectorAll('input, select, textarea');
  formFields.forEach(field => {
    field.addEventListener('blur', function() {
      validateField(field);
    });
    
    // Add aria-describedby for fields with help text
    const fieldId = field.id;
    if (fieldId) {
      const helpText = form.querySelector(`.form-text[for="${fieldId}"]`) || 
                      form.querySelector(`.form-text#${fieldId}-help`);
      
      if (helpText) {
        const helpTextId = helpText.id || `${fieldId}-help`;
        helpText.id = helpTextId;
        
        // Set aria-describedby to connect field with help text
        const currentDescribedBy = field.getAttribute('aria-describedby') || '';
        if (!currentDescribedBy.includes(helpTextId)) {
          field.setAttribute('aria-describedby', 
            currentDescribedBy ? `${currentDescribedBy} ${helpTextId}` : helpTextId);
        }
      }
    }
  });
}

/**
 * Add accessible error handling to a form
 * @param {HTMLFormElement} form - The form element
 */
function addAccessibleErrorHandling(form) {
  // Create or get error container for screen reader announcements
  let srErrorContainer = form.querySelector('.sr-validation-errors');
  if (!srErrorContainer) {
    srErrorContainer = document.createElement('div');
    srErrorContainer.className = 'sr-validation-errors sr-only';
    srErrorContainer.setAttribute('aria-live', 'assertive');
    form.appendChild(srErrorContainer);
  }
  
  // Create error containers for each field if they don't exist
  const formFields = form.querySelectorAll('input, select, textarea');
  formFields.forEach(field => {
    const fieldId = field.id;
    if (!fieldId) return;
    
    // Check if error container already exists
    let errorContainer = form.querySelector(`#${fieldId}-error`);
    if (!errorContainer) {
      errorContainer = document.createElement('div');
      errorContainer.id = `${fieldId}-error`;
      errorContainer.className = 'invalid-feedback';
      
      // Insert error container after the field
      field.parentNode.insertBefore(errorContainer, field.nextSibling);
    }
    
    // Connect field to error container with aria-describedby
    const currentDescribedBy = field.getAttribute('aria-describedby') || '';
    if (!currentDescribedBy.includes(errorContainer.id)) {
      field.setAttribute('aria-describedby', 
        currentDescribedBy ? `${currentDescribedBy} ${errorContainer.id}` : errorContainer.id);
    }
  });
}

/**
 * Validate a form
 * @param {HTMLFormElement} form - The form element
 * @returns {boolean} - Whether the form is valid
 */
function validateForm(form) {
  let isValid = true;
  
  // Validate each field
  const formFields = form.querySelectorAll('input, select, textarea');
  formFields.forEach(field => {
    if (!validateField(field)) {
      isValid = false;
    }
  });
  
  return isValid;
}

/**
 * Validate a field
 * @param {HTMLElement} field - The field element
 * @returns {boolean} - Whether the field is valid
 */
function validateField(field) {
  // Skip disabled or hidden fields
  if (field.disabled || field.type === 'hidden') {
    return true;
  }
  
  // Check validity
  const isValid = field.validity.valid;
  
  // Update field styling and ARIA attributes
  if (isValid) {
    field.classList.remove('is-invalid');
    field.classList.add('is-valid');
    field.setAttribute('aria-invalid', 'false');
  } else {
    field.classList.remove('is-valid');
    field.classList.add('is-invalid');
    field.setAttribute('aria-invalid', 'true');
    
    // Update error message
    const errorContainer = document.getElementById(`${field.id}-error`);
    if (errorContainer) {
      errorContainer.textContent = getValidationMessage(field);
    }
  }
  
  return isValid;
}

/**
 * Get validation message for a field
 * @param {HTMLElement} field - The field element
 * @returns {string} - The validation message
 */
function getValidationMessage(field) {
  // Return custom validation message if set
  if (field.validationMessage) {
    return field.validationMessage;
  }
  
  // Generate message based on validity state
  if (field.validity.valueMissing) {
    return `${getFieldLabel(field)} is required.`;
  } else if (field.validity.typeMismatch) {
    return `Please enter a valid ${field.type} for ${getFieldLabel(field)}.`;
  } else if (field.validity.tooShort) {
    return `${getFieldLabel(field)} must be at least ${field.minLength} characters.`;
  } else if (field.validity.tooLong) {
    return `${getFieldLabel(field)} must be no more than ${field.maxLength} characters.`;
  } else if (field.validity.rangeUnderflow) {
    return `${getFieldLabel(field)} must be at least ${field.min}.`;
  } else if (field.validity.rangeOverflow) {
    return `${getFieldLabel(field)} must be no more than ${field.max}.`;
  } else if (field.validity.patternMismatch) {
    return `${getFieldLabel(field)} does not match the required pattern.`;
  }
  
  return 'This field is invalid.';
}

/**
 * Get label text for a field
 * @param {HTMLElement} field - The field element
 * @returns {string} - The label text
 */
function getFieldLabel(field) {
  // Try to find associated label
  const fieldId = field.id;
  if (fieldId) {
    const label = document.querySelector(`label[for="${fieldId}"]`);
    if (label) {
      return label.textContent.trim();
    }
  }
  
  // Use placeholder as fallback
  if (field.placeholder) {
    return field.placeholder;
  }
  
  // Use name as last resort
  if (field.name) {
    // Convert camelCase or snake_case to words
    return field.name
      .replace(/_/g, ' ')
      .replace(/([A-Z])/g, ' $1')
      .replace(/^./, str => str.toUpperCase())
      .trim();
  }
  
  return 'This field';
}

/**
 * Clear all form errors
 * @param {HTMLFormElement} form - The form element
 */
function clearFormErrors(form) {
  // Clear field errors
  const formFields = form.querySelectorAll('input, select, textarea');
  formFields.forEach(field => {
    field.classList.remove('is-invalid');
    field.setAttribute('aria-invalid', 'false');
    
    const errorContainer = document.getElementById(`${field.id}-error`);
    if (errorContainer) {
      errorContainer.textContent = '';
    }
  });
  
  // Clear screen reader error container
  const srErrorContainer = form.querySelector('.sr-validation-errors');
  if (srErrorContainer) {
    srErrorContainer.textContent = '';
  }
}

/**
 * Add validation requirement indicators to form fields
 * @param {HTMLFormElement} form - The form element
 */
function addValidationRequirements(form) {
  const formFields = form.querySelectorAll('input, select, textarea');
  
  formFields.forEach(field => {
    // Skip fields that don't need requirement indicators
    if (field.type === 'hidden' || field.type === 'submit' || field.type === 'button' ||
        field.type === 'reset' || field.type === 'checkbox' || field.type === 'radio') {
      return;
    }
    
    // Get field container (usually a div with class form-group)
    const fieldContainer = field.closest('.form-group') || field.parentNode;
    
    // Create requirements container if it doesn't exist
    let requirementsContainer = fieldContainer.querySelector('.field-requirements');
    if (!requirementsContainer) {
      requirementsContainer = document.createElement('div');
      requirementsContainer.className = 'field-requirements small text-muted mt-1';
      requirementsContainer.id = `${field.id}-requirements`;
      requirementsContainer.setAttribute('aria-live', 'polite');
      
      // Insert after the field or after the error container
      const errorContainer = document.getElementById(`${field.id}-error`);
      if (errorContainer) {
        errorContainer.parentNode.insertBefore(requirementsContainer, errorContainer.nextSibling);
      } else {
        field.parentNode.insertBefore(requirementsContainer, field.nextSibling);
      }
      
      // Connect field to requirements container with aria-describedby
      const currentDescribedBy = field.getAttribute('aria-describedby') || '';
      if (!currentDescribedBy.includes(requirementsContainer.id)) {
        field.setAttribute('aria-describedby',
          currentDescribedBy ? `${currentDescribedBy} ${requirementsContainer.id}` : requirementsContainer.id);
      }
    }
    
    // Add requirement indicators based on field attributes
    const requirements = [];
    
    // Required field
    if (field.required) {
      requirements.push('<span class="requirement" data-requirement="required"><i class="bi bi-asterisk text-danger"></i> Required</span>');
    }
    
    // Min length
    if (field.minLength && field.minLength > 0) {
      requirements.push(`<span class="requirement" data-requirement="minlength"><i class="bi bi-rulers"></i> At least ${field.minLength} characters</span>`);
    }
    
    // Max length
    if (field.maxLength && field.maxLength > 0 && field.maxLength < 524288) { // Ignore default max length
      requirements.push(`<span class="requirement" data-requirement="maxlength"><i class="bi bi-rulers"></i> Maximum ${field.maxLength} characters</span>`);
    }
    
    // Min/max for number inputs
    if (field.type === 'number') {
      if (field.min) {
        requirements.push(`<span class="requirement" data-requirement="min"><i class="bi bi-chevron-double-right"></i> Minimum value: ${field.min}</span>`);
      }
      if (field.max) {
        requirements.push(`<span class="requirement" data-requirement="max"><i class="bi bi-chevron-double-left"></i> Maximum value: ${field.max}</span>`);
      }
    }
    
    // Pattern
    if (field.pattern) {
      const patternDescription = field.getAttribute('data-pattern-description') || 'Must match the required pattern';
      requirements.push(`<span class="requirement" data-requirement="pattern"><i class="bi bi-braces"></i> ${patternDescription}</span>`);
    }
    
    // Email
    if (field.type === 'email') {
      requirements.push('<span class="requirement" data-requirement="email"><i class="bi bi-envelope"></i> Valid email address required</span>');
    }
    
    // URL
    if (field.type === 'url') {
      requirements.push('<span class="requirement" data-requirement="url"><i class="bi bi-link"></i> Valid URL required</span>');
    }
    
    // Password
    if (field.type === 'password') {
      // Add password strength requirements
      const passwordRequirements = field.getAttribute('data-password-requirements');
      if (passwordRequirements) {
        const requirements = passwordRequirements.split(',');
        
        requirements.forEach(req => {
          const [type, value] = req.trim().split(':');
          
          switch (type) {
            case 'minlength':
              field.setAttribute('minlength', value);
              break;
            case 'uppercase':
              field.setAttribute('data-require-uppercase', value || 'true');
              requirements.push('<span class="requirement" data-requirement="uppercase"><i class="bi bi-type-uppercase"></i> At least one uppercase letter</span>');
              break;
            case 'lowercase':
              field.setAttribute('data-require-lowercase', value || 'true');
              requirements.push('<span class="requirement" data-requirement="lowercase"><i class="bi bi-type-lowercase"></i> At least one lowercase letter</span>');
              break;
            case 'number':
              field.setAttribute('data-require-number', value || 'true');
              requirements.push('<span class="requirement" data-requirement="number"><i class="bi bi-123"></i> At least one number</span>');
              break;
            case 'special':
              field.setAttribute('data-require-special', value || 'true');
              requirements.push('<span class="requirement" data-requirement="special"><i class="bi bi-asterisk"></i> At least one special character</span>');
              break;
          }
        });
      }
      
      // Add password strength meter
      const strengthMeter = document.createElement('div');
      strengthMeter.className = 'password-strength-meter mt-2';
      strengthMeter.id = `${field.id}-strength-meter`;
      strengthMeter.innerHTML = `
        <div class="progress" role="progressbar" aria-label="Password strength" aria-valuenow="0" aria-valuemin="0" aria-valuemax="100">
          <div class="progress-bar" style="width: 0%"></div>
        </div>
        <div class="password-strength-text small mt-1">Password strength: <span>Not entered</span></div>
      `;
      
      requirementsContainer.parentNode.insertBefore(strengthMeter, requirementsContainer.nextSibling);
    }
    
    // Set requirements HTML
    if (requirements.length > 0) {
      requirementsContainer.innerHTML = requirements.join('<br>');
    } else {
      requirementsContainer.style.display = 'none';
    }
  });
}

/**
 * Add real-time validation to form fields
 * @param {HTMLFormElement} form - The form element
 */
function addRealTimeValidation(form) {
  const formFields = form.querySelectorAll('input, select, textarea');
  
  formFields.forEach(field => {
    // Skip fields that don't need real-time validation
    if (field.type === 'hidden' || field.type === 'submit' || field.type === 'button' || field.type === 'reset') {
      return;
    }
    
    // Add input event listener for real-time validation
    field.addEventListener('input', function() {
      // Don't validate on every keystroke for performance reasons
      clearTimeout(field.validationTimeout);
      field.validationTimeout = setTimeout(() => {
        validateFieldRealTime(field);
      }, 300);
    });
    
    // Add change event listener for select fields
    if (field.tagName === 'SELECT') {
      field.addEventListener('change', function() {
        validateFieldRealTime(field);
      });
    }
    
    // Special handling for password fields
    if (field.type === 'password') {
      field.addEventListener('input', function() {
        updatePasswordStrength(field);
      });
    }
  });
  
  // Add special handling for password confirmation fields
  const passwordConfirmFields = form.querySelectorAll('input[data-confirm-password]');
  passwordConfirmFields.forEach(confirmField => {
    const passwordFieldId = confirmField.getAttribute('data-confirm-password');
    const passwordField = document.getElementById(passwordFieldId);
    
    if (passwordField) {
      // Validate confirmation when password changes
      passwordField.addEventListener('input', function() {
        if (confirmField.value) {
          validatePasswordConfirmation(confirmField, passwordField);
        }
      });
      
      // Validate confirmation when confirmation field changes
      confirmField.addEventListener('input', function() {
        validatePasswordConfirmation(confirmField, passwordField);
      });
    }
  });
}

/**
 * Validate a field in real-time
 * @param {HTMLElement} field - The field element
 */
function validateFieldRealTime(field) {
  // Skip validation if field is empty and not required
  if (!field.value && !field.required) {
    // Reset validation state
    field.classList.remove('is-invalid', 'is-valid');
    field.removeAttribute('aria-invalid');
    
    // Clear error message
    const errorContainer = document.getElementById(`${field.id}-error`);
    if (errorContainer) {
      errorContainer.textContent = '';
    }
    
    return true;
  }
  
  // Perform validation
  const isValid = field.validity.valid;
  
  // Update field styling and ARIA attributes
  if (isValid) {
    field.classList.remove('is-invalid');
    field.classList.add('is-valid');
    field.setAttribute('aria-invalid', 'false');
    
    // Clear error message
    const errorContainer = document.getElementById(`${field.id}-error`);
    if (errorContainer) {
      errorContainer.textContent = '';
    }
    
    // Update requirements
    updateRequirements(field, true);
  } else {
    field.classList.remove('is-valid');
    field.classList.add('is-invalid');
    field.setAttribute('aria-invalid', 'true');
    
    // Update error message
    const errorContainer = document.getElementById(`${field.id}-error`);
    if (errorContainer) {
      errorContainer.textContent = getValidationMessage(field);
    }
    
    // Update requirements
    updateRequirements(field, false);
  }
  
  return isValid;
}

/**
 * Update validation requirements visual state
 * @param {HTMLElement} field - The field element
 * @param {boolean} isValid - Whether the field is valid
 */
function updateRequirements(field, isValid) {
  const requirementsContainer = document.getElementById(`${field.id}-requirements`);
  if (!requirementsContainer) return;
  
  const requirements = requirementsContainer.querySelectorAll('.requirement');
  
  requirements.forEach(requirement => {
    const requirementType = requirement.getAttribute('data-requirement');
    let requirementMet = true;
    
    // Check if this specific requirement is met
    switch (requirementType) {
      case 'required':
        requirementMet = field.value.trim() !== '';
        break;
      case 'minlength':
        requirementMet = field.value.length >= field.minLength;
        break;
      case 'maxlength':
        requirementMet = field.value.length <= field.maxLength;
        break;
      case 'min':
        requirementMet = parseFloat(field.value) >= parseFloat(field.min);
        break;
      case 'max':
        requirementMet = parseFloat(field.value) <= parseFloat(field.max);
        break;
      case 'pattern':
        requirementMet = new RegExp(field.pattern).test(field.value);
        break;
      case 'email':
        requirementMet = /^[^\s@]+@[^\s@]+\.[^\s@]+$/.test(field.value);
        break;
      case 'url':
        requirementMet = /^(https?:\/\/)?([\da-z\.-]+)\.([a-z\.]{2,6})([\/\w \.-]*)*\/?$/.test(field.value);
        break;
      case 'uppercase':
        requirementMet = /[A-Z]/.test(field.value);
        break;
      case 'lowercase':
        requirementMet = /[a-z]/.test(field.value);
        break;
      case 'number':
        requirementMet = /[0-9]/.test(field.value);
        break;
      case 'special':
        requirementMet = /[^A-Za-z0-9]/.test(field.value);
        break;
    }
    
    // Update requirement styling
    if (requirementMet) {
      requirement.classList.remove('text-danger');
      requirement.classList.add('text-success');
      const icon = requirement.querySelector('i');
      if (icon) {
        icon.className = icon.className.replace('bi-asterisk', 'bi-check-circle');
        icon.className = icon.className.replace('text-danger', 'text-success');
      }
    } else {
      requirement.classList.remove('text-success');
      requirement.classList.add('text-danger');
      const icon = requirement.querySelector('i');
      if (icon) {
        icon.className = icon.className.replace('bi-check-circle', 'bi-asterisk');
        icon.className = icon.className.replace('text-success', 'text-danger');
      }
    }
  });
}

/**
 * Validate password confirmation field
 * @param {HTMLElement} confirmField - The confirmation field
 * @param {HTMLElement} passwordField - The password field
 */
function validatePasswordConfirmation(confirmField, passwordField) {
  const isValid = confirmField.value === passwordField.value;
  
  if (isValid) {
    confirmField.classList.remove('is-invalid');
    confirmField.classList.add('is-valid');
    confirmField.setAttribute('aria-invalid', 'false');
    
    // Clear error message
    const errorContainer = document.getElementById(`${confirmField.id}-error`);
    if (errorContainer) {
      errorContainer.textContent = '';
    }
  } else {
    confirmField.classList.remove('is-valid');
    confirmField.classList.add('is-invalid');
    confirmField.setAttribute('aria-invalid', 'true');
    
    // Update error message
    const errorContainer = document.getElementById(`${confirmField.id}-error`);
    if (errorContainer) {
      errorContainer.textContent = 'Passwords do not match.';
    }
  }
  
  return isValid;
}

/**
 * Initialize password strength meters
 */
function initPasswordStrengthMeters() {
  const passwordFields = document.querySelectorAll('input[type="password"]:not([data-confirm-password])');
  
  passwordFields.forEach(field => {
    // Initial update
    updatePasswordStrength(field);
  });
}

/**
 * Update password strength meter
 * @param {HTMLElement} field - The password field
 */
function updatePasswordStrength(field) {
  const strengthMeter = document.getElementById(`${field.id}-strength-meter`);
  if (!strengthMeter) return;
  
  const progressBar = strengthMeter.querySelector('.progress-bar');
  const strengthText = strengthMeter.querySelector('.password-strength-text span');
  
  if (!progressBar || !strengthText) return;
  
  // Calculate password strength
  const strength = calculatePasswordStrength(field.value);
  
  // Update progress bar
  progressBar.style.width = `${strength.score}%`;
  progressBar.setAttribute('aria-valuenow', strength.score);
  
  // Update strength class
  progressBar.className = 'progress-bar';
  if (strength.score < 25) {
    progressBar.classList.add('bg-danger');
  } else if (strength.score < 50) {
    progressBar.classList.add('bg-warning');
  } else if (strength.score < 75) {
    progressBar.classList.add('bg-info');
  } else {
    progressBar.classList.add('bg-success');
  }
  
  // Update strength text
  strengthText.textContent = strength.label;
  
  // Announce to screen readers
  const srAnnouncements = document.getElementById('sr-announcements');
  if (srAnnouncements && field.value) {
    srAnnouncements.textContent = `Password strength: ${strength.label}`;
  }
}

/**
 * Calculate password strength
 * @param {string} password - The password to evaluate
 * @returns {Object} Strength score and label
 */
function calculatePasswordStrength(password) {
  if (!password) {
    return { score: 0, label: 'Not entered' };
  }
  
  let score = 0;
  
  // Length
  if (password.length >= 8) score += 25;
  if (password.length >= 12) score += 15;
  if (password.length >= 16) score += 10;
  
  // Complexity
  if (/[A-Z]/.test(password)) score += 10;
  if (/[a-z]/.test(password)) score += 10;
  if (/[0-9]/.test(password)) score += 10;
  if (/[^A-Za-z0-9]/.test(password)) score += 15;
  
  // Variety
  const uniqueChars = new Set(password).size;
  score += Math.min(uniqueChars * 2, 15);
  
  // Cap at 100
  score = Math.min(score, 100);
  
  // Determine label
  let label;
  if (score < 25) {
    label = 'Very weak';
  } else if (score < 50) {
    label = 'Weak';
  } else if (score < 75) {
    label = 'Good';
  } else if (score < 90) {
    label = 'Strong';
  } else {
    label = 'Very strong';
  }
  
  return { score, label };
}

/**
 * Announce validation errors to screen readers
 * @param {HTMLFormElement} form - The form element
 */
function announceValidationErrors(form) {
  // Get all error messages
  const errorMessages = [];
  const invalidFields = form.querySelectorAll('.is-invalid');
  
  invalidFields.forEach(field => {
    const errorContainer = document.getElementById(`${field.id}-error`);
    if (errorContainer && errorContainer.textContent) {
      errorMessages.push(errorContainer.textContent);
    }
  });
  
  // Announce errors to screen readers
  if (errorMessages.length > 0) {
    const srErrorContainer = form.querySelector('.sr-validation-errors');
    if (srErrorContainer) {
      srErrorContainer.textContent = `Form has ${errorMessages.length} error${errorMessages.length > 1 ? 's' : ''}: ${errorMessages.join('. ')}`;
    }
  }
}
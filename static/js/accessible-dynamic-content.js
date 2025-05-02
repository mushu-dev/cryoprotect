/**
 * CryoProtect Analyzer - Accessible Dynamic Content
 * This file contains functions to improve accessibility of dynamically loaded content
 */

document.addEventListener('DOMContentLoaded', function() {
  // Initialize accessible dynamic content handling
  initAccessibleDynamicContent();
});

/**
 * Initialize accessible dynamic content handling
 */
function initAccessibleDynamicContent() {
  // Set up observers for dynamic content containers
  setupDynamicContentObservers();
  
  // Enhance existing dynamic content containers
  enhanceDynamicContentContainers();
  
  // Add accessibility enhancements to charts
  enhanceChartAccessibility();
  
  // Add accessibility enhancements to molecular viewers
  enhanceMolecularViewerAccessibility();
}

/**
 * Set up mutation observers for dynamic content containers
 */
function setupDynamicContentObservers() {
  // Get all containers that will receive dynamic content
  const dynamicContainers = document.querySelectorAll('[id$="-container"]');
  
  // Create a mutation observer
  const observer = new MutationObserver(function(mutations) {
    mutations.forEach(function(mutation) {
      if (mutation.type === 'childList' && mutation.addedNodes.length > 0) {
        // Process added nodes for accessibility
        mutation.addedNodes.forEach(function(node) {
          if (node.nodeType === Node.ELEMENT_NODE) {
            enhanceDynamicElement(node);
          }
        });
      }
    });
  });
  
  // Observe each container
  dynamicContainers.forEach(function(container) {
    observer.observe(container, { childList: true, subtree: true });
    
    // Ensure container has aria-live attribute if it doesn't already
    if (!container.hasAttribute('aria-live')) {
      container.setAttribute('aria-live', 'polite');
    }
  });
}

/**
 * Enhance existing dynamic content containers
 */
function enhanceDynamicContentContainers() {
  // Get all containers that will receive dynamic content
  const dynamicContainers = document.querySelectorAll('[id$="-container"]');
  
  dynamicContainers.forEach(function(container) {
    // Ensure container has aria-live attribute
    if (!container.hasAttribute('aria-live')) {
      container.setAttribute('aria-live', 'polite');
    }
    
    // Process existing content
    Array.from(container.children).forEach(function(child) {
      if (child.nodeType === Node.ELEMENT_NODE) {
        enhanceDynamicElement(child);
      }
    });
  });
}

/**
 * Enhance a dynamically added element for accessibility
 * @param {HTMLElement} element - The element to enhance
 */
function enhanceDynamicElement(element) {
  // Add appropriate ARIA attributes based on element type
  if (element.classList.contains('card')) {
    enhanceCard(element);
  } else if (element.classList.contains('alert')) {
    enhanceAlert(element);
  } else if (element.tagName === 'TABLE') {
    enhanceTable(element);
  } else if (element.classList.contains('spinner-border')) {
    enhanceSpinner(element);
  }
  
  // Process child elements
  Array.from(element.children).forEach(function(child) {
    if (child.nodeType === Node.ELEMENT_NODE) {
      enhanceDynamicElement(child);
    }
  });
}

/**
 * Enhance a card element for accessibility
 * @param {HTMLElement} card - The card element
 */
function enhanceCard(card) {
  // Ensure card header has appropriate role
  const cardHeader = card.querySelector('.card-header');
  if (cardHeader) {
    // If card header contains a heading, make sure it's properly associated
    const heading = cardHeader.querySelector('h1, h2, h3, h4, h5, h6');
    if (heading) {
      const headingId = heading.id || `heading-${Date.now()}-${Math.floor(Math.random() * 1000)}`;
      heading.id = headingId;
      card.setAttribute('aria-labelledby', headingId);
    }
  }
  
  // Ensure interactive elements within card are accessible
  const interactiveElements = card.querySelectorAll('button, a, input, select, textarea');
  interactiveElements.forEach(function(element) {
    if (!element.hasAttribute('aria-label') && !element.hasAttribute('aria-labelledby')) {
      // Try to find a label from nearby text
      const buttonText = element.textContent.trim();
      if (buttonText) {
        element.setAttribute('aria-label', buttonText);
      }
    }
  });
}

/**
 * Enhance an alert element for accessibility
 * @param {HTMLElement} alert - The alert element
 */
function enhanceAlert(alert) {
  // Ensure alert has appropriate role and live region attributes
  if (!alert.hasAttribute('role')) {
    alert.setAttribute('role', 'alert');
  }
  
  if (!alert.hasAttribute('aria-live')) {
    alert.setAttribute('aria-live', 'assertive');
  }
}

/**
 * Enhance a table element for accessibility
 * @param {HTMLElement} table - The table element
 */
function enhanceTable(table) {
  // Ensure table has appropriate ARIA attributes
  if (!table.hasAttribute('role')) {
    table.setAttribute('role', 'table');
  }
  
  // Add caption if missing
  if (!table.querySelector('caption')) {
    const tableId = table.id;
    if (tableId) {
      // Try to find a heading that might describe this table
      const heading = document.querySelector(`h1[id$="${tableId}"], h2[id$="${tableId}"], h3[id$="${tableId}"], h4[id$="${tableId}"], h5[id$="${tableId}"], h6[id$="${tableId}"]`);
      if (heading) {
        const caption = document.createElement('caption');
        caption.textContent = heading.textContent;
        table.prepend(caption);
      }
    }
  }
  
  // Ensure thead has appropriate role
  const thead = table.querySelector('thead');
  if (thead && !thead.hasAttribute('role')) {
    thead.setAttribute('role', 'rowgroup');
  }
  
  // Ensure tbody has appropriate role
  const tbody = table.querySelector('tbody');
  if (tbody && !tbody.hasAttribute('role')) {
    tbody.setAttribute('role', 'rowgroup');
  }
  
  // Ensure rows have appropriate role
  const rows = table.querySelectorAll('tr');
  rows.forEach(function(row) {
    if (!row.hasAttribute('role')) {
      row.setAttribute('role', 'row');
    }
  });
  
  // Ensure header cells have appropriate role
  const headerCells = table.querySelectorAll('th');
  headerCells.forEach(function(cell) {
    if (!cell.hasAttribute('role')) {
      cell.setAttribute('role', 'columnheader');
    }
    
    // Add scope attribute if missing
    if (!cell.hasAttribute('scope')) {
      cell.setAttribute('scope', 'col');
    }
  });
  
  // Ensure data cells have appropriate role
  const dataCells = table.querySelectorAll('td');
  dataCells.forEach(function(cell) {
    if (!cell.hasAttribute('role')) {
      cell.setAttribute('role', 'cell');
    }
  });
}

/**
 * Enhance a spinner element for accessibility
 * @param {HTMLElement} spinner - The spinner element
 */
function enhanceSpinner(spinner) {
  // Ensure spinner has appropriate role
  if (!spinner.hasAttribute('role')) {
    spinner.setAttribute('role', 'status');
  }
  
  // Ensure spinner has a visually hidden text
  let hiddenText = spinner.querySelector('.visually-hidden');
  if (!hiddenText) {
    hiddenText = document.createElement('span');
    hiddenText.className = 'visually-hidden';
    hiddenText.textContent = 'Loading...';
    spinner.appendChild(hiddenText);
  }
  
  // Ensure parent container has aria-live
  const container = spinner.parentElement;
  if (container && !container.hasAttribute('aria-live')) {
    container.setAttribute('aria-live', 'polite');
  }
}

/**
 * Enhance chart accessibility
 */
function enhanceChartAccessibility() {
  // Get all chart containers
  const chartContainers = document.querySelectorAll('.chart-container');
  
  chartContainers.forEach(function(container) {
    // Ensure chart container has appropriate role
    if (!container.hasAttribute('role')) {
      container.setAttribute('role', 'img');
    }
    
    // Ensure chart has an accessible name
    if (!container.hasAttribute('aria-label') && !container.hasAttribute('aria-labelledby')) {
      // Try to find a heading that might describe this chart
      const heading = container.previousElementSibling;
      if (heading && /^h[1-6]$/i.test(heading.tagName)) {
        const headingId = heading.id || `heading-${Date.now()}-${Math.floor(Math.random() * 1000)}`;
        heading.id = headingId;
        container.setAttribute('aria-labelledby', headingId);
      } else {
        // Default label
        container.setAttribute('aria-label', 'Data visualization chart');
      }
    }
    
    // Check for description
    const descId = container.getAttribute('aria-describedby');
    if (!descId) {
      // Look for a description element with the same ID + "-description"
      const containerId = container.id;
      if (containerId) {
        const descElement = document.getElementById(`${containerId}-description`);
        if (descElement) {
          container.setAttribute('aria-describedby', descElement.id);
        }
      }
    }
  });
}

/**
 * Enhance molecular viewer accessibility
 */
function enhanceMolecularViewerAccessibility() {
  // Get all molecular viewer containers
  const viewerContainers = document.querySelectorAll('#molecule-viewer, .molecule-viewer');
  
  viewerContainers.forEach(function(container) {
    // Ensure viewer container has appropriate role
    if (!container.hasAttribute('role')) {
      container.setAttribute('role', 'img');
    }
    
    // Ensure viewer has an accessible name
    if (!container.hasAttribute('aria-label') && !container.hasAttribute('aria-labelledby')) {
      container.setAttribute('aria-label', 'Molecular structure 3D visualization');
    }
    
    // Add keyboard support for viewer controls
    const viewerControls = document.querySelectorAll('#molecule-viewer-toolbar button, .molecule-viewer-toolbar button');
    viewerControls.forEach(function(control) {
      if (!control.hasAttribute('aria-label')) {
        const buttonText = control.textContent.trim();
        if (buttonText) {
          control.setAttribute('aria-label', buttonText);
        } else {
          // Try to determine function from icon
          const icon = control.querySelector('i');
          if (icon) {
            const iconClass = Array.from(icon.classList).find(cls => cls.startsWith('bi-'));
            if (iconClass) {
              const controlName = iconClass.replace('bi-', '').replace(/-/g, ' ');
              control.setAttribute('aria-label', controlName);
            }
          }
        }
      }
    });
  });
}
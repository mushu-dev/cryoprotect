/**
 * CryoProtect Analyzer Web Interface
 * Authentication UI Module - Handles UI updates based on authentication state
 */

document.addEventListener('DOMContentLoaded', function() {
  // Initialize authentication state
  initAuthUI();
  
  // Add event listener for logout button
  const logoutButton = document.getElementById('logout-button');
  if (logoutButton) {
    logoutButton.addEventListener('click', function(e) {
      e.preventDefault();
      handleLogout();
    });
  }
});

/**
 * Initialize the authentication UI
 */
function initAuthUI() {
  // Check if user is logged in
  Auth.checkSession().then(user => {
    updateAuthUI(user);
  }).catch(error => {
    console.error('Error checking authentication:', error);
    updateAuthUI(null);
  });
}

/**
 * Update the UI based on authentication state
 * @param {Object|null} user - The authenticated user or null if not authenticated
 */
function updateAuthUI(user) {
  // Elements that should be shown when authenticated
  const authenticatedElements = document.querySelectorAll('.nav-authenticated');
  
  // Elements that should be shown when not authenticated
  const unauthenticatedElements = document.querySelectorAll('.nav-unauthenticated');
  
  // Elements that require authentication to be enabled
  const requiresAuthElements = document.querySelectorAll('.requires-auth');
  
  if (user) {
    // User is authenticated
    console.log('User authenticated:', user.email);
    
    // Update user name display
    const userNameElements = document.querySelectorAll('.user-name');
    userNameElements.forEach(element => {
      element.textContent = user.user_metadata?.full_name || user.email;
    });
    
    // Show authenticated elements
    authenticatedElements.forEach(element => {
      element.classList.remove('d-none');
    });
    
    // Hide unauthenticated elements
    unauthenticatedElements.forEach(element => {
      element.classList.add('d-none');
    });
    
    // Enable elements that require authentication
    requiresAuthElements.forEach(element => {
      element.classList.remove('disabled');
    });
    
    // Load user-specific data if needed
    loadUserData(user);
  } else {
    // User is not authenticated
    console.log('User not authenticated');
    
    // Hide authenticated elements
    authenticatedElements.forEach(element => {
      element.classList.add('d-none');
    });
    
    // Show unauthenticated elements
    unauthenticatedElements.forEach(element => {
      element.classList.remove('d-none');
    });
    
    // Disable elements that require authentication
    requiresAuthElements.forEach(element => {
      element.classList.add('disabled');
    });
  }
}

/**
 * Handle user logout
 */
function handleLogout() {
  Auth.signOut().then(() => {
    console.log('User logged out');
    // Redirect to login page
    window.location.href = '/login';
  }).catch(error => {
    console.error('Error logging out:', error);
    alert('Failed to log out. Please try again.');
  });
}

/**
 * Load user-specific data
 * @param {Object} user - The authenticated user
 */
function loadUserData(user) {
  // This function can be used to load user-specific data
  // For example, loading dashboard statistics, user preferences, etc.
  
  // Example: Load dashboard statistics
  loadDashboardStats(user);
}

/**
 * Load dashboard statistics for the authenticated user
 * @param {Object} user - The authenticated user
 */
function loadDashboardStats(user) {
  // This is a placeholder function that would typically make API calls
  // to fetch user-specific dashboard statistics
  
  // For demonstration purposes, we'll just update the counts with random values
  const moleculeCount = document.getElementById('molecule-count');
  const mixtureCount = document.getElementById('mixture-count');
  const predictionCount = document.getElementById('prediction-count');
  
  if (moleculeCount) {
    // In a real application, this would be an API call
    // For now, just set a random number for demonstration
    setTimeout(() => {
      moleculeCount.textContent = Math.floor(Math.random() * 100) + 50;
    }, 500);
  }
  
  if (mixtureCount) {
    setTimeout(() => {
      mixtureCount.textContent = Math.floor(Math.random() * 50) + 20;
    }, 700);
  }
  
  if (predictionCount) {
    setTimeout(() => {
      predictionCount.textContent = Math.floor(Math.random() * 30) + 10;
    }, 900);
  }
}
/**
 * CryoProtect Analyzer Web Interface
 * Authentication Module - Handles user authentication with Supabase
 */

const Auth = (function() {
  // Supabase client configuration
  let supabaseClient = null;
  
  // Token storage key
  const TOKEN_KEY = 'cryoprotect_auth_token';
  const USER_KEY = 'cryoprotect_user';
  
  /**
   * Initialize Supabase client
   * @returns {Object} Supabase client instance
   */
  function initClient() {
    if (!supabaseClient) {
      // Get Supabase URL and key from meta tags
      const supabaseUrl = document.querySelector('meta[name="supabase-url"]')?.content;
      const supabaseKey = document.querySelector('meta[name="supabase-key"]')?.content;
      
      if (!supabaseUrl || !supabaseKey) {
        console.error('Supabase configuration not found');
        return null;
      }
      
      // Initialize Supabase client
      supabaseClient = supabase.createClient(supabaseUrl, supabaseKey);
    }
    
    return supabaseClient;
  }
  
  /**
   * Sign in with email and password
   * @param {string} email - User email
   * @param {string} password - User password
   * @returns {Promise} Promise that resolves with the user data
   */
  async function signIn(email, password) {
    const client = initClient();
    if (!client) {
      throw new Error('Supabase client not initialized');
    }
    
    try {
      const { data, error } = await client.auth.signInWithPassword({
        email,
        password
      });
      
      if (error) {
        throw new Error(error.message);
      }
      
      // Store authentication data
      storeSession(data.session);
      storeUser(data.user);

      // Sync user profile with backend
      await syncUserProfile(data.user);

      return data.user;
    } catch (error) {
      console.error('Sign in error:', error);
      throw error;
    }
  }

  /**
   * Sign in with a social provider (Google, GitHub)
   * @param {string} provider - 'google' or 'github'
   */
  async function signInWithProvider(provider) {
    const client = initClient();
    if (!client) {
      throw new Error('Supabase client not initialized');
    }
    try {
      const { error } = await client.auth.signInWithOAuth({
        provider,
        options: {
          redirectTo: window.location.origin + '/login'
        }
      });
      if (error) {
        throw new Error(error.message);
      }
    } catch (error) {
      console.error('Social sign in error:', error);
      throw error;
    }
  }

  /**
   * Sync user profile with backend user_profile table
   * @param {Object} user - Supabase Auth user object
   */
  async function syncUserProfile(user) {
    if (!user) return;
    try {
      const token = getToken();
      await fetch('/api/v1/user_profile', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
          'Authorization': `Bearer ${token}`
        },
        body: JSON.stringify({
          email: user.email,
          name: user.user_metadata?.full_name || ''
        })
      });
    } catch (error) {
      console.error('Failed to sync user profile:', error);
    }
  }
  
  /**
   * Sign up a new user
   * @param {string} email - User email
   * @param {string} password - User password
   * @returns {Promise} Promise that resolves with the user data
   */
  async function signUp(email, password) {
    const client = initClient();
    if (!client) {
      throw new Error('Supabase client not initialized');
    }
    
    try {
      const { data, error } = await client.auth.signUp({
        email,
        password
      });
      
      if (error) {
        throw new Error(error.message);
      }
      
      // If email confirmation is required, data.user will be available
      // but data.session might be null
      if (data.session) {
        storeSession(data.session);
        storeUser(data.user);
        // Sync user profile with backend
        await syncUserProfile(data.user);
      }
      
      return data;
    } catch (error) {
      console.error('Sign up error:', error);
      throw error;
    }
  }
  
  /**
   * Sign out the current user
   * @returns {Promise} Promise that resolves when sign out is complete
   */
  async function signOut() {
    const client = initClient();
    if (!client) {
      throw new Error('Supabase client not initialized');
    }
    
    try {
      const { error } = await client.auth.signOut();
      
      if (error) {
        throw new Error(error.message);
      }
      
      // Clear stored authentication data
      clearSession();
      clearUser();
      
      return true;
    } catch (error) {
      console.error('Sign out error:', error);
      throw error;
    }
  }
  
  /**
   * Check if the user is signed in
   * @returns {Promise} Promise that resolves with the user data if signed in, or null if not
   */
  async function checkSession() {
    const client = initClient();
    if (!client) {
      return null;
    }
    
    try {
      // First check if we have a stored session
      const storedUser = getStoredUser();
      if (storedUser) {
        // Verify the session with Supabase
        const { data, error } = await client.auth.getUser();
        
        if (error) {
          clearSession();
          clearUser();
          return null;
        }
        
        return data.user;
      }
      
      return null;
    } catch (error) {
      console.error('Check session error:', error);
      clearSession();
      clearUser();
      return null;
    }
  }
  
  /**
   * Request a password reset for the given email
   * @param {string} email - The email address to send the reset link to
   * @returns {Promise} Promise that resolves when the reset email is sent
   */
  async function resetPassword(email) {
    const client = initClient();
    if (!client) {
      throw new Error('Supabase client not initialized');
    }
    
    try {
      const { error } = await client.auth.resetPasswordForEmail(email, {
        redirectTo: window.location.origin + '/reset-password'
      });
      
      if (error) {
        throw new Error(error.message);
      }
      
      return true;
    } catch (error) {
      console.error('Reset password error:', error);
      throw error;
    }
  }
  
  /**
   * Update the user's password
   * @param {string} newPassword - The new password
   * @returns {Promise} Promise that resolves when the password is updated
   */
  async function updatePassword(newPassword) {
    const client = initClient();
    if (!client) {
      throw new Error('Supabase client not initialized');
    }
    
    try {
      const { error } = await client.auth.updateUser({
        password: newPassword
      });
      
      if (error) {
        throw new Error(error.message);
      }
      
      return true;
    } catch (error) {
      console.error('Update password error:', error);
      throw error;
    }
  }
  
  /**
   * Update the user's profile data
   * @param {Object} userData - The user data to update
   * @returns {Promise} Promise that resolves with the updated user data
   */
  async function updateProfile(userData) {
    const client = initClient();
    if (!client) {
      throw new Error('Supabase client not initialized');
    }
    
    try {
      const { data, error } = await client.auth.updateUser({
        data: userData
      });
      
      if (error) {
        throw new Error(error.message);
      }
      
      // Update stored user data
      storeUser(data.user);
      
      return data.user;
    } catch (error) {
      console.error('Update profile error:', error);
      throw error;
    }
  }
  
  /**
   * Get the current authentication token
   * @returns {string|null} The authentication token or null if not signed in
   */
  function getToken() {
    return localStorage.getItem(TOKEN_KEY);
  }
  
  /**
   * Store the authentication session
   * @param {Object} session - The authentication session
   */
  function storeSession(session) {
    if (session && session.access_token) {
      localStorage.setItem(TOKEN_KEY, session.access_token);
    }
  }
  
  /**
   * Clear the stored authentication session
   */
  function clearSession() {
    localStorage.removeItem(TOKEN_KEY);
  }
  
  /**
   * Store the user data
   * @param {Object} user - The user data
   */
  function storeUser(user) {
    if (user) {
      localStorage.setItem(USER_KEY, JSON.stringify(user));
    }
  }
  
  /**
   * Get the stored user data
   * @returns {Object|null} The user data or null if not stored
   */
  function getStoredUser() {
    const userData = localStorage.getItem(USER_KEY);
    if (userData) {
      try {
        return JSON.parse(userData);
      } catch (error) {
        console.error('Error parsing stored user data:', error);
        clearUser();
      }
    }
    return null;
  }
  
  /**
   * Clear the stored user data
   */
  function clearUser() {
    localStorage.removeItem(USER_KEY);
  }
  
  // Public API
  return {
    signIn,
    signUp,
    signOut,
    checkSession,
    resetPassword,
    updatePassword,
    updateProfile,
    getToken,
    getStoredUser
  };
})();

// Initialize authentication on page load
document.addEventListener('DOMContentLoaded', function() {
  // Check if user is logged in
  Auth.checkSession().then(user => {
    if (user) {
      // Update UI for logged in user
      updateAuthUI(user);
    } else {
      // Update UI for logged out user
      updateAuthUI(null);
    }
  });
  
  // Set up login form
  const loginForm = document.getElementById('login-form');
  if (loginForm) {
    // Social login buttons
    const googleBtn = document.getElementById('google-login-btn');
    const githubBtn = document.getElementById('github-login-btn');
    if (googleBtn) {
      googleBtn.addEventListener('click', async function(e) {
        e.preventDefault();
        try {
          await Auth.signInWithProvider('google');
        } catch (error) {
          alert(error.message || 'Google login failed.');
        }
      });
    }
    if (githubBtn) {
      githubBtn.addEventListener('click', async function(e) {
        e.preventDefault();
        try {
          await Auth.signInWithProvider('github');
        } catch (error) {
          alert(error.message || 'GitHub login failed.');
        }
      });
    }
    loginForm.addEventListener('submit', async function(e) {
      e.preventDefault();

      const email = document.getElementById('email').value;
      const password = document.getElementById('password').value;

      // Remove any existing error messages
      const existingError = document.querySelector('.alert-danger');
      if (existingError) {
        existingError.remove();
      }
      // Remove any existing MFA prompt
      const existingMfa = document.getElementById('mfa-prompt');
      if (existingMfa) {
        existingMfa.remove();
      }

      try {
        // Use fetch to POST to /auth/login to handle MFA-required responses
        const response = await fetch('/auth/login', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ email, password })
        });

        if (response.status === 200) {
          // Login successful
          window.location.href = '/';
          return;
        } else if (response.status === 206) {
          // MFA required
          const data = await response.json();
          // Show MFA prompt
          const mfaDiv = document.createElement('div');
          mfaDiv.id = 'mfa-prompt';
          mfaDiv.className = 'mt-3';
          mfaDiv.innerHTML = `
            <label for="totp-code" class="form-label">Multi-Factor Authentication Code</label>
            <input type="text" id="totp-code" class="form-control mb-2" placeholder="Enter TOTP code" autocomplete="one-time-code" />
            <button id="mfa-submit-btn" class="btn btn-primary">Verify</button>
            <div id="mfa-error" class="alert alert-danger mt-2 d-none"></div>
          `;
          loginForm.appendChild(mfaDiv);

          document.getElementById('mfa-submit-btn').addEventListener('click', async function(ev) {
            ev.preventDefault();
            const totpCode = document.getElementById('totp-code').value;
            const mfaError = document.getElementById('mfa-error');
            mfaError.classList.add('d-none');
            try {
              const mfaResp = await fetch('/auth/mfa/verify', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({
                  mfa_token: data.mfa.factor_id || data.mfa.mfa_token,
                  totp_code: totpCode
                })
              });
              if (mfaResp.status === 200) {
                window.location.href = '/';
              } else {
                const mfaData = await mfaResp.json();
                mfaError.textContent = mfaData.message || 'MFA verification failed.';
                mfaError.classList.remove('d-none');
              }
            } catch (err) {
              mfaError.textContent = err.message || 'MFA verification failed.';
              mfaError.classList.remove('d-none');
            }
          });
          return;
        } else {
          // Other error
          const data = await response.json();
          throw new Error(data.message || 'Login failed. Please check your credentials.');
        }
      } catch (error) {
        const errorElement = document.createElement('div');
        errorElement.className = 'alert alert-danger mt-3';
        errorElement.textContent = error.message || 'Login failed. Please check your credentials.';
        loginForm.appendChild(errorElement);
      }
    });
  }
  
  // Set up registration form
  const registerForm = document.getElementById('register-form');
  if (registerForm) {
    registerForm.addEventListener('submit', async function(e) {
      e.preventDefault();
      
      const email = document.getElementById('email').value;
      const password = document.getElementById('password').value;
      const confirmPassword = document.getElementById('confirm-password').value;
      
      // Validate passwords match
      if (password !== confirmPassword) {
        const errorElement = document.getElementById('register-error');
        errorElement.textContent = 'Passwords do not match.';
        errorElement.style.display = 'block';
        return;
      }
      
      // Validate password length
      if (password.length < 8) {
        const errorElement = document.getElementById('register-error');
        errorElement.textContent = 'Password must be at least 8 characters long.';
        errorElement.style.display = 'block';
        return;
      }
      
      try {
        const { user, session } = await Auth.signUp(email, password);
        
        if (session) {
          // User is automatically signed in
          window.location.href = '/';
        } else {
          // Email confirmation required
          const errorElement = document.getElementById('register-error');
          errorElement.className = 'alert alert-info mt-3';
          errorElement.textContent = 'Please check your email to confirm your registration.';
          errorElement.style.display = 'block';
        }
      } catch (error) {
        const errorElement = document.getElementById('register-error');
        errorElement.className = 'alert alert-danger mt-3';
        errorElement.textContent = error.message || 'Registration failed. Please try again.';
        errorElement.style.display = 'block';
      }
    });
  }
  
  // Set up password reset request form
  const requestResetForm = document.getElementById('request-reset-form');
  if (requestResetForm) {
    requestResetForm.addEventListener('submit', async function(e) {
      e.preventDefault();
      
      const email = document.getElementById('reset-email').value;
      
      try {
        await Auth.resetPassword(email);
        
        // Show success message
        document.getElementById('reset-success').style.display = 'block';
        document.getElementById('reset-error').style.display = 'none';
        
        // Clear form
        requestResetForm.reset();
      } catch (error) {
        // Show error message
        const errorElement = document.getElementById('reset-error');
        errorElement.textContent = error.message || 'Failed to send reset email. Please try again.';
        errorElement.style.display = 'block';
        document.getElementById('reset-success').style.display = 'none';
      }
    });
  }
  
  // Set up update password form
  const updatePasswordForm = document.getElementById('update-password-form');
  if (updatePasswordForm) {
    // Check if we have a password reset token in the URL
    const hash = window.location.hash;
    if (hash && hash.includes('type=recovery')) {
      // Show the update password form
      document.getElementById('request-reset-form-container').style.display = 'none';
      document.getElementById('update-password-form-container').style.display = 'block';
      
      updatePasswordForm.addEventListener('submit', async function(e) {
        e.preventDefault();
        
        const newPassword = document.getElementById('new-password').value;
        const confirmPassword = document.getElementById('confirm-password').value;
        
        // Validate passwords match
        if (newPassword !== confirmPassword) {
          const errorElement = document.getElementById('reset-error');
          errorElement.textContent = 'Passwords do not match.';
          errorElement.style.display = 'block';
          return;
        }
        
        // Validate password length
        if (newPassword.length < 8) {
          const errorElement = document.getElementById('reset-error');
          errorElement.textContent = 'Password must be at least 8 characters long.';
          errorElement.style.display = 'block';
          return;
        }
        
        try {
          await Auth.updatePassword(newPassword);
          
          // Show success message and redirect to login
          const successElement = document.getElementById('reset-success');
          successElement.textContent = 'Password updated successfully. Redirecting to login...';
          successElement.style.display = 'block';
          document.getElementById('reset-error').style.display = 'none';
          
          // Redirect to login after a delay
          setTimeout(() => {
            window.location.href = '/login';
          }, 2000);
        } catch (error) {
          // Show error message
          const errorElement = document.getElementById('reset-error');
          errorElement.textContent = error.message || 'Failed to update password. Please try again.';
          errorElement.style.display = 'block';
          document.getElementById('reset-success').style.display = 'none';
        }
      });
    }
  }
  
  // Set up profile form
  const profileForm = document.getElementById('profile-form');
  if (profileForm) {
    // Load user data
    Auth.checkSession().then(async user => {
      if (user) {
        // Fetch profile from backend
        try {
          const token = Auth.getToken();
          const resp = await fetch('/api/v1/user_profile', {
            headers: { 'Authorization': `Bearer ${token}` }
          });
          if (resp.ok) {
            const profile = await resp.json();
            document.getElementById('profile-email').value = profile.email || user.email;
            document.getElementById('profile-name').value = profile.name || user.user_metadata?.full_name || '';
            // Optionally handle organization if present in profile
            document.getElementById('profile-organization').value = profile.organization || user.user_metadata?.organization || '';
          } else {
            // Fallback to Supabase user data
            document.getElementById('profile-email').value = user.email;
            document.getElementById('profile-name').value = user.user_metadata?.full_name || '';
            document.getElementById('profile-organization').value = user.user_metadata?.organization || '';
          }
        } catch (err) {
          // Fallback to Supabase user data
          document.getElementById('profile-email').value = user.email;
          document.getElementById('profile-name').value = user.user_metadata?.full_name || '';
          document.getElementById('profile-organization').value = user.user_metadata?.organization || '';
        }
        // Show profile content
        document.getElementById('profile-loading').style.display = 'none';
        document.getElementById('profile-content').style.display = 'block';
      } else {
        // Redirect to login if not logged in
        window.location.href = '/login';
      }
    });
    
    // Handle profile form submission
    profileForm.addEventListener('submit', async function(e) {
      e.preventDefault();
      
      const fullName = document.getElementById('profile-name').value;
      const organization = document.getElementById('profile-organization').value;
      
      try {
        // Update Supabase Auth profile
        await Auth.updateProfile({
          full_name: fullName,
          organization: organization
        });
        // Update backend user_profile
        const token = Auth.getToken();
        await fetch('/api/v1/user_profile', {
          method: 'PUT',
          headers: {
            'Content-Type': 'application/json',
            'Authorization': `Bearer ${token}`
          },
          body: JSON.stringify({
            name: fullName,
            organization: organization
          })
        });
        // Show success message
        const successElement = document.getElementById('profile-success');
        successElement.textContent = 'Profile updated successfully.';
        successElement.style.display = 'block';
        document.getElementById('profile-error').style.display = 'none';
      } catch (error) {
        // Show error message
        const errorElement = document.getElementById('profile-error');
        errorElement.textContent = error.message || 'Failed to update profile. Please try again.';
        errorElement.style.display = 'block';
        document.getElementById('profile-success').style.display = 'none';
      }
    });
  }
  
  // Set up password change form
  const passwordForm = document.getElementById('password-form');
  if (passwordForm) {
    passwordForm.addEventListener('submit', async function(e) {
      e.preventDefault();
      
      const currentPassword = document.getElementById('current-password').value;
      const newPassword = document.getElementById('new-password').value;
      const confirmNewPassword = document.getElementById('confirm-new-password').value;
      
      // Validate passwords match
      if (newPassword !== confirmNewPassword) {
        const errorElement = document.getElementById('profile-error');
        errorElement.textContent = 'New passwords do not match.';
        errorElement.style.display = 'block';
        document.getElementById('profile-success').style.display = 'none';
        return;
      }
      
      // Validate password length
      if (newPassword.length < 8) {
        const errorElement = document.getElementById('profile-error');
        errorElement.textContent = 'Password must be at least 8 characters long.';
        errorElement.style.display = 'block';
        document.getElementById('profile-success').style.display = 'none';
        return;
      }
      
      try {
        // First verify the current password by trying to sign in
        const email = document.getElementById('profile-email').value;
        await Auth.signIn(email, currentPassword);
        
        // Then update the password
        await Auth.updatePassword(newPassword);
        
        // Show success message
        const successElement = document.getElementById('profile-success');
        successElement.textContent = 'Password changed successfully.';
        successElement.style.display = 'block';
        document.getElementById('profile-error').style.display = 'none';
        
        // Clear form
        passwordForm.reset();
      } catch (error) {
        // Show error message
        const errorElement = document.getElementById('profile-error');
        errorElement.textContent = error.message || 'Failed to change password. Please verify your current password.';
        errorElement.style.display = 'block';
        document.getElementById('profile-success').style.display = 'none';
      }
    });
  }
  
  // Set up logout link
  const logoutLink = document.getElementById('logout-link');
  if (logoutLink) {
    logoutLink.addEventListener('click', async function(e) {
      e.preventDefault();
      
      try {
        await Auth.signOut();
        window.location.href = '/login';
      } catch (error) {
        console.error('Logout error:', error);
        alert('Failed to log out. Please try again.');
      }
    });
  }
});

/**
 * Update the UI based on authentication state
 * @param {Object|null} user - The authenticated user or null if not authenticated
 */
function updateAuthUI(user) {
  // Update username display if it exists
  const usernameDisplay = document.getElementById('username-display');
  if (usernameDisplay) {
    if (user) {
      usernameDisplay.textContent = user.user_metadata?.full_name || user.email;
    } else {
      usernameDisplay.textContent = 'Guest';
    }
  }
  
  // Update login/logout links
  const loginLinks = document.querySelectorAll('a[href="/login"]');
  const logoutLinks = document.querySelectorAll('#logout-link');
  const profileLinks = document.querySelectorAll('a[href="/profile"]');
  
  if (user) {
    // User is logged in
    loginLinks.forEach(link => {
      link.style.display = 'none';
    });
    
    logoutLinks.forEach(link => {
      link.style.display = 'block';
    });
    
    profileLinks.forEach(link => {
      link.style.display = 'block';
    });
  } else {
    // User is logged out
    loginLinks.forEach(link => {
      link.style.display = 'block';
    });
    
    logoutLinks.forEach(link => {
      link.style.display = 'none';
    });
    
    profileLinks.forEach(link => {
      link.style.display = 'none';
    });
  }
}
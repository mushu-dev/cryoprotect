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
      
      return data.user;
    } catch (error) {
      console.error('Sign in error:', error);
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
    signOut,
    checkSession,
    getToken
  };
})();
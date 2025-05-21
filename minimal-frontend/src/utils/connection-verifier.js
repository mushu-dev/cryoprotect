/**
 * Connection Verifier Utility
 * 
 * A utility for verifying connections between frontend and backend systems,
 * as well as database connections.
 */

import axios from 'axios';

// Configuration with environment variable support
const config = {
  apiUrl: process.env.NEXT_PUBLIC_API_URL || 'https://cryoprotect-8030e4025428.herokuapp.com/api',
  convexUrl: process.env.NEXT_PUBLIC_CONVEX_URL || 'https://upbeat-parrot-866.convex.cloud',
  useConvex: process.env.NEXT_PUBLIC_USE_CONVEX === 'true',
  timeout: 10000, // 10 seconds
};

/**
 * Verify API connection
 * @returns {Promise<{success: boolean, message: string, details: object}>}
 */
export async function verifyApiConnection() {
  try {
    // Try to connect to the health endpoint
    // This should be a public, unauthenticated endpoint
    const response = await axios.get(`${config.apiUrl}/health`, {
      timeout: config.timeout
    });
    
    if (response.status === 200) {
      // Extract useful information from the response
      const version = response.data.version || 'unknown';
      const uptime = response.data.uptime_formatted || 'unknown';
      
      return {
        success: true,
        message: `API connection successful (v${version}, uptime: ${uptime})`,
        details: response.data
      };
    } else {
      return {
        success: false,
        message: `API returned unexpected status: ${response.status}`,
        details: response.data
      };
    }
  } catch (error) {
    // Try a fallback to an alternate route in case the health endpoint is not available
    try {
      const fallbackResponse = await axios.get(`${config.apiUrl}/v1/system/status`, {
        timeout: config.timeout
      });
      
      if (fallbackResponse.status === 200) {
        return {
          success: true,
          message: 'API connection successful (fallback route)',
          details: {
            fallback: true,
            ...fallbackResponse.data
          }
        };
      }
    } catch (fallbackError) {
      // Both attempts failed
      return {
        success: false,
        message: `API connection failed (both routes): ${error.message}`,
        details: {
          primary: {
            name: error.name,
            code: error.code,
            config: error.config ? {
              url: error.config.url,
              method: error.config.method,
              timeout: error.config.timeout
            } : null
          },
          fallback: {
            name: fallbackError.name,
            code: fallbackError.code,
            config: fallbackError.config ? {
              url: fallbackError.config.url,
              method: fallbackError.config.method,
              timeout: fallbackError.config.timeout
            } : null
          }
        }
      };
    }
    
    return {
      success: false,
      message: `API connection failed: ${error.message}`,
      details: {
        name: error.name,
        code: error.code,
        config: error.config ? {
          url: error.config.url,
          method: error.config.method,
          timeout: error.config.timeout
        } : null
      }
    };
  }
}

/**
 * Verify Convex connection (if enabled)
 * @param {object} convexClient - The Convex client instance
 * @returns {Promise<{success: boolean, message: string, details: object}>}
 */
export async function verifyConvexConnection(convexClient) {
  // First, check if Convex is enabled in configuration
  if (!config.useConvex) {
    return {
      success: false,
      message: 'Convex is not enabled',
      details: { enabled: false }
    };
  }
  
  // Check if there's a client instance provided
  if (!convexClient) {
    // We'll still try to check the backend Convex connection even without a client
    try {
      // Try to connect to the Convex connection check endpoint in the API
      const response = await axios.get(`${config.apiUrl}/v1/system/convex-connection-check`, {
        timeout: config.timeout
      });
      
      if (response.status === 200) {
        return {
          success: response.data.status === 'ok',
          message: response.data.message,
          details: {
            backendCheck: true,
            connection_details: response.data.connection_details,
            tables_available: response.data.tables_available,
            client_provided: false
          }
        };
      } else {
        return {
          success: false,
          message: `Convex backend check returned unexpected status: ${response.status}`,
          details: {
            backendCheck: true,
            client_provided: false,
            ...response.data
          }
        };
      }
    } catch (error) {
      return {
        success: false,
        message: `Convex backend check failed: ${error.message}`,
        details: {
          backendCheck: true,
          client_provided: false,
          name: error.name,
          code: error.code,
          status: error.response?.status,
          data: error.response?.data
        }
      };
    }
  }
  
  // If Convex client is provided, try both client-side and backend checks
  try {
    // Try client-side check
    let clientStatus = { success: false, details: null };
    try {
      // Try to query a basic function if available
      if (convexClient.query && typeof convexClient.query === 'function') {
        // Check if the healthCheck function exists - if not, try a simple table query
        try {
          const result = await convexClient.query('system.healthCheck')();
          clientStatus = {
            success: true,
            details: result
          };
        } catch (clientError) {
          // Fallback: try a simple table query
          const result = await convexClient.query('molecules.getFirstMolecule')() || 
                         { id: 'test', success: true };
          clientStatus = {
            success: true,
            details: result
          };
        }
      } else if (convexClient.table && typeof convexClient.table === 'function') {
        // Use the table API if available (for adapter client)
        const result = await convexClient.table('molecules').select('id').limit(1).execute();
        clientStatus = {
          success: !result.error,
          details: result.data || result.error
        };
      } else {
        clientStatus = {
          success: false,
          details: { error: 'Client does not have query or table methods' }
        };
      }
    } catch (clientErr) {
      clientStatus = {
        success: false,
        details: {
          name: clientErr.name,
          message: clientErr.message,
          code: clientErr.code
        }
      };
    }
    
    // Also try server-side check
    let serverStatus = { success: false, details: null };
    try {
      const response = await axios.get(`${config.apiUrl}/v1/system/convex-connection-check`, {
        timeout: config.timeout
      });
      
      serverStatus = {
        success: response.data.status === 'ok',
        details: {
          connection_details: response.data.connection_details,
          tables_available: response.data.tables_available
        }
      };
    } catch (serverErr) {
      serverStatus = {
        success: false,
        details: {
          name: serverErr.name,
          message: serverErr.message,
          code: serverErr.code,
          status: serverErr.response?.status,
          data: serverErr.response?.data
        }
      };
    }
    
    // Combine results - consider successful if either check passes
    const isSuccess = clientStatus.success || serverStatus.success;
    
    return {
      success: isSuccess,
      message: isSuccess 
        ? 'Convex connection successful'
        : 'Convex connection failed - both client and server checks failed',
      details: {
        clientCheck: clientStatus,
        serverCheck: serverStatus
      }
    };
  } catch (error) {
    return {
      success: false,
      message: `Convex connection verification failed: ${error.message}`,
      details: {
        name: error.name,
        code: error.code,
        message: error.message
      }
    };
  }
}

/**
 * Verify authentication connection
 * @param {string} token - JWT token to verify
 * @returns {Promise<{success: boolean, message: string, details: object}>}
 */
export async function verifyAuthConnection(token) {
  if (!token) {
    return {
      success: false,
      message: 'No authentication token provided',
      details: { token: null }
    };
  }
  
  try {
    // Try to connect to the verification endpoint
    const response = await axios.get(`${config.apiUrl}/auth/verify`, {
      headers: {
        Authorization: `Bearer ${token}`
      },
      timeout: config.timeout
    });
    
    if (response.status === 200) {
      return {
        success: true,
        message: 'Authentication verified',
        details: {
          user: response.data.user,
          roles: response.data.roles,
          expires: response.data.expires
        }
      };
    } else {
      return {
        success: false,
        message: `Authentication verification returned unexpected status: ${response.status}`,
        details: response.data
      };
    }
  } catch (error) {
    return {
      success: false,
      message: `Authentication verification failed: ${error.message}`,
      details: {
        name: error.name,
        code: error.code,
        status: error.response?.status,
        data: error.response?.data
      }
    };
  }
}

/**
 * Verify service role connection
 * @returns {Promise<{success: boolean, message: string, details: object}>}
 */
export async function verifyServiceRoleConnection() {
  try {
    // Try to connect to the service role verification endpoint
    // Use the full API path as we know it from the API definition
    const response = await axios.get(`${config.apiUrl}/v1/system/service-role-check`, {
      timeout: config.timeout
    });
    
    if (response.status === 200) {
      return {
        success: true,
        message: 'Service role connection successful',
        details: response.data
      };
    } else {
      return {
        success: false,
        message: `Service role check returned unexpected status: ${response.status}`,
        details: response.data
      };
    }
  } catch (error) {
    return {
      success: false,
      message: `Service role check failed: ${error.message}`,
      details: {
        name: error.name,
        code: error.code,
        status: error.response?.status,
        data: error.response?.data
      }
    };
  }
}

/**
 * Verify database connection
 * @returns {Promise<{success: boolean, message: string, details: object}>}
 */
export async function verifyDatabaseConnection() {
  try {
    // Try to connect to the database health check endpoint
    // Use the health/database endpoint that we know exists from the API
    const response = await axios.get(`${config.apiUrl}/health/database`, {
      timeout: config.timeout
    });
    
    if (response.status === 200) {
      // Extract key information from the response
      const dbStatus = response.data.status === 'ok';
      const schemaStatus = response.data.schema_status || 'unknown';
      const integrityStatus = response.data.integrity_status || 'unknown';
      
      return {
        success: dbStatus,
        message: dbStatus 
          ? `Database connection successful (Schema: ${schemaStatus}, Integrity: ${integrityStatus})` 
          : 'Database connection issues detected',
        details: response.data
      };
    } else {
      return {
        success: false,
        message: `Database health check returned unexpected status: ${response.status}`,
        details: response.data
      };
    }
  } catch (error) {
    // Try fallback to the basic health endpoint which also checks database connectivity
    try {
      const fallbackResponse = await axios.get(`${config.apiUrl}/health`, {
        timeout: config.timeout
      });
      
      if (fallbackResponse.status === 200) {
        const dbStatus = fallbackResponse.data.database_status === 'connected';
        
        return {
          success: dbStatus,
          message: dbStatus 
            ? 'Database connection successful (basic check)' 
            : `Database disconnected: ${fallbackResponse.data.database_status}`,
          details: {
            fallback: true,
            ...fallbackResponse.data
          }
        };
      }
    } catch (fallbackError) {
      // Both checks failed
      return {
        success: false,
        message: `Database health checks failed: Primary (${error.message}), Fallback (${fallbackError.message})`,
        details: {
          primary: {
            name: error.name,
            code: error.code,
            status: error.response?.status,
            data: error.response?.data
          },
          fallback: {
            name: fallbackError.name,
            code: fallbackError.code,
            status: fallbackError.response?.status,
            data: fallbackError.response?.data
          }
        }
      };
    }
    
    return {
      success: false,
      message: `Database health check failed: ${error.message}`,
      details: {
        name: error.name,
        code: error.code,
        status: error.response?.status,
        data: error.response?.data
      }
    };
  }
}

/**
 * Run all connection verifications
 * @param {object} convexClient - The Convex client instance (optional)
 * @param {string} token - JWT token for auth verification (optional)
 * @returns {Promise<{api: object, convex: object, auth: object, serviceRole: object, database: object}>}
 */
export async function verifyAllConnections(convexClient = null, token = null) {
  // Run all verifications in parallel
  const [apiResult, convexResult, authResult, serviceRoleResult, dbResult] = await Promise.all([
    verifyApiConnection(),
    convexClient ? verifyConvexConnection(convexClient) : { success: false, message: 'Convex client not provided' },
    token ? verifyAuthConnection(token) : { success: false, message: 'Auth token not provided' },
    verifyServiceRoleConnection(),
    verifyDatabaseConnection()
  ]);
  
  return {
    api: apiResult,
    convex: convexResult,
    auth: authResult,
    serviceRole: serviceRoleResult,
    database: dbResult
  };
}

export default {
  verifyApiConnection,
  verifyConvexConnection,
  verifyAuthConnection,
  verifyServiceRoleConnection,
  verifyDatabaseConnection,
  verifyAllConnections
};
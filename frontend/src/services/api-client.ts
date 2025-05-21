import axios, { AxiosError, AxiosRequestConfig, AxiosResponse } from 'axios';
import config from '../config/config';

// Create axios instance with default config
const apiClient = axios.create({
  baseURL: process.env.NEXT_PUBLIC_API_URL || 'https://cryoprotect-8030e4025428.herokuapp.com/v1',
  timeout: 15000,
  headers: {
    'Content-Type': 'application/json',
  },
  withCredentials: true, // Include cookies in cross-site requests
});

// Add request interceptor for logging and authentication
apiClient.interceptors.request.use(
  (config) => {
    // Log API requests if enabled
    if (process.env.NEXT_PUBLIC_ENABLE_API_LOGGING === 'true') {
      console.log(`API Request: ${config.method?.toUpperCase()} ${config.url}`);
    }
    
    // Add any authentication headers if needed
    const token = typeof window !== 'undefined' ? localStorage.getItem('auth_token') : null;
    if (token) {
      config.headers.Authorization = `Bearer ${token}`;
    }
    
    // Add protection bypass token from config
    const protectionBypass = process.env.PROTECTION_BYPASS || config.frontend.protection_bypass;
    if (protectionBypass) {
      config.headers['x-protection-bypass'] = protectionBypass;
    }
    
    return config;
  },
  (error) => {
    return Promise.reject(error);
  }
);

// Add response interceptor for error handling
apiClient.interceptors.response.use(
  (response: AxiosResponse) => {
    // Log API responses if enabled
    if (process.env.NEXT_PUBLIC_ENABLE_API_LOGGING === 'true') {
      console.log(`API Response: ${response.status} ${response.config.url}`);
    }
    
    return response;
  },
  (error: AxiosError) => {
    // Handle API errors
    if (error.response) {
      // Log error responses if enabled
      if (process.env.NEXT_PUBLIC_ENABLE_API_LOGGING === 'true') {
        console.error(
          `API Error: ${error.response.status} ${error.config?.url}`,
          error.response.data
        );
      }
      
      // Handle 401 Unauthorized error and redirect to login
      if (error.response.status === 401 && typeof window !== 'undefined') {
        localStorage.removeItem('auth_token');
        // Redirect to login page if we're not already there
        if (!window.location.pathname.includes('/login')) {
          window.location.href = `/login?redirect=${encodeURIComponent(window.location.pathname)}`;
        }
      }
      
      // Handle 403 Forbidden errors
      if (error.response.status === 403) {
        console.error('Permission denied');
        // You can show a permission denied message or redirect
      }
    } else if (error.request) {
      // The request was made but no response was received
      console.error('API Network Error:', error.message);
      
      // If API is unreachable, set flag to use mock data
      if (typeof window !== 'undefined') {
        localStorage.setItem('api_connection_status', 'failed');
      }
    } else {
      // Something happened in setting up the request
      console.error('API Error:', error.message);
    }
    
    return Promise.reject(error);
  }
);

// Helper function to test API connectivity
export const testApiConnectivity = async (): Promise<boolean> => {
  try {
    const response = await apiClient.get('/health/connectivity');
    return response.data.status === 'connected';
  } catch (error) {
    console.error('API connectivity test failed:', error);
    return false;
  }
};

// Export for use in services
export default apiClient;
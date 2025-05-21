/**
 * CryoProtect Configuration
 * 
 * This file is auto-generated from the unified configuration schema.
 * It provides TypeScript types and defaults for the application configuration.
 */

// Main configuration interfaces
export interface AppConfig {
  /**
   * Application name
   * @default "CryoProtect"
   */
  name: string;
  /**
   * Application environment
   * @default "development"
   */
  environment: 'development' | 'testing' | 'staging' | 'production';
  /**
   * Enable debug mode
   * @default false
   */
  debug?: boolean;
  /**
   * Enable testing mode
   * @default false
   */
  testing?: boolean;
  /**
   * Secret key for cryptographic signing
   */
  secret_key: string;
  /**
   * Feature flag settings
   */
  feature_flags?: {
  enable_experimental_features?: boolean;
  use_mock_data?: boolean;
};
}

export interface ApiConfig {
  /**
   * API title for documentation
   * @default "CryoProtect Analyzer API"
   */
  title?: string;
  /**
   * API version
   * @default "v1"
   */
  version?: string;
  /**
   * Base URL for the API
   */
  base_url: string;
  /**
   * OpenAPI specification version
   * @default "3.0.2"
   */
  openapi_version?: string;
  /**
   * URL prefix for OpenAPI documentation
   * @default "/"
   */
  openapi_url_prefix?: string;
  /**
   * URL path for Swagger UI
   * @default "/swagger"
   */
  swagger_ui_path?: string;
  /**
   * URL for Swagger UI assets
   * @default "https://cdn.jsdelivr.net/npm/swagger-ui-dist/"
   */
  swagger_ui_url?: string;
  /**
   * CORS configuration
   */
  cors?: {
  origins?: string[];
};
}

export interface DatabaseConfig {
  /**
   * Database connection mode
   * @default "auto"
   */
  connection_mode?: 'local' | 'supabase' | 'mcp' | 'auto';
  /**
   * Local database configuration
   */
  local?: {
  host?: string;
  port?: string;
  database?: string;
  user?: string;
  password?: string;
  min_connections?: number;
  max_connections?: number;
  use_ssl?: boolean;
};
  /**
   * Supabase configuration
   */
  supabase?: {
  url: string;
  key: string;
  service_key?: string;
  project_id?: string;
  host?: string;
  port?: string;
  database?: string;
  user?: string;
  password?: string;
  ip_address?: string;
  min_connections?: number;
  max_connections?: number;
};
  /**
   * MCP configuration
   */
  mcp?: {
  project_id: string;
};
  /**
   * Connection timeout in seconds
   * @default 30
   */
  connection_timeout?: number;
  /**
   * Connection lifetime in seconds
   * @default 3600
   */
  connection_lifetime?: number;
  /**
   * Idle timeout in seconds
   * @default 300
   */
  idle_timeout?: number;
  /**
   * Application name for database connections
   * @default "CryoProtect"
   */
  application_name?: string;
}

export interface AuthConfig {
  /**
   * Authentication URL
   */
  url: string;
  /**
   * Enable NextAuth.js debug mode
   * @default false
   */
  debug?: boolean;
}

export interface LoggingConfig {
  /**
   * Logging level
   * @default "INFO"
   */
  level?: 'DEBUG' | 'INFO' | 'WARNING' | 'ERROR' | 'CRITICAL';
  /**
   * Log file name
   * @default "cryoprotectant_analysis.log"
   */
  file?: string;
  /**
   * Enable logging to file
   * @default true
   */
  log_to_file?: boolean;
}

export interface CachingConfig {
  /**
   * Cache type
   * @default "memory"
   */
  type?: 'memory' | 'redis' | 'filesystem';
  /**
   * Redis URL for caching
   */
  redis_url?: string;
  /**
   * Default cache timeout in seconds
   * @default 300
   */
  default_timeout?: number;
}

export interface Rate_limitingConfig {
  /**
   * Enable rate limiting
   * @default true
   */
  enabled?: boolean;
  /**
   * Rate limit storage URL
   * @default "memory://"
   */
  storage_url?: string;
  /**
   * Rate limiting strategy
   * @default "fixed-window"
   */
  strategy?: 'fixed-window' | 'moving-window';
  /**
   * Rate limiting identifier
   * @default "hybrid"
   */
  by?: 'ip' | 'user' | 'hybrid';
  /**
   * Enable rate limit headers
   * @default true
   */
  headers_enabled?: boolean;
  /**
   * Retry after seconds
   * @default 60
   */
  retry_after?: number;
  /**
   * Rate limit rates by role
   * @default {"admin":["5000 per day","500 per hour","100 per minute"],"premium":["2000 per day","200 per hour","40 per minute"],"basic":["1000 per day","100 per hour","20 per minute"]}
   */
  roles?: Record<string, any>;
}

export interface ChemblConfig {
  /**
   * ChEMBL API URL
   * @default "https://www.ebi.ac.uk/chembl/api/data"
   */
  api_url?: string;
  /**
   * ChEMBL API key
   */
  api_key?: string;
  /**
   * Whether ChEMBL API key is required
   * @default false
   */
  api_key_required?: boolean;
  /**
   * Delay between API requests in seconds
   * @default 0.3
   */
  api_delay?: number;
  /**
   * Path to file containing ChEMBL IDs
   */
  id_file?: string;
  /**
   * ChEMBL cache directory
   * @default "cache/chembl"
   */
  cache_dir?: string;
  /**
   * Maximum requests per second
   * @default 5
   */
  requests_per_second?: number;
  /**
   * Maximum number of retries
   * @default 5
   */
  max_retries?: number;
  /**
   * Failure threshold
   * @default 3
   */
  failure_threshold?: number;
  /**
   * Recovery timeout in seconds
   * @default 60
   */
  recovery_timeout?: number;
  /**
   * Cache time-to-live in seconds (default: 30 days)
   * @default 2592000
   */
  cache_ttl?: number;
  /**
   * Memory cache size
   * @default 1000
   */
  memory_cache_size?: number;
  /**
   * Memory usage threshold (percentage)
   * @default 80
   */
  memory_threshold?: number;
  /**
   * Batch size for ChEMBL API requests
   * @default 100
   */
  batch_size?: number;
  /**
   * Memory check frequency
   * @default 10
   */
  memory_check_frequency?: number;
  /**
   * Checkpoint directory
   * @default "checkpoints"
   */
  checkpoint_dir?: string;
}

export interface FrontendConfig {
  /**
   * Vercel deployment configuration
   */
  vercel?: {
  project_id?: string;
  org_id?: string;
  deploy_hook?: string;
};
  /**
   * Google Analytics ID
   */
  google_analytics_id?: string;
  /**
   * Sentry DSN for error tracking
   */
  sentry_dsn?: string;
  /**
   * Protection bypass token for development/preview environments
   */
  protection_bypass?: string;
}

// Default configuration values
export const defaultAppConfig: AppConfig = {
  name: "CryoProtect",
  environment: "development",
  debug: false,
  testing: false,
};

export const defaultApiConfig: ApiConfig = {
  title: "CryoProtect Analyzer API",
  version: "v1",
  openapi_version: "3.0.2",
  openapi_url_prefix: "/",
  swagger_ui_path: "/swagger",
  swagger_ui_url: "https://cdn.jsdelivr.net/npm/swagger-ui-dist/",
};

export const defaultDatabaseConfig: DatabaseConfig = {
  connection_mode: "auto",
  connection_timeout: 30,
  connection_lifetime: 3600,
  idle_timeout: 300,
  application_name: "CryoProtect",
};

export const defaultAuthConfig: AuthConfig = {
  debug: false,
};

export const defaultLoggingConfig: LoggingConfig = {
  level: "INFO",
  file: "cryoprotectant_analysis.log",
  log_to_file: true,
};

export const defaultCachingConfig: CachingConfig = {
  type: "memory",
  default_timeout: 300,
};

export const defaultRate_limitingConfig: Rate_limitingConfig = {
  enabled: true,
  storage_url: "memory://",
  strategy: "fixed-window",
  by: "hybrid",
  headers_enabled: true,
  retry_after: 60,
  roles: {
  "admin": [
    "5000 per day",
    "500 per hour",
    "100 per minute"
  ],
  "premium": [
    "2000 per day",
    "200 per hour",
    "40 per minute"
  ],
  "basic": [
    "1000 per day",
    "100 per hour",
    "20 per minute"
  ]
},
};

export const defaultChemblConfig: ChemblConfig = {
  api_url: "https://www.ebi.ac.uk/chembl/api/data",
  api_key_required: false,
  api_delay: 0.3,
  cache_dir: "cache/chembl",
  requests_per_second: 5,
  max_retries: 5,
  failure_threshold: 3,
  recovery_timeout: 60,
  cache_ttl: 2592000,
  memory_cache_size: 1000,
  memory_threshold: 80,
  batch_size: 100,
  memory_check_frequency: 10,
  checkpoint_dir: "checkpoints",
};

export const defaultFrontendConfig: FrontendConfig = {
};

export interface AppConfig {
  app: AppConfig;
  api: ApiConfig;
  database: DatabaseConfig;
  auth: AuthConfig;
  logging: LoggingConfig;
  caching: CachingConfig;
  rate_limiting: Rate_limitingConfig;
  chembl: ChemblConfig;
  frontend: FrontendConfig;
}

// Environment-specific configurations
export const developmentConfig = {
  app: {
    debug: true,
    feature_flags: {
      enable_experimental_features: true,
      use_mock_data: true,
    },
  },
  auth: {
    debug: true,
  },
  logging: {
    level: 'DEBUG',
  },
}};

export const testingConfig = {
{
  app: {
    debug: true,
    testing: true,
    feature_flags: {
      use_mock_data: true,
    },
  },
  logging: {
    level: 'DEBUG',
  },
}};

export const stagingConfig = {
{
  app: {
    debug: false,
    feature_flags: {
      enable_experimental_features: true,
    },
  },
  api: {
    cors: {
      origins: {
        0: 'https://staging.cryoprotect-analyzer.com',
      },
    },
  },
  logging: {
    level: 'INFO',
  },
}};

export const productionConfig = {
{
  app: {
    debug: false,
    feature_flags: {
      enable_experimental_features: false,
    },
  },
  api: {
    cors: {
      origins: {
        0: 'https://cryoprotect-analyzer.com',
      },
    },
  },
  logging: {
    level: 'WARNING',
  },
  caching: {
    type: 'redis',
  },
}};

// Helper functions

/**
 * Get the active configuration based on current environment
 * @returns The active configuration for the current environment
 */
export function getConfig() {
  // Determine environment
  const env = typeof process !== 'undefined' && process.env
    ? (process.env.NODE_ENV || 'development')
    : 'development';
  
  // Load base configuration
  const config = {
    app: { ...defaultAppConfig },
    api: { ...defaultApiConfig },
    database: { ...defaultDatabaseConfig },
    auth: { ...defaultAuthConfig },
    logging: { ...defaultLoggingConfig },
    caching: { ...defaultCachingConfig },
    rateLimiting: { ...defaultRateLimitingConfig },
    chembl: { ...defaultChemblConfig },
    frontend: { ...defaultFrontendConfig }
  };
  
  // Apply environment-specific configuration
  if (env === 'development') {
    mergeConfig(config, developmentConfig);
  } else if (env === 'test' || env === 'testing') {
    mergeConfig(config, testingConfig);
  } else if (env === 'staging') {
    mergeConfig(config, stagingConfig);
  } else if (env === 'production') {
    mergeConfig(config, productionConfig);
  }
  
  // Load environment variables
  loadFromEnv(config);
  
  return config;
}

/**
 * Deep merge objects
 * @param target - Target object
 * @param source - Source object
 */
function mergeConfig(target, source) {
  for (const key of Object.keys(source)) {
    if (source[key] instanceof Object && key in target) {
      mergeConfig(target[key], source[key]);
    } else {
      target[key] = source[key];
    }
  }
}

/**
 * Load configuration from environment variables
 * @param config - Configuration object to populate
 */
function loadFromEnv(config) {
  if (typeof process === 'undefined' || !process.env) {
    return;
  }
  
  // Process environment variables
  // Environment variables should be in the format: SECTION_PROPERTY
  for (const [key, value] of Object.entries(process.env)) {
    if (key.startsWith('NEXT_PUBLIC_')) {
      const envKey = key.replace('NEXT_PUBLIC_', '');
      const parts = envKey.split('_');
      
      if (parts.length >= 2) {
        // Convert to camelCase section name
        const section = parts[0].toLowerCase();
        
        // Join the rest for the property name in camelCase
        let propParts = parts.slice(1).map(p => p.toLowerCase());
        const propName = propParts[0] + propParts.slice(1).map(p => p.charAt(0).toUpperCase() + p.slice(1)).join('');
        
        // Set the value if the section and property exist
        if (config[section] && propName in config[section]) {
          try {
            // Convert value to proper type
            let typedValue = value;
            
            // Check if it's a boolean
            if (value.toLowerCase() === 'true') {
              typedValue = true;
            } else if (value.toLowerCase() === 'false') {
              typedValue = false;
            }
            // Check if it's a number
            else if (!isNaN(value) && value.trim() !== '') {
              typedValue = Number(value);
            }
            // Check if it's JSON
            else if ((value.startsWith('{') && value.endsWith('}')) || 
                     (value.startsWith('[') && value.endsWith(']'))) {
              try {
                typedValue = JSON.parse(value);
              } catch (e) {
                // Keep as string if not valid JSON
              }
            }
            
            config[section][propName] = typedValue;
          } catch (error) {
            console.error(`Error setting config from env var ${key}: ${error.message}`);
          }
        }
      }
    }
  }
}

// Export a singleton instance of the configuration
export const config = getConfig();
export default config;

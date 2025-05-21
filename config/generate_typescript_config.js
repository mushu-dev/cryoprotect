#!/usr/bin/env node
/**
 * Configuration Generator for TypeScript
 *
 * This script generates TypeScript configuration classes and types from the unified schema.json file.
 * It creates a config.ts file with proper type definitions and validation.
 */

const fs = require('fs');
const path = require('path');

/**
 * Convert a JSON schema type to TypeScript type
 * @param {Object} propertyDef - JSON schema property definition
 * @returns {string} TypeScript type
 */
function getTypeScriptType(propertyDef) {
  const type = propertyDef.type;
  
  if (type === 'string') {
    if (propertyDef.enum) {
      return propertyDef.enum.map(val => `'${val}'`).join(' | ');
    }
    return 'string';
  } else if (type === 'integer' || type === 'number') {
    return 'number';
  } else if (type === 'boolean') {
    return 'boolean';
  } else if (type === 'array') {
    const itemType = propertyDef.items ? getTypeScriptType(propertyDef.items) : 'any';
    return `${itemType}[]`;
  } else if (type === 'object') {
    if (propertyDef.properties) {
      const properties = Object.entries(propertyDef.properties)
        .map(([key, prop]) => {
          const isRequired = propertyDef.required && propertyDef.required.includes(key);
          return `  ${key}${isRequired ? '' : '?'}: ${getTypeScriptType(prop)};`;
        })
        .join('\n');
      return `{\n${properties}\n}`;
    }
    return 'Record<string, any>';
  }
  
  return 'any';
}

/**
 * Generate TypeScript interface from schema section
 * @param {string} name - Section name
 * @param {Object} section - Schema section
 * @returns {string} TypeScript interface
 */
function generateInterface(name, section) {
  const interfaceName = name.charAt(0).toUpperCase() + name.slice(1) + 'Config';
  let output = `export interface ${interfaceName} {\n`;
  
  // Add properties
  for (const [propName, propDef] of Object.entries(section.properties || {})) {
    // Add JSDoc comment if description exists
    if (propDef.description) {
      output += `  /**\n   * ${propDef.description}\n`;
      
      if (propDef.default !== undefined) {
        output += `   * @default ${JSON.stringify(propDef.default)}\n`;
      }
      
      output += `   */\n`;
    }
    
    const isRequired = section.required && section.required.includes(propName);
    const type = getTypeScriptType(propDef);
    
    output += `  ${propName}${isRequired ? '' : '?'}: ${type};\n`;
  }
  
  output += '}\n\n';
  return output;
}

/**
 * Generate TypeScript default values from schema section
 * @param {string} name - Section name
 * @param {Object} section - Schema section
 * @returns {string} TypeScript default values
 */
function generateDefaults(name, section) {
  const defaultName = `default${name.charAt(0).toUpperCase() + name.slice(1)}Config`;
  let output = `export const ${defaultName}: ${name.charAt(0).toUpperCase() + name.slice(1)}Config = {\n`;
  
  // Add default properties
  for (const [propName, propDef] of Object.entries(section.properties || {})) {
    if (propDef.default !== undefined) {
      const defaultValue = JSON.stringify(propDef.default, null, 2);
      output += `  ${propName}: ${defaultValue},\n`;
    }
  }
  
  output += '};\n\n';
  return output;
}

/**
 * Generate TypeScript environment specific configurations
 * @param {Object} envMapping - Environment mapping from schema
 * @returns {string} TypeScript environment configurations
 */
function generateEnvironmentConfigs(envMapping) {
  let output = '';
  
  for (const [envName, values] of Object.entries(envMapping)) {
    const configName = `${envName}Config`;
    
    output += `export const ${configName} = {\n`;
    
    // Process dot notation paths to create nested objects
    const processedValues = {};
    
    for (const [dotPath, value] of Object.entries(values)) {
      const parts = dotPath.split('.');
      let current = processedValues;
      
      // Build the nested structure
      for (let i = 0; i < parts.length - 1; i++) {
        const part = parts[i];
        if (!current[part]) {
          current[part] = {};
        }
        current = current[part];
      }
      
      // Set the value at the leaf
      current[parts[parts.length - 1]] = value;
    }
    
    // Convert to TypeScript object notation
    const stringifyObject = (obj, indent = 2) => {
      const spaces = ' '.repeat(indent);
      let result = '{\n';
      
      for (const [key, value] of Object.entries(obj)) {
        if (typeof value === 'object' && value !== null) {
          result += `${spaces}${key}: ${stringifyObject(value, indent + 2)},\n`;
        } else if (typeof value === 'string') {
          result += `${spaces}${key}: '${value}',\n`;
        } else {
          result += `${spaces}${key}: ${value},\n`;
        }
      }
      
      result += ' '.repeat(indent - 2) + '}';
      return result;
    };
    
    output += stringifyObject(processedValues);
    output += '};\n\n';
  }
  
  return output;
}

/**
 * Generate TypeScript configuration helper functions
 * @returns {string} TypeScript helper functions
 */
function generateHelpers() {
  return `
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
            console.error(\`Error setting config from env var \${key}: \${error.message}\`);
          }
        }
      }
    }
  }
}
`;
}

/**
 * Generate complete TypeScript configuration
 * @param {Object} schema - Configuration schema
 * @returns {string} Complete TypeScript configuration
 */
function generateTypeScriptConfig(schema) {
  let output = `/**
 * CryoProtect Configuration
 * 
 * This file is auto-generated from the unified configuration schema.
 * It provides TypeScript types and defaults for the application configuration.
 */

// Main configuration interfaces
`;

  // Generate interfaces for each section
  for (const [sectionName, section] of Object.entries(schema.properties)) {
    output += generateInterface(sectionName, section);
  }
  
  // Generate default values
  output += '// Default configuration values\n';
  for (const [sectionName, section] of Object.entries(schema.properties)) {
    output += generateDefaults(sectionName, section);
  }
  
  // Generate root config interface
  output += 'export interface AppConfig {\n';
  for (const sectionName of Object.keys(schema.properties)) {
    const interfaceName = sectionName.charAt(0).toUpperCase() + sectionName.slice(1) + 'Config';
    output += `  ${sectionName}: ${interfaceName};\n`;
  }
  output += '}\n\n';
  
  // Generate environment-specific configurations
  output += '// Environment-specific configurations\n';
  output += generateEnvironmentConfigs(schema.environmentMapping || {});
  
  // Generate helper functions
  output += '// Helper functions\n';
  output += generateHelpers();
  
  // Export a default config instance
  output += `
// Export a singleton instance of the configuration
export const config = getConfig();
export default config;
`;

  return output;
}

/**
 * Main function to generate TypeScript configuration
 */
function main() {
  const schemaPath = path.resolve(__dirname, 'schema.json');
  const outputPath = path.resolve(__dirname, '../frontend/src/config/config.ts');
  const outputDir = path.dirname(outputPath);
  
  // Ensure schema file exists
  if (!fs.existsSync(schemaPath)) {
    console.error(`Error: Schema file not found at ${schemaPath}`);
    process.exit(1);
  }
  
  // Ensure output directory exists
  if (!fs.existsSync(outputDir)) {
    fs.mkdirSync(outputDir, { recursive: true });
  }
  
  // Read schema
  const schema = JSON.parse(fs.readFileSync(schemaPath, 'utf-8'));
  
  // Generate TypeScript configuration
  const tsConfig = generateTypeScriptConfig(schema);
  
  // Write output
  fs.writeFileSync(outputPath, tsConfig);
  
  console.log(`Generated TypeScript configuration file at ${outputPath}`);
}

// Run main function
if (require.main === module) {
  main();
}

module.exports = {
  generateTypeScriptConfig
};
#!/usr/bin/env node
/**
 * scan_js_eslint.js - JavaScript security scanner using ESLint
 * 
 * This script scans JavaScript code for security vulnerabilities using ESLint
 * with security plugins. It can be run as a standalone script or integrated
 * into CI/CD pipelines.
 * 
 * Usage:
 *   node scan_js_eslint.js [options]
 * 
 * Options:
 *   --path PATH           Path to scan (default: ./static/js)
 *   --format FORMAT       Output format (json, stylish, html) (default: stylish)
 *   --output FILE         Output file (default: eslint-results.{format})
 *   --exit-on-critical    Exit with code 1 if critical issues found
 *   --config FILE         Custom ESLint config file
 *   --help                Show this help message and exit
 */

const fs = require('fs');
const path = require('path');
const { execSync, spawnSync } = require('child_process');

// Parse command line arguments
const args = parseArgs();

// Main function
async function main() {
  try {
    // Check if ESLint and required plugins are installed
    if (!checkEslintInstalled()) {
      if (!await installEslint()) {
        process.exit(1);
      }
    }

    // Create default ESLint config if needed
    if (!args.config) {
      createDefaultConfig();
    }

    // Run ESLint scan
    const [scanSuccess, outputFile] = runEslintScan(args);
    
    if (scanSuccess && outputFile) {
      const summary = generateSummary(outputFile, args.format, scanSuccess);
      
      if (args.exitOnCritical && checkForCriticalIssues(outputFile, args.format)) {
        console.log("Exiting with code 1 due to critical issues found.");
        process.exit(1);
      }
    }
    
    process.exit(scanSuccess ? 0 : 1);
  } catch (error) {
    console.error(`Error: ${error.message}`);
    process.exit(1);
  }
}

// Parse command line arguments
function parseArgs() {
  const args = {
    path: './static/js',
    format: 'stylish',
    output: null,
    exitOnCritical: false,
    config: null,
    help: false
  };

  for (let i = 2; i < process.argv.length; i++) {
    const arg = process.argv[i];
    
    if (arg === '--path' && i + 1 < process.argv.length) {
      args.path = process.argv[++i];
    } else if (arg === '--format' && i + 1 < process.argv.length) {
      args.format = process.argv[++i];
    } else if (arg === '--output' && i + 1 < process.argv.length) {
      args.output = process.argv[++i];
    } else if (arg === '--exit-on-critical') {
      args.exitOnCritical = true;
    } else if (arg === '--config' && i + 1 < process.argv.length) {
      args.config = process.argv[++i];
    } else if (arg === '--help') {
      args.help = true;
      showHelp();
      process.exit(0);
    }
  }

  // Set default output file if not provided
  if (!args.output) {
    const ext = args.format === 'html' ? 'html' : (args.format === 'json' ? 'json' : 'txt');
    args.output = `eslint-results.${ext}`;
  }

  return args;
}

// Show help message
function showHelp() {
  console.log(`
Usage: node scan_js_eslint.js [options]

Options:
  --path PATH           Path to scan (default: ./static/js)
  --format FORMAT       Output format (json, stylish, html) (default: stylish)
  --output FILE         Output file (default: eslint-results.{format})
  --exit-on-critical    Exit with code 1 if critical issues found
  --config FILE         Custom ESLint config file
  --help                Show this help message and exit
  `);
}

// Check if ESLint and required plugins are installed
function checkEslintInstalled() {
  try {
    execSync('npx eslint --version', { stdio: 'ignore' });
    return true;
  } catch (error) {
    return false;
  }
}

// Install ESLint and required plugins
async function installEslint() {
  console.log("ESLint not found. Installing ESLint and security plugins...");
  
  try {
    execSync('npm install --save-dev eslint eslint-plugin-security eslint-plugin-no-unsanitized', { stdio: 'inherit' });
    console.log("ESLint and security plugins installed successfully.");
    return true;
  } catch (error) {
    console.error("Failed to install ESLint and security plugins. Please install them manually:");
    console.error("npm install --save-dev eslint eslint-plugin-security eslint-plugin-no-unsanitized");
    return false;
  }
}

// Create default ESLint config if needed
function createDefaultConfig() {
  const configFile = '.eslintrc.json';
  
  if (!fs.existsSync(configFile)) {
    const config = {
      "env": {
        "browser": true,
        "es2021": true,
        "node": true
      },
      "extends": "eslint:recommended",
      "plugins": ["security", "no-unsanitized"],
      "parserOptions": {
        "ecmaVersion": "latest",
        "sourceType": "module"
      },
      "rules": {
        // Security rules
        "security/detect-buffer-noassert": "error",
        "security/detect-child-process": "warn",
        "security/detect-disable-mustache-escape": "error",
        "security/detect-eval-with-expression": "error",
        "security/detect-new-buffer": "error",
        "security/detect-no-csrf-before-method-override": "error",
        "security/detect-non-literal-fs-filename": "warn",
        "security/detect-non-literal-regexp": "warn",
        "security/detect-non-literal-require": "warn",
        "security/detect-object-injection": "warn",
        "security/detect-possible-timing-attacks": "warn",
        "security/detect-pseudoRandomBytes": "error",
        "security/detect-unsafe-regex": "error",
        
        // No unsanitized rules
        "no-unsanitized/method": "error",
        "no-unsanitized/property": "error"
      }
    };
    
    fs.writeFileSync(configFile, JSON.stringify(config, null, 2));
    console.log(`Created default ESLint config: ${configFile}`);
  }
}

// Run ESLint scan
function runEslintScan(args) {
  console.log(`Running ESLint scan on ${args.path}...`);
  
  const eslintArgs = [
    args.path,
    '--format', args.format,
    '--output-file', args.output
  ];
  
  if (args.config) {
    eslintArgs.push('--config', args.config);
  }
  
  try {
    // ESLint returns exit code 1 if it finds issues, which we don't want to treat as an error
    const result = spawnSync('npx', ['eslint', ...eslintArgs], { encoding: 'utf8' });
    
    if (result.status !== 0 && result.status !== 1) {
      console.error(`Error running ESLint: ${result.stderr}`);
      return [false, null];
    }
    
    console.log(`Scan completed. Results saved to ${args.output}`);
    return [true, args.output];
  } catch (error) {
    console.error(`Error running ESLint: ${error.message}`);
    return [false, null];
  }
}

// Check for critical issues in the scan results
function checkForCriticalIssues(outputFile, format) {
  if (!fs.existsSync(outputFile)) {
    console.log(`Output file ${outputFile} not found.`);
    return false;
  }
  
  if (format === 'json') {
    try {
      const results = JSON.parse(fs.readFileSync(outputFile, 'utf8'));
      
      let criticalCount = 0;
      for (const result of results) {
        for (const message of result.messages) {
          // Check for security-related rules with error severity
          if (message.severity === 2 && 
              (message.ruleId?.startsWith('security/') || 
               message.ruleId?.startsWith('no-unsanitized/'))) {
            criticalCount++;
          }
        }
      }
      
      if (criticalCount > 0) {
        console.log(`Found ${criticalCount} critical security issues!`);
        return true;
      }
    } catch (error) {
      console.error(`Error parsing JSON results: ${error.message}`);
    }
  }
  
  return false;
}

// Generate a summary of the scan results
function generateSummary(outputFile, format, scanSuccess) {
  const summary = {
    scanner: "eslint-security",
    timestamp: new Date().toISOString(),
    success: scanSuccess,
    output_file: outputFile,
    output_format: format
  };
  
  if (scanSuccess && format === 'json' && fs.existsSync(outputFile)) {
    try {
      const results = JSON.parse(fs.readFileSync(outputFile, 'utf8'));
      
      // Count issues by severity and rule
      const issueCounts = {
        error: 0,
        warning: 0,
        total: 0
      };
      
      const securityIssues = {
        error: 0,
        warning: 0,
        total: 0
      };
      
      for (const result of results) {
        for (const message of result.messages) {
          issueCounts.total++;
          
          if (message.severity === 2) {
            issueCounts.error++;
          } else if (message.severity === 1) {
            issueCounts.warning++;
          }
          
          // Count security-specific issues
          if (message.ruleId?.startsWith('security/') || 
              message.ruleId?.startsWith('no-unsanitized/')) {
            securityIssues.total++;
            
            if (message.severity === 2) {
              securityIssues.error++;
            } else if (message.severity === 1) {
              securityIssues.warning++;
            }
          }
        }
      }
      
      summary.issues = issueCounts;
      summary.security_issues = securityIssues;
    } catch (error) {
      console.error(`Error generating summary: ${error.message}`);
    }
  }
  
  // Save summary to file
  const timestamp = new Date().toISOString().replace(/[:.]/g, '-');
  const summaryFile = `eslint-summary-${timestamp}.json`;
  
  fs.writeFileSync(summaryFile, JSON.stringify(summary, null, 2));
  console.log(`Summary saved to ${summaryFile}`);
  
  return summary;
}

// Run the main function
main();
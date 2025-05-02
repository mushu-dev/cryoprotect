#!/usr/bin/env node

/**
 * CryoProtect Analyzer - Performance Indexes Verification Utility
 * 
 * This script verifies that the performance indexes have been successfully created.
 * It requires the Supabase CLI to be installed.
 */

const { execSync } = require('child_process');
const readline = require('readline');

const rl = readline.createInterface({
  input: process.stdin,
  output: process.stdout
});

// Colors for console output
const colors = {
  reset: '\x1b[0m',
  bright: '\x1b[1m',
  dim: '\x1b[2m',
  red: '\x1b[31m',
  green: '\x1b[32m',
  yellow: '\x1b[33m',
  blue: '\x1b[34m',
  cyan: '\x1b[36m'
};

// Print banner
console.log(`
${colors.cyan}${colors.bright}=================================================
 CryoProtect Analyzer - Performance Indexes Verification
=================================================
${colors.reset}
This utility will verify that the performance indexes
have been successfully created in your Supabase database.
`);

// Check if Supabase CLI is installed
try {
  execSync('supabase --version', { stdio: 'ignore' });
  console.log(`${colors.green}✓ Supabase CLI is installed${colors.reset}`);
} catch (error) {
  console.error(`
${colors.red}✗ Supabase CLI is not installed or not in PATH${colors.reset}

Please install the Supabase CLI first:
  npm install -g supabase
  
Then run this script again.
`);
  process.exit(1);
}

// Ask for Supabase project reference
rl.question(`
${colors.yellow}Enter your Supabase project reference:${colors.reset} `, (projectRef) => {
  if (!projectRef.trim()) {
    console.error(`
${colors.red}✗ Project reference is required${colors.reset}

You can find your project reference in the Supabase dashboard:
1. Go to https://app.supabase.io/
2. Select your project
3. Go to Project Settings > General
4. Copy the "Reference ID"
`);
    rl.close();
    process.exit(1);
  }

  console.log(`
${colors.yellow}Linking to Supabase project: ${projectRef}${colors.reset}
`);

  try {
    // Link to Supabase project
    execSync(`supabase link --project-ref ${projectRef}`, { stdio: 'inherit' });
    
    console.log(`
${colors.green}✓ Successfully linked to Supabase project${colors.reset}
    
${colors.yellow}Verifying performance indexes...${colors.reset}
`);

    // SQL queries to check for the existence of each index
    const verificationQueries = [
      {
        name: "Composite index on predictions(mixture_id, property_type_id)",
        query: "SELECT EXISTS (SELECT 1 FROM pg_indexes WHERE indexname = 'idx_predictions_mixture_property')"
      },
      {
        name: "Text search index on molecule names",
        query: "SELECT EXISTS (SELECT 1 FROM pg_indexes WHERE indexname = 'idx_molecule_name_trgm')"
      },
      {
        name: "Text search index on mixture names",
        query: "SELECT EXISTS (SELECT 1 FROM pg_indexes WHERE indexname = 'idx_mixture_name_trgm')"
      },
      {
        name: "pg_trgm extension",
        query: "SELECT EXISTS (SELECT 1 FROM pg_extension WHERE extname = 'pg_trgm')"
      }
    ];

    // Create a temporary file with the verification queries
    const fs = require('fs');
    const path = require('path');
    const tempFile = path.join(__dirname, 'temp_verification.sql');
    
    let allIndexesExist = true;
    
    for (const query of verificationQueries) {
      // Write the query to a temporary file
      fs.writeFileSync(tempFile, query.query);
      
      // Execute the query
      const result = execSync(`supabase db execute --file ${tempFile}`, { encoding: 'utf8' });
      
      // Check if the index exists
      const exists = result.trim().toLowerCase().includes('true');
      
      if (exists) {
        console.log(`${colors.green}✓ ${query.name} exists${colors.reset}`);
      } else {
        console.log(`${colors.red}✗ ${query.name} does not exist${colors.reset}`);
        allIndexesExist = false;
      }
    }
    
    // Clean up the temporary file
    fs.unlinkSync(tempFile);
    
    if (allIndexesExist) {
      console.log(`
${colors.green}${colors.bright}✓ All performance indexes have been successfully created!${colors.reset}

Your CryoProtect Analyzer database has been optimized with the following indexes:
- Composite index on predictions(mixture_id, property_type_id)
- Text search index on molecule names using pg_trgm
- Text search index on mixture names using pg_trgm

These indexes should significantly improve query performance.
`);
    } else {
      console.log(`
${colors.yellow}${colors.bright}⚠ Some performance indexes are missing.${colors.reset}

Please run the apply_performance_indexes script to create all the required indexes.
`);
    }
  } catch (error) {
    console.error(`
${colors.red}✗ An error occurred while verifying the indexes${colors.reset}

Error details:
${error.message}

Please check your project reference and try again.
`);
  }
  
  rl.close();
});
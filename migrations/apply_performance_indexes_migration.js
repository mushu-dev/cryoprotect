#!/usr/bin/env node

/**
 * CryoProtect Analyzer - Performance Indexes Migration Utility
 * 
 * This script helps apply the performance indexes migration to a Supabase project.
 * It requires the Supabase CLI to be installed.
 */

const { execSync } = require('child_process');
const fs = require('fs');
const path = require('path');
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
 CryoProtect Analyzer - Performance Indexes Migration
=================================================
${colors.reset}
This utility will apply the performance indexes migration
to your Supabase project to optimize database performance.
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

// Check if migration file exists
const migrationFile = path.join(__dirname, '010_performance_indexes.sql');
if (!fs.existsSync(migrationFile)) {
  console.error(`
${colors.red}✗ Migration file not found: ${migrationFile}${colors.reset}

Please make sure the migration file exists in the same directory as this script.
`);
  process.exit(1);
}

console.log(`${colors.green}✓ Migration file found: ${migrationFile}${colors.reset}`);

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
    
${colors.yellow}Applying performance indexes migration...${colors.reset}
`);

    // Apply migration directly using SQL
    const migrationContent = fs.readFileSync(migrationFile, 'utf8');
    
    // Create a temporary file with the migration content
    const tempFile = path.join(__dirname, 'temp_migration.sql');
    fs.writeFileSync(tempFile, migrationContent);
    
    // Execute the SQL file against the Supabase database
    execSync(`supabase db execute --file ${tempFile}`, { stdio: 'inherit' });
    
    // Clean up the temporary file
    fs.unlinkSync(tempFile);
    
    console.log(`
${colors.green}${colors.bright}✓ Performance indexes migration successfully applied!${colors.reset}

Your CryoProtect Analyzer database has been optimized with the following indexes:
- Composite index on predictions(mixture_id, property_type_id)
- Text search index on molecule names using pg_trgm
- Text search index on mixture names using pg_trgm

These indexes should significantly improve query performance for:
- Queries that filter on both mixture_id and property_type_id
- Text search operations on molecule and mixture names

${colors.dim}For more information, see the README_Performance_Testing.md file.${colors.reset}
`);
  } catch (error) {
    console.error(`
${colors.red}✗ An error occurred while applying the migration${colors.reset}

Error details:
${error.message}

Please check your project reference and try again.
You can also apply the migration manually using the Supabase dashboard:
1. Go to https://app.supabase.io/
2. Select your project
3. Go to SQL Editor
4. Copy the contents of 010_performance_indexes.sql
5. Paste into the SQL Editor and run the query
`);
  }
  
  rl.close();
});
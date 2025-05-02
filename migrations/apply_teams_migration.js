/**
 * Apply Teams Migration Script
 * 
 * This script applies the teams schema migration to the Supabase database.
 * It reads the SQL file and executes it using the Supabase client.
 */

require('dotenv').config();
const fs = require('fs');
const path = require('path');
const { createClient } = require('@supabase/supabase-js');

// Supabase configuration
const supabaseUrl = process.env.SUPABASE_URL;
const supabaseKey = process.env.SUPABASE_KEY;

if (!supabaseUrl || !supabaseKey) {
  console.error('Error: SUPABASE_URL and SUPABASE_KEY must be set in .env file');
  process.exit(1);
}

// Create Supabase client
const supabase = createClient(supabaseUrl, supabaseKey);

// Path to the migration file
const migrationFilePath = path.join(__dirname, '003_teams_schema.sql');

// Read the migration file
fs.readFile(migrationFilePath, 'utf8', async (err, sql) => {
  if (err) {
    console.error(`Error reading migration file: ${err.message}`);
    process.exit(1);
  }

  console.log('Applying teams schema migration...');

  try {
    // Execute the SQL
    const { error } = await supabase.rpc('exec_sql', { sql });

    if (error) {
      console.error(`Error applying migration: ${error.message}`);
      process.exit(1);
    }

    console.log('Teams schema migration applied successfully!');
  } catch (error) {
    console.error(`Error applying migration: ${error.message}`);
    process.exit(1);
  }
});
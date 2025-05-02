/**
 * Apply Export and Sharing Schema Migration
 * 
 * This script applies the export and sharing schema migration to the Supabase database.
 */

require('dotenv').config();
const fs = require('fs');
const path = require('path');
const { createClient } = require('@supabase/supabase-js');

// Get Supabase credentials from environment variables
const supabaseUrl = process.env.SUPABASE_URL;
const supabaseKey = process.env.SUPABASE_SERVICE_KEY;

if (!supabaseUrl || !supabaseKey) {
    console.error('Error: SUPABASE_URL and SUPABASE_SERVICE_KEY environment variables are required.');
    process.exit(1);
}

// Create Supabase client
const supabase = createClient(supabaseUrl, supabaseKey);

// Read migration SQL file
const migrationFile = path.join(__dirname, '004_export_sharing_schema.sql');
const migrationSql = fs.readFileSync(migrationFile, 'utf8');

// Apply migration
async function applyMigration() {
    console.log('Applying export and sharing schema migration...');
    
    try {
        // Execute the SQL migration
        const { error } = await supabase.rpc('exec_sql', { sql: migrationSql });
        
        if (error) {
            console.error('Error applying migration:', error);
            process.exit(1);
        }
        
        console.log('Migration applied successfully!');
        
        // Verify tables were created
        const tables = ['shares', 'report_templates', 'visualization_templates', 'export_logs', 'share_access_logs'];
        
        for (const table of tables) {
            const { data, error } = await supabase
                .from('information_schema.tables')
                .select('table_name')
                .eq('table_name', table)
                .eq('table_schema', 'public');
                
            if (error) {
                console.error(`Error verifying table ${table}:`, error);
                continue;
            }
            
            if (data && data.length > 0) {
                console.log(`Table ${table} created successfully.`);
            } else {
                console.warn(`Warning: Table ${table} may not have been created.`);
            }
        }
        
        // Verify functions were created
        const functions = ['update_share_access_count', 'is_share_expired', 'get_shared_item_data'];
        
        for (const func of functions) {
            const { data, error } = await supabase
                .from('information_schema.routines')
                .select('routine_name')
                .eq('routine_name', func)
                .eq('routine_schema', 'public');
                
            if (error) {
                console.error(`Error verifying function ${func}:`, error);
                continue;
            }
            
            if (data && data.length > 0) {
                console.log(`Function ${func} created successfully.`);
            } else {
                console.warn(`Warning: Function ${func} may not have been created.`);
            }
        }
        
    } catch (err) {
        console.error('Unexpected error:', err);
        process.exit(1);
    }
}

// Run the migration
applyMigration();
/**
 * Sample Data Population Script
 * 
 * This script populates the Convex database with sample data for development and testing.
 * It uses the sampleDataGenerator utilities to create realistic scientific data.
 */

// Import Convex client and necessary modules
const { ConvexClient } = require("convex/browser");
const { generateCompleteSampleDataset } = require("../utils/sampleDataGenerator");
const { createClient } = require("@supabase/supabase-js");

// Default to development environment variables
// Consider using dotenv or similar for more robust env handling
const CONVEX_URL = process.env.CONVEX_URL || "https://eager-cow-39.convex.cloud";
const SUPABASE_URL = process.env.SUPABASE_URL || "https://your-project.supabase.co";
const SUPABASE_KEY = process.env.SUPABASE_KEY || "your-anon-key";

// Parse command line arguments
const args = process.argv.slice(2);
const sampleSize = parseInt(args[0]) || 10;
const skipConfirmation = args.includes("--force");

// Create Convex client
const convexClient = new ConvexClient(CONVEX_URL);

/**
 * Main execution function
 */
async function main() {
  try {
    if (!skipConfirmation) {
      // Check with user before proceeding
      console.log("\n⚠️  WARNING: This will populate your Convex database with sample data.");
      console.log(`⚠️  Number of records per type to generate: ${sampleSize}`);
      console.log(`⚠️  Database URL: ${CONVEX_URL}`);
      console.log("⚠️  This operation may overwrite existing data.\n");
      
      const response = await promptUser("Do you want to proceed? (yes/no): ");
      if (response.toLowerCase() !== "yes") {
        console.log("Operation cancelled by user.");
        process.exit(0);
      }
    }
    
    console.log("Starting sample data population...");
    
    // Generate data using our utility
    const results = await generateCompleteSampleDataset(convexClient.api, sampleSize);
    
    // Log results
    console.log("\n✅ Sample data population complete!");
    console.log(`✅ Generated ${Object.values(results).flat().length} total records`);
    console.log("✅ Generated the following record counts:");
    for (const [type, ids] of Object.entries(results)) {
      console.log(`   - ${type}: ${ids.length} records`);
    }
    
    process.exit(0);
  } catch (error) {
    console.error("❌ Error populating sample data:", error);
    process.exit(1);
  }
}

/**
 * Simple utility to prompt user for input
 */
function promptUser(question) {
  const readline = require("readline").createInterface({
    input: process.stdin,
    output: process.stdout
  });
  
  return new Promise((resolve) => {
    readline.question(question, (answer) => {
      readline.close();
      resolve(answer);
    });
  });
}

/**
 * Alternative: Run data migration from Supabase
 */
async function runMigration() {
  const { runFullMigration } = require("../utils/migrationUtility");
  
  try {
    console.log("\n⚠️  WARNING: This will migrate data from Supabase to Convex.");
    console.log(`⚠️  Supabase URL: ${SUPABASE_URL}`);
    console.log(`⚠️  Convex URL: ${CONVEX_URL}`);
    console.log("⚠️  This operation may overwrite existing data in Convex.\n");
    
    const response = await promptUser("Do you want to proceed with migration? (yes/no): ");
    if (response.toLowerCase() !== "yes") {
      console.log("Migration cancelled by user.");
      return;
    }
    
    console.log("Starting migration from Supabase to Convex...");
    
    // Run the migration
    const ctx = await runFullMigration(
      convexClient,
      SUPABASE_URL,
      SUPABASE_KEY,
      { logLevel: 'info', batchSize: 50 }
    );
    
    console.log("\n✅ Migration complete!");
    
    // Generate report
    let report = "\n=== Migration Report ===\n\n";
    for (const tableName in ctx.stats) {
      const stats = ctx.stats[tableName];
      report += `${tableName}: ${stats.successfulRecords}/${stats.totalRecords} records migrated`;
      report += ` (${stats.failedRecords} failures)\n`;
    }
    
    console.log(report);
  } catch (error) {
    console.error("❌ Error during migration:", error);
  }
}

// Execute the main function if this script is run directly
if (require.main === module) {
  if (args.includes("--migrate")) {
    runMigration();
  } else {
    main();
  }
}

module.exports = {
  main,
  runMigration
};
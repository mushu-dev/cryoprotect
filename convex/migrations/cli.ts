/**
 * Command-line interface for managing Convex migrations.
 * 
 * This module provides a CLI for applying and rolling back migrations,
 * as well as checking migration status.
 */

import { ConvexClient } from "convex/browser";
import fs from "fs";
import path from "path";
import { fileURLToPath } from "url";
import { program } from "commander";

// Get the directory name
const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

// Import environment variables
import dotenv from "dotenv";
dotenv.config();

// Default paths
const DEFAULT_MIGRATIONS_DIR = path.join(__dirname, "migrations");

// Configure the CLI
program
  .name("convex-migrations")
  .description("Manage Convex database migrations for CryoProtect")
  .version("1.0.0");

// Status command
program
  .command("status")
  .description("Show the status of migrations")
  .option("--format <format>", "Output format (text or json)", "text")
  .option("--dir <dir>", "Directory containing migration files")
  .action(async (options) => {
    await runStatus(options);
  });

// Apply command
program
  .command("apply")
  .description("Apply pending migrations")
  .option("--target <version>", "Target migration version")
  .option("--dry-run", "Perform a dry run without making changes")
  .option("--dir <dir>", "Directory containing migration files")
  .action(async (options) => {
    await runApply(options);
  });

// Rollback command
program
  .command("rollback")
  .description("Roll back applied migrations")
  .option("--target <version>", "Target migration version to roll back to")
  .option("--dry-run", "Perform a dry run without making changes")
  .action(async (options) => {
    await runRollback(options);
  });

// Create command
program
  .command("create <name>")
  .description("Create a new migration file")
  .option("--dir <dir>", "Directory to create the migration file in")
  .action(async (name, options) => {
    await createMigration(name, options);
  });

// Parse command-line arguments
program.parse(process.argv);

// Function to initialize the Convex client
function initializeClient() {
  const convexUrl = process.env.CONVEX_URL;
  
  if (!convexUrl) {
    console.error("Error: CONVEX_URL environment variable is required");
    process.exit(1);
  }
  
  return new ConvexClient(convexUrl);
}

// Function to run the status command
async function runStatus(options) {
  try {
    const client = initializeClient();
    
    // Get migration status
    const status = await client.query("migrations.runner:getMigrationStatus", {
      migrationsDir: options.dir || DEFAULT_MIGRATIONS_DIR
    });
    
    if (options.format === "json") {
      console.log(JSON.stringify(status, null, 2));
    } else {
      console.log("\nMigration Status:\n");
      
      const migrations = Object.values(status.migrations).sort((a, b) => 
        a.version.localeCompare(b.version)
      );
      
      if (migrations.length === 0) {
        console.log("No migrations found.");
      } else {
        for (const migration of migrations) {
          const appliedText = migration.applied ? "Applied" : "Pending";
          const appliedAtText = migration.applied ? ` at ${new Date(migration.appliedAt).toISOString()}` : "";
          const sourceText = migration.source === "database" ? "[DB]" : "[File]";
          
          console.log(`${migration.version}: ${migration.name} - ${appliedText}${appliedAtText} ${sourceText}`);
          
          if (migration.error) {
            console.log(`  Error: ${migration.error}`);
          }
        }
      }
      
      console.log();
    }
  } catch (error) {
    console.error("Error checking migration status:", error);
    process.exit(1);
  }
}

// Function to run the apply command
async function runApply(options) {
  try {
    const client = initializeClient();
    
    console.log(`Applying migrations${options.target ? ` to ${options.target}` : ""} (${options.dryRun ? "dry run" : "live"})`);
    
    // Apply migrations
    const result = await client.action("migrations.runner:applyMigrations", {
      targetVersion: options.target,
      dryRun: options.dryRun === true,
      migrationsDir: options.dir || DEFAULT_MIGRATIONS_DIR
    });
    
    if (options.dryRun) {
      console.log("\nDry run - would apply the following migrations:\n");
    } else {
      console.log("\nApplied migrations:\n");
    }
    
    for (const migration of result.applied) {
      console.log(`${migration.version}: ${migration.name} - ${migration.status}`);
      
      if (migration.error) {
        console.log(`  Error: ${migration.error}`);
      }
    }
    
    console.log(`\n${result.message}`);
  } catch (error) {
    console.error("Error applying migrations:", error);
    process.exit(1);
  }
}

// Function to run the rollback command
async function runRollback(options) {
  try {
    const client = initializeClient();
    
    console.log(`Rolling back migrations${options.target ? ` to ${options.target}` : ""} (${options.dryRun ? "dry run" : "live"})`);
    
    // Roll back migrations
    const result = await client.action("migrations.runner:rollbackMigrations", {
      targetVersion: options.target,
      dryRun: options.dryRun === true
    });
    
    if (options.dryRun) {
      console.log("\nDry run - would roll back the following migrations:\n");
    } else {
      console.log("\nRolled back migrations:\n");
    }
    
    for (const migration of result.rolledBack) {
      console.log(`${migration.version}: ${migration.name} - ${migration.status}`);
      
      if (migration.error) {
        console.log(`  Error: ${migration.error}`);
      }
    }
    
    console.log(`\n${result.message}`);
  } catch (error) {
    console.error("Error rolling back migrations:", error);
    process.exit(1);
  }
}

// Function to create a new migration file
async function createMigration(name, options) {
  try {
    const migrationsDir = options.dir || DEFAULT_MIGRATIONS_DIR;
    const examplesDir = path.join(migrationsDir, "examples");
    
    // Ensure directories exist
    if (!fs.existsSync(migrationsDir)) {
      fs.mkdirSync(migrationsDir, { recursive: true });
    }
    
    if (!fs.existsSync(examplesDir)) {
      fs.mkdirSync(examplesDir, { recursive: true });
    }
    
    // Get the next version number
    const files = fs.readdirSync(examplesDir);
    const versions = files
      .filter(file => file.match(/^\d{3}_.*\.json$/))
      .map(file => parseInt(file.substring(0, 3), 10));
    
    const nextVersion = versions.length > 0 
      ? Math.max(...versions) + 1 
      : 1;
    
    const paddedVersion = nextVersion.toString().padStart(3, "0");
    const snakeCaseName = name
      .replace(/[^a-zA-Z0-9]+/g, "_")
      .replace(/([A-Z])/g, "_$1")
      .toLowerCase()
      .replace(/^_/, "");
    
    const filename = `${paddedVersion}_${snakeCaseName}.json`;
    const filePath = path.join(examplesDir, filename);
    
    // Create the migration file
    const migration = {
      version: paddedVersion,
      name: snakeCaseName,
      description: name,
      schema: {},
      data: [],
      indexes: []
    };
    
    fs.writeFileSync(filePath, JSON.stringify(migration, null, 2));
    
    console.log(`Created migration file: ${filePath}`);
  } catch (error) {
    console.error("Error creating migration:", error);
    process.exit(1);
  }
}

// If no command is provided, show help
if (!process.argv.slice(2).length) {
  program.outputHelp();
}
#!/usr/bin/env python3
"""
verify_database_integrity_enhanced.py

This script performs a comprehensive verification of the CryoProtect database integrity,
checking foreign key constraints, data consistency, and data quality across all tables.

Usage:
    python verify_database_integrity_enhanced.py [--output=<filename>] [--verbose]
    
Options:
    --output=<filename>    Output file for the JSON report (default: database_integrity_report.json)
    --verbose              Enable verbose output with detailed progress information
"""

import os
import json
import argparse
import logging
import psycopg2
import psycopg2.extras
from datetime import datetime
from dotenv import load_dotenv
from typing import Dict, List, Any, Tuple, Optional, Set, Union
import concurrent.futures
import time
import re

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

# Database connection parameters
DB_HOST = os.getenv("SUPABASE_DB_HOST")
DB_PORT = os.getenv("SUPABASE_DB_PORT")
DB_NAME = os.getenv("SUPABASE_DB_NAME")
DB_USER = os.getenv("SUPABASE_DB_USER")
DB_PASSWORD = os.getenv("SUPABASE_DB_PASSWORD")
DB_SERVICE_ROLE = os.getenv("SUPABASE_SERVICE_ROLE")

# Chemical data validation patterns
SMILES_PATTERN = r'^([^J][a-z0-9@+\-\[\]\(\)\\\/%=#$]{6,})$'
INCHI_PATTERN = r'^InChI=1S?\/[A-Za-z0-9.+\-\\\/()]+$'
INCHIKEY_PATTERN = r'^[A-Z]{14}-[A-Z]{10}-[A-Z]$'

# Known valid sources for data
VALID_DATA_SOURCES = {
    'pubchem', 'chembl', 'manual', 'prediction', 'experiment', 
    'reference', 'calculated', 'import', 'api', 'auto'
}

class DatabaseIntegrityVerifier:
    """Class to verify the integrity of the CryoProtect database."""
    
    def __init__(self, verbose=False):
        """Initialize the verifier with DB connection and settings."""
        self.verbose = verbose
        self.conn = None
        self.cursor = None
        self.report = {
            "timestamp": datetime.now().isoformat(),
            "summary": {
                "tables_checked": 0,
                "rows_checked": 0,
                "foreign_keys_checked": 0,
                "issues_found": 0,
            },
            "issues": [],
            "data_counts": {},
            "status": "PENDING"
        }
        self.tables_info = {}
        self.foreign_keys = []
        
    def log(self, message, level="info"):
        """Log message if verbose mode is enabled."""
        if self.verbose:
            getattr(logger, level)(message)
            
    def connect_to_database(self) -> bool:
        """Establish connection to the database."""
        try:
            self.conn = psycopg2.connect(
                host=DB_HOST,
                port=DB_PORT,
                dbname=DB_NAME,
                user=DB_USER,
                password=DB_PASSWORD
            )
            self.cursor = self.conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
            self.log("Connected to database successfully")
            return True
        except Exception as e:
            logger.error(f"Database connection error: {e}")
            return False

    def close_connection(self):
        """Close the database connection."""
        if self.cursor:
            self.cursor.close()
        if self.conn:
            self.conn.close()
        self.log("Database connection closed")

    def get_schema_info(self):
        """Retrieve schema information about tables, columns, and constraints."""
        try:
            # Get all tables in public schema
            self.cursor.execute("""
                SELECT table_name 
                FROM information_schema.tables 
                WHERE table_schema = 'public'
                ORDER BY table_name
            """)
            tables = [row[0] for row in self.cursor.fetchall()]
            
            # Get columns for each table
            for table_name in tables:
                self.cursor.execute("""
                    SELECT column_name, data_type, is_nullable
                    FROM information_schema.columns
                    WHERE table_schema = 'public' AND table_name = %s
                    ORDER BY ordinal_position
                """, (table_name,))
                columns = self.cursor.fetchall()
                
                self.tables_info[table_name] = {
                    "columns": [{
                        "name": col[0],
                        "type": col[1],
                        "nullable": col[2] == "YES"
                    } for col in columns],
                    "row_count": 0
                }
            
            # Get foreign key constraints
            self.cursor.execute("""
                SELECT
                    tc.constraint_name,
                    tc.table_name,
                    kcu.column_name,
                    ccu.table_name AS foreign_table_name,
                    ccu.column_name AS foreign_column_name
                FROM
                    information_schema.table_constraints AS tc
                    JOIN information_schema.key_column_usage AS kcu
                      ON tc.constraint_name = kcu.constraint_name
                      AND tc.table_schema = kcu.table_schema
                    JOIN information_schema.constraint_column_usage AS ccu
                      ON ccu.constraint_name = tc.constraint_name
                      AND ccu.table_schema = tc.table_schema
                WHERE
                    tc.constraint_type = 'FOREIGN KEY'
                    AND tc.table_schema = 'public'
                ORDER BY tc.table_name, kcu.column_name
            """)
            self.foreign_keys = self.cursor.fetchall()
            
            self.log(f"Retrieved schema information for {len(tables)} tables and {len(self.foreign_keys)} foreign key constraints")
            return True
        except Exception as e:
            logger.error(f"Error retrieving schema info: {e}")
            return False

    def count_rows(self):
        """Count rows in each table."""
        for table_name in self.tables_info:
            try:
                self.cursor.execute(f"SELECT COUNT(*) FROM public.{table_name}")
                count = self.cursor.fetchone()[0]
                self.tables_info[table_name]["row_count"] = count
                self.report["data_counts"][table_name] = count
                self.log(f"Table {table_name}: {count} rows")
            except Exception as e:
                logger.error(f"Error counting rows in {table_name}: {e}")
                self.add_issue("error", f"Error counting rows in {table_name}", str(e))
        
        # Calculate totals
        total_rows = sum(info["row_count"] for info in self.tables_info.values())
        self.report["summary"]["tables_checked"] = len(self.tables_info)
        self.report["summary"]["rows_checked"] = total_rows
        self.log(f"Total rows across all tables: {total_rows}")

    def verify_foreign_keys(self):
        """Verify all foreign key relationships."""
        start_time = time.time()
        fk_issues = []
        
        for constraint in self.foreign_keys:
            constraint_name, table_name, column_name, ref_table, ref_column = constraint
            self.log(f"Checking foreign key: {table_name}.{column_name} -> {ref_table}.{ref_column}")
            
            # Skip if table has no data
            if self.tables_info[table_name]["row_count"] == 0:
                self.log(f"Skipping {table_name}.{column_name} check - table is empty")
                continue
            
            try:
                # Get distinct values for faster checking
                self.cursor.execute(f"""
                    SELECT DISTINCT {column_name} 
                    FROM public.{table_name} 
                    WHERE {column_name} IS NOT NULL
                """)
                fk_values = [row[0] for row in self.cursor.fetchall()]
                
                if not fk_values:
                    self.log(f"No non-null values in {table_name}.{column_name} to verify")
                    continue
                
                # Process in batches to avoid massive IN() clauses
                batch_size = 1000
                invalid_values = []
                
                for i in range(0, len(fk_values), batch_size):
                    batch = fk_values[i:i+batch_size]
                    if not batch:
                        continue
                        
                    # Format values properly according to data type
                    if isinstance(batch[0], str):
                        # For string types, need to escape and quote
                        placeholders = ', '.join(["%s"] * len(batch))
                        query = f"""
                            SELECT a.{column_name} 
                            FROM (VALUES {placeholders}) AS a({column_name})
                            WHERE a.{column_name} NOT IN (
                                SELECT {ref_column} FROM public.{ref_table}
                            )
                        """
                        self.cursor.execute(query, batch)
                    else:
                        # For numeric/uuid types
                        placeholders = ','.join(["%s"] * len(batch))
                        query = f"""
                            SELECT a.value 
                            FROM (SELECT unnest(%s::uuid[]) AS value) a
                            WHERE a.value NOT IN (
                                SELECT {ref_column} FROM public.{ref_table}
                            )
                        """
                        self.cursor.execute(query, (batch,))
                    
                    invalid_batch = [row[0] for row in self.cursor.fetchall()]
                    invalid_values.extend(invalid_batch)
                
                if invalid_values:
                    # Limit the number of values shown in report
                    display_values = invalid_values[:10]
                    more_text = f" (and {len(invalid_values) - 10} more)" if len(invalid_values) > 10 else ""
                    issue_detail = f"Values: {', '.join(str(v) for v in display_values)}{more_text}"
                    issue = {
                        "type": "foreign_key_violation",
                        "table": table_name,
                        "column": column_name,
                        "foreign_table": ref_table,
                        "foreign_column": ref_column,
                        "invalid_count": len(invalid_values),
                        "details": issue_detail
                    }
                    fk_issues.append(issue)
                    self.add_issue("error", f"Foreign key violation in {table_name}.{column_name}", issue_detail)
                    
            except Exception as e:
                logger.error(f"Error checking foreign key {table_name}.{column_name}: {e}")
                self.add_issue("error", f"Error checking foreign key {table_name}.{column_name}", str(e))
        
        self.report["summary"]["foreign_keys_checked"] = len(self.foreign_keys)
        self.log(f"Foreign key verification completed in {time.time() - start_time:.2f} seconds")
        return fk_issues

    def verify_required_fields(self):
        """Check that non-nullable fields have values in all tables."""
        for table_name, info in self.tables_info.items():
            if info["row_count"] == 0:
                continue
                
            required_columns = [col["name"] for col in info["columns"] if not col["nullable"]]
            if not required_columns:
                continue
                
            self.log(f"Checking required fields in {table_name}")
            
            for column in required_columns:
                try:
                    self.cursor.execute(f"""
                        SELECT COUNT(*) FROM public.{table_name}
                        WHERE {column} IS NULL
                    """)
                    null_count = self.cursor.fetchone()[0]
                    
                    if null_count > 0:
                        self.add_issue(
                            "error", 
                            f"NULL values in required field", 
                            f"Table {table_name} has {null_count} NULL values in required column {column}"
                        )
                except Exception as e:
                    logger.error(f"Error checking required field {table_name}.{column}: {e}")
                    self.add_issue("error", f"Error checking required field {table_name}.{column}", str(e))

    def verify_chemical_data_quality(self):
        """Verify quality of chemical data in molecules table."""
        try:
            # Skip if molecules table is empty
            if self.tables_info.get("molecules", {}).get("row_count", 0) == 0:
                self.log("Skipping chemical data quality check - molecules table is empty")
                return []
                
            # 1. Check SMILES format
            self.cursor.execute("""
                SELECT id, name, smiles FROM public.molecules
                WHERE smiles IS NOT NULL 
                AND smiles NOT SIMILAR TO '[A-Za-z0-9@+\\-\\[\\]\\(\\)\\\\\\/%=#$]+'
            """)
            invalid_smiles = self.cursor.fetchall()
            if invalid_smiles:
                invalid_count = len(invalid_smiles)
                examples = [f"{row['name']} ({row['smiles']})" for row in invalid_smiles[:5]]
                self.add_issue(
                    "warning",
                    "Invalid SMILES notation",
                    f"Found {invalid_count} molecules with potentially invalid SMILES. Examples: {', '.join(examples)}"
                )
            
            # 2. Check InChI format
            self.cursor.execute("""
                SELECT id, name, inchi FROM public.molecules
                WHERE inchi IS NOT NULL 
                AND inchi NOT LIKE 'InChI=1%'
            """)
            invalid_inchi = self.cursor.fetchall()
            if invalid_inchi:
                invalid_count = len(invalid_inchi)
                examples = [f"{row['name']} ({row['inchi']})" for row in invalid_inchi[:5]]
                self.add_issue(
                    "warning",
                    "Invalid InChI notation",
                    f"Found {invalid_count} molecules with potentially invalid InChI. Examples: {', '.join(examples)}"
                )
                
            # 3. Check InChIKey format
            self.cursor.execute("""
                SELECT id, name, inchikey FROM public.molecules
                WHERE inchikey IS NOT NULL 
                AND inchikey NOT SIMILAR TO '[A-Z]{14}-[A-Z]{10}-[A-Z]'
            """)
            invalid_inchikey = self.cursor.fetchall()
            if invalid_inchikey:
                invalid_count = len(invalid_inchikey)
                examples = [f"{row['name']} ({row['inchikey']})" for row in invalid_inchikey[:5]]
                self.add_issue(
                    "warning",
                    "Invalid InChIKey format",
                    f"Found {invalid_count} molecules with invalid InChIKey format. Examples: {', '.join(examples)}"
                )
                
            # 4. Check molecular weight is within reasonable range
            self.cursor.execute("""
                SELECT id, name, molecular_weight FROM public.molecules
                WHERE molecular_weight IS NOT NULL 
                AND (molecular_weight <= 0 OR molecular_weight > 5000)
            """)
            invalid_weight = self.cursor.fetchall()
            if invalid_weight:
                invalid_count = len(invalid_weight)
                examples = [f"{row['name']} ({row['molecular_weight']})" for row in invalid_weight[:5]]
                self.add_issue(
                    "warning",
                    "Suspicious molecular weight",
                    f"Found {invalid_count} molecules with suspicious molecular weight. Examples: {', '.join(examples)}"
                )
                
            # 5. Check for duplicate InChIKeys (possible duplicate molecules)
            self.cursor.execute("""
                SELECT inchikey, COUNT(*) as count, array_agg(name) as names
                FROM public.molecules
                WHERE inchikey IS NOT NULL
                GROUP BY inchikey
                HAVING COUNT(*) > 1
            """)
            duplicate_inchikeys = self.cursor.fetchall()
            if duplicate_inchikeys:
                duplicate_count = len(duplicate_inchikeys)
                examples = [f"{row['inchikey']} ({row['count']} entries: {', '.join(row['names'][:3])}" 
                           for row in duplicate_inchikeys[:5]]
                self.add_issue(
                    "warning",
                    "Duplicate molecules by InChIKey",
                    f"Found {duplicate_count} InChIKeys with multiple molecule entries. Examples: {', '.join(examples)}"
                )
                
        except Exception as e:
            logger.error(f"Error verifying chemical data quality: {e}")
            self.add_issue("error", "Error verifying chemical data quality", str(e))

    def verify_mixture_integrity(self):
        """Verify mixture data integrity and composition."""
        try:
            # Skip if mixtures table is empty
            if self.tables_info.get("mixtures", {}).get("row_count", 0) == 0:
                self.log("Skipping mixture integrity check - mixtures table is empty")
                return
                
            # 1. Check mixtures without components
            self.cursor.execute("""
                SELECT m.id, m.name
                FROM public.mixtures m
                LEFT JOIN public.mixture_components c ON m.id = c.mixture_id
                WHERE c.id IS NULL
            """)
            empty_mixtures = self.cursor.fetchall()
            if empty_mixtures:
                empty_count = len(empty_mixtures)
                examples = [f"{row['name']}" for row in empty_mixtures[:5]]
                self.add_issue(
                    "warning",
                    "Mixtures without components",
                    f"Found {empty_count} mixtures without any components. Examples: {', '.join(examples)}"
                )
                
            # 2. Check concentration units consistency
            self.cursor.execute("""
                SELECT mixture_id, COUNT(DISTINCT concentration_unit) as unit_count, 
                       array_agg(DISTINCT concentration_unit) as units
                FROM public.mixture_components
                GROUP BY mixture_id
                HAVING COUNT(DISTINCT concentration_unit) > 1
            """)
            inconsistent_units = self.cursor.fetchall()
            if inconsistent_units:
                inconsistent_count = len(inconsistent_units)
                examples = [f"Mixture ID {row['mixture_id']} uses units: {', '.join(row['units'])}" 
                           for row in inconsistent_units[:5]]
                self.add_issue(
                    "warning",
                    "Inconsistent concentration units",
                    f"Found {inconsistent_count} mixtures with inconsistent concentration units. Examples: {', '.join(examples)}"
                )
                
            # 3. Check total concentration (for percentage-based units, should sum to ~100%)
            self.cursor.execute("""
                SELECT m.id, m.name, SUM(c.concentration) as total_concentration, 
                       c.concentration_unit
                FROM public.mixtures m
                JOIN public.mixture_components c ON m.id = c.mixture_id
                WHERE c.concentration_unit IN ('%', 'percent', 'pct')
                GROUP BY m.id, m.name, c.concentration_unit
                HAVING SUM(c.concentration) < 95 OR SUM(c.concentration) > 105
            """)
            invalid_percentages = self.cursor.fetchall()
            if invalid_percentages:
                invalid_count = len(invalid_percentages)
                examples = [f"{row['name']} (total: {row['total_concentration']}%)" 
                           for row in invalid_percentages[:5]]
                self.add_issue(
                    "warning",
                    "Invalid percentage concentrations",
                    f"Found {invalid_count} mixtures with percentage concentrations not summing to ~100%. Examples: {', '.join(examples)}"
                )
                
        except Exception as e:
            logger.error(f"Error verifying mixture integrity: {e}")
            self.add_issue("error", "Error verifying mixture integrity", str(e))

    def verify_property_data_consistency(self):
        """Verify consistency of property data between related tables."""
        try:
            # Check if properties are consistent between molecules and molecular_properties
            if (self.tables_info.get("molecules", {}).get("row_count", 0) > 0 and 
                self.tables_info.get("molecular_properties", {}).get("row_count", 0) > 0):
                
                # Check molecules with JSONB properties but no matching entries in molecular_properties
                self.cursor.execute("""
                    SELECT COUNT(*) FROM public.molecules m
                    WHERE m.properties != '{}'::jsonb 
                    AND NOT EXISTS (
                        SELECT 1 FROM public.molecular_properties mp 
                        WHERE mp.molecule_id = m.id
                    )
                """)
                inconsistent_count = self.cursor.fetchone()[0]
                if inconsistent_count > 0:
                    self.add_issue(
                        "warning",
                        "Inconsistent property data",
                        f"Found {inconsistent_count} molecules with properties JSON but no entries in molecular_properties table"
                    )
                    
            # Check experiment properties consistency
            if (self.tables_info.get("experiments", {}).get("row_count", 0) > 0 and 
                self.tables_info.get("experiment_properties", {}).get("row_count", 0) > 0):
                
                # Check experiments with no properties
                self.cursor.execute("""
                    SELECT COUNT(*) FROM public.experiments e
                    WHERE NOT EXISTS (
                        SELECT 1 FROM public.experiment_properties ep 
                        WHERE ep.experiment_id = e.id
                    )
                """)
                no_props_count = self.cursor.fetchone()[0]
                if no_props_count > 0:
                    self.add_issue(
                        "warning",
                        "Experiments without properties",
                        f"Found {no_props_count} experiments with no entries in experiment_properties table"
                    )
                
        except Exception as e:
            logger.error(f"Error verifying property data consistency: {e}")
            self.add_issue("error", "Error verifying property data consistency", str(e))

    def verify_team_project_access(self):
        """Verify team and project access control integrity."""
        try:
            # Check if users belong to teams
            self.cursor.execute("""
                SELECT COUNT(*) FROM public.user_profile 
                WHERE team_id IS NULL
            """)
            no_team_count = self.cursor.fetchone()[0]
            if no_team_count > 0:
                self.add_issue(
                    "warning",
                    "Users without team",
                    f"Found {no_team_count} users without an assigned team"
                )
                
            # Check if projects have associated teams
            self.cursor.execute("""
                SELECT COUNT(*) FROM public.projects
                WHERE team_id IS NULL
            """)
            no_team_projects = self.cursor.fetchone()[0]
            if no_team_projects > 0:
                self.add_issue(
                    "warning",
                    "Projects without team",
                    f"Found {no_team_projects} projects without an associated team"
                )
                
            # Check if there are team members for each team
            self.cursor.execute("""
                SELECT t.id, t.name
                FROM public.teams t
                LEFT JOIN public.team_members tm ON t.id = tm.team_id
                WHERE tm.id IS NULL
            """)
            empty_teams = self.cursor.fetchall()
            if empty_teams:
                empty_count = len(empty_teams)
                examples = [f"{row['name']}" for row in empty_teams[:5]]
                self.add_issue(
                    "warning",
                    "Teams without members",
                    f"Found {empty_count} teams without any members. Examples: {', '.join(examples)}"
                )
                
        except Exception as e:
            logger.error(f"Error verifying team/project access: {e}")
            self.add_issue("error", "Error verifying team/project access", str(e))

    def add_issue(self, severity, title, description):
        """Add an issue to the report."""
        issue = {
            "severity": severity,
            "title": title,
            "description": description,
            "timestamp": datetime.now().isoformat()
        }
        self.report["issues"].append(issue)
        self.report["summary"]["issues_found"] += 1
        
        if severity == "error":
            logger.error(f"{title}: {description}")
        else:
            logger.warning(f"{title}: {description}")

    def execute(self) -> Dict[str, Any]:
        """Execute all verification checks and generate a report."""
        start_time = time.time()
        
        logger.info("Starting CryoProtect Database Integrity Verification")
        
        # Connect to the database
        if not self.connect_to_database():
            self.report["status"] = "FAILED"
            self.report["error"] = "Failed to connect to database"
            return self.report
            
        try:
            # Get schema information
            if not self.get_schema_info():
                self.report["status"] = "FAILED"
                self.report["error"] = "Failed to retrieve schema information"
                return self.report
                
            # Count rows in all tables
            self.count_rows()
            
            # Perform verification checks
            self.verify_foreign_keys()
            self.verify_required_fields()
            self.verify_chemical_data_quality()
            self.verify_mixture_integrity()
            self.verify_property_data_consistency()
            self.verify_team_project_access()
            
            # Set final status
            self.report["status"] = "PASS" if self.report["summary"]["issues_found"] == 0 else "FAIL"
            self.report["execution_time_seconds"] = round(time.time() - start_time, 2)
            
            # Print summary
            logger.info("=" * 60)
            logger.info("CryoProtect Database Integrity Verification Summary")
            logger.info("=" * 60)
            logger.info(f"Status: {self.report['status']}")
            logger.info(f"Tables checked: {self.report['summary']['tables_checked']}")
            logger.info(f"Rows checked: {self.report['summary']['rows_checked']}")
            logger.info(f"Foreign key constraints checked: {self.report['summary']['foreign_keys_checked']}")
            logger.info(f"Issues found: {self.report['summary']['issues_found']}")
            logger.info(f"Execution time: {self.report['execution_time_seconds']} seconds")
            logger.info("=" * 60)
            
            return self.report
            
        except Exception as e:
            logger.error(f"Error during verification: {e}")
            self.report["status"] = "FAILED"
            self.report["error"] = str(e)
            return self.report
            
        finally:
            self.close_connection()


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Verify CryoProtect database integrity")
    parser.add_argument("--output", default="database_integrity_report.json", 
                        help="Output file for the JSON report")
    parser.add_argument("--verbose", action="store_true", 
                        help="Enable verbose output")
    return parser.parse_args()


def main():
    """Main entry point."""
    args = parse_args()
    
    # Initialize and run the verifier
    verifier = DatabaseIntegrityVerifier(verbose=args.verbose)
    report = verifier.execute()
    
    # Save report to file
    with open(args.output, "w") as f:
        json.dump(report, f, indent=2)
    
    logger.info(f"Verification complete. Report saved to {args.output}")
    
    # Return exit code based on status
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    exit_code = main()
    exit(exit_code)
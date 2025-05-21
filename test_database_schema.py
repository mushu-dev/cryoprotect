#!/usr/bin/env python3
"""
CryoProtect v2 - Database Schema and Data Integrity Tests

This script tests the database schema and data integrity of the CryoProtect v2 database.
It verifies table structures, foreign key relationships, data consistency, and validation constraints.
"""

import os
import sys
import json
import uuid
import logging
from datetime import datetime
from typing import Dict, List, Any, Optional, Tuple

# Add the parent directory to the path so we can import the app modules
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("database_schema_test.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Import the Supabase client
try:
    from config import SUPABASE_URL, SUPABASE_KEY, SUPABASE_SERVICE_KEY
    from supabase import create_client, Client
except ImportError as e:
    logger.error(f"Failed to import required modules: {str(e)}")
    logger.error("Please ensure that the required modules are installed.")
    sys.exit(1)

class DatabaseSchemaTest:
    """Test class for database schema and data integrity."""

    def __init__(self):
        """Initialize the test class."""
        self.supabase: Optional[Client] = None
        self.test_results: Dict[str, Any] = {
            "status": "Not Started",
            "total_tests": 0,
            "passed_tests": 0,
            "failed_tests": 0,
            "skipped_tests": 0,
            "test_cases": []
        }
        self.expected_tables = [
            "molecules",
            "property_types",
            "molecular_properties",
            "mixtures",
            "mixture_components",
            "calculation_methods",
            "predictions",
            "experiments",
            "molecule_proteins",
            "molecule_experiments",
            "projects",
            "teams",
            "user_profile"
        ]
        self.expected_foreign_keys = [
            ("molecular_properties", "molecule_id", "molecules", "id"),
            ("molecular_properties", "property_type_id", "property_types", "id"),
            ("mixture_components", "mixture_id", "mixtures", "id"),
            ("mixture_components", "molecule_id", "molecules", "id"),
            ("predictions", "molecule_id", "molecules", "id"),
            ("predictions", "mixture_id", "mixtures", "id"),
            ("predictions", "property_type_id", "property_types", "id"),
            ("predictions", "calculation_method_id", "calculation_methods", "id"),
            ("experiments", "mixture_id", "mixtures", "id"),
            ("experiments", "molecule_id", "molecules", "id"),
            ("experiments", "property_type_id", "property_types", "id"),
            ("molecule_proteins", "molecule_id", "molecules", "id"),
            ("molecule_proteins", "protein_id", "proteins", "id"),
            ("molecule_experiments", "molecule_id", "molecules", "id"),
            ("molecule_experiments", "experiment_id", "experiments", "id"),
            ("projects", "team_id", "teams", "id"),
            ("user_profile", "team_id", "teams", "id")
        ]

    def connect_to_database(self) -> bool:
        """Connect to the Supabase database."""
        try:
            self.supabase = create_client(SUPABASE_URL, SUPABASE_SERVICE_KEY)
            logger.info("Connected to Supabase database")
            return True
        except Exception as e:
            logger.error(f"Failed to connect to Supabase database: {str(e)}")
            return False

    def run_tests(self) -> Dict[str, Any]:
        """Run all database schema and data integrity tests."""
        self.test_results["status"] = "Running"
        self.test_results["start_time"] = datetime.now().isoformat()

        # Connect to the database
        if not self.connect_to_database():
            self.test_results["status"] = "Failed"
            self.test_results["end_time"] = datetime.now().isoformat()
            return self.test_results

        # Run the tests
        self.test_table_existence()
        self.test_table_structures()
        self.test_foreign_key_relationships()
        self.test_data_consistency()
        self.test_data_validation()
        self.test_rls_enablement()

        # Calculate test results
        self.test_results["total_tests"] = len(self.test_results["test_cases"])
        self.test_results["passed_tests"] = sum(1 for tc in self.test_results["test_cases"] if tc["status"] == "Passed")
        self.test_results["failed_tests"] = sum(1 for tc in self.test_results["test_cases"] if tc["status"] == "Failed")
        self.test_results["skipped_tests"] = sum(1 for tc in self.test_results["test_cases"] if tc["status"] == "Skipped")
        
        if self.test_results["failed_tests"] == 0:
            self.test_results["status"] = "Passed"
        else:
            self.test_results["status"] = "Failed"
        
        self.test_results["end_time"] = datetime.now().isoformat()
        
        return self.test_results

    def add_test_result(self, test_id: str, test_name: str, status: str, message: str) -> None:
        """Add a test result to the test results."""
        self.test_results["test_cases"].append({
            "id": test_id,
            "name": test_name,
            "status": status,
            "message": message
        })
        
        if status == "Passed":
            logger.info(f"Test {test_id} - {test_name}: PASSED")
        elif status == "Failed":
            logger.error(f"Test {test_id} - {test_name}: FAILED - {message}")
        else:
            logger.warning(f"Test {test_id} - {test_name}: SKIPPED - {message}")

    def test_table_existence(self) -> None:
        """Test that all expected tables exist in the database."""
        try:
            # Get the list of tables in the public schema
            response = self.supabase.table("information_schema.tables").select("table_name").eq("table_schema", "public").execute()
            tables = [table["table_name"] for table in response.data]
            
            # Check that all expected tables exist
            missing_tables = [table for table in self.expected_tables if table not in tables]
            
            if not missing_tables:
                self.add_test_result("DB-1.1", "Table Existence", "Passed", "All expected tables exist")
            else:
                self.add_test_result("DB-1.1", "Table Existence", "Failed", f"Missing tables: {', '.join(missing_tables)}")
        except Exception as e:
            self.add_test_result("DB-1.1", "Table Existence", "Failed", f"Error checking table existence: {str(e)}")

    def test_table_structures(self) -> None:
        """Test that all tables have the correct structure."""
        try:
            # For each table, check its structure
            for table in self.expected_tables:
                # Skip tables that don't exist
                response = self.supabase.table("information_schema.tables").select("table_name").eq("table_schema", "public").eq("table_name", table).execute()
                if not response.data:
                    self.add_test_result(f"DB-1.2.{table}", f"Table Structure - {table}", "Skipped", f"Table {table} does not exist")
                    continue
                
                # Get the columns for the table
                response = self.supabase.table("information_schema.columns").select("column_name, data_type, is_nullable").eq("table_schema", "public").eq("table_name", table).execute()
                columns = {col["column_name"]: {"data_type": col["data_type"], "is_nullable": col["is_nullable"]} for col in response.data}
                
                # Check for required columns
                required_columns = ["id"]
                if table not in ["user_profile"]:  # Tables that might not have created_by
                    required_columns.append("created_by")
                
                missing_columns = [col for col in required_columns if col not in columns]
                
                if not missing_columns:
                    self.add_test_result(f"DB-1.2.{table}", f"Table Structure - {table}", "Passed", f"Table {table} has all required columns")
                else:
                    self.add_test_result(f"DB-1.2.{table}", f"Table Structure - {table}", "Failed", f"Table {table} is missing columns: {', '.join(missing_columns)}")
        except Exception as e:
            self.add_test_result("DB-1.2", "Table Structures", "Failed", f"Error checking table structures: {str(e)}")

    def test_foreign_key_relationships(self) -> None:
        """Test that all foreign key relationships are correctly defined."""
        try:
            # Get all foreign key constraints
            query = """
            SELECT
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
            """
            
            # Execute the query using the Supabase client
            response = self.supabase.rpc("exec_sql", {"query": query}).execute()
            
            # Process the results
            foreign_keys = []
            for row in response.data:
                foreign_keys.append((
                    row["table_name"],
                    row["column_name"],
                    row["foreign_table_name"],
                    row["foreign_column_name"]
                ))
            
            # Check that all expected foreign keys exist
            missing_fks = []
            for fk in self.expected_foreign_keys:
                if fk not in foreign_keys:
                    missing_fks.append(f"{fk[0]}.{fk[1]} -> {fk[2]}.{fk[3]}")
            
            if not missing_fks:
                self.add_test_result("DB-2", "Foreign Key Relationships", "Passed", "All expected foreign key relationships exist")
            else:
                self.add_test_result("DB-2", "Foreign Key Relationships", "Failed", f"Missing foreign key relationships: {', '.join(missing_fks)}")
        except Exception as e:
            self.add_test_result("DB-2", "Foreign Key Relationships", "Failed", f"Error checking foreign key relationships: {str(e)}")

    def test_data_consistency(self) -> None:
        """Test that data is consistent across related tables."""
        try:
            # Check for orphaned records in molecular_properties
            query = """
            SELECT COUNT(*) AS orphaned_count
            FROM molecular_properties mp
            LEFT JOIN molecules m ON mp.molecule_id = m.id
            WHERE m.id IS NULL
            """
            response = self.supabase.rpc("exec_sql", {"query": query}).execute()
            orphaned_count = response.data[0]["orphaned_count"] if response.data else 0
            
            if orphaned_count == 0:
                self.add_test_result("DB-3.1", "Orphaned Records - molecular_properties", "Passed", "No orphaned records found")
            else:
                self.add_test_result("DB-3.1", "Orphaned Records - molecular_properties", "Failed", f"Found {orphaned_count} orphaned records")
            
            # Check for orphaned records in mixture_components
            query = """
            SELECT COUNT(*) AS orphaned_count
            FROM mixture_components mc
            LEFT JOIN mixtures m ON mc.mixture_id = m.id
            WHERE m.id IS NULL
            """
            response = self.supabase.rpc("exec_sql", {"query": query}).execute()
            orphaned_count = response.data[0]["orphaned_count"] if response.data else 0
            
            if orphaned_count == 0:
                self.add_test_result("DB-3.2", "Orphaned Records - mixture_components", "Passed", "No orphaned records found")
            else:
                self.add_test_result("DB-3.2", "Orphaned Records - mixture_components", "Failed", f"Found {orphaned_count} orphaned records")
            
            # Check for orphaned records in predictions
            query = """
            SELECT COUNT(*) AS orphaned_count
            FROM predictions p
            LEFT JOIN mixtures m ON p.mixture_id = m.id
            LEFT JOIN molecules mol ON p.molecule_id = mol.id
            WHERE (p.mixture_id IS NOT NULL AND m.id IS NULL)
               OR (p.molecule_id IS NOT NULL AND mol.id IS NULL)
            """
            response = self.supabase.rpc("exec_sql", {"query": query}).execute()
            orphaned_count = response.data[0]["orphaned_count"] if response.data else 0
            
            if orphaned_count == 0:
                self.add_test_result("DB-3.3", "Orphaned Records - predictions", "Passed", "No orphaned records found")
            else:
                self.add_test_result("DB-3.3", "Orphaned Records - predictions", "Failed", f"Found {orphaned_count} orphaned records")
            
            # Check for orphaned records in experiments
            query = """
            SELECT COUNT(*) AS orphaned_count
            FROM experiments e
            LEFT JOIN mixtures m ON e.mixture_id = m.id
            LEFT JOIN molecules mol ON e.molecule_id = mol.id
            WHERE (e.mixture_id IS NOT NULL AND m.id IS NULL)
               OR (e.molecule_id IS NOT NULL AND mol.id IS NULL)
            """
            response = self.supabase.rpc("exec_sql", {"query": query}).execute()
            orphaned_count = response.data[0]["orphaned_count"] if response.data else 0
            
            if orphaned_count == 0:
                self.add_test_result("DB-3.4", "Orphaned Records - experiments", "Passed", "No orphaned records found")
            else:
                self.add_test_result("DB-3.4", "Orphaned Records - experiments", "Failed", f"Found {orphaned_count} orphaned records")
        except Exception as e:
            self.add_test_result("DB-3", "Data Consistency", "Failed", f"Error checking data consistency: {str(e)}")

    def test_data_validation(self) -> None:
        """Test that data validation constraints are enforced."""
        try:
            # Test that we can't insert a molecule with a duplicate inchikey
            # First, get an existing molecule
            response = self.supabase.table("molecules").select("*").limit(1).execute()
            if not response.data:
                self.add_test_result("DB-4.1", "Data Validation - Duplicate InChIKey", "Skipped", "No molecules found to test with")
                return
            
            existing_molecule = response.data[0]
            
            # Try to insert a molecule with the same inchikey
            try:
                new_molecule = {
                    "id": str(uuid.uuid4()),
                    "name": "Test Molecule",
                    "inchikey": existing_molecule["inchikey"],
                    "created_by": existing_molecule["created_by"]
                }
                self.supabase.table("molecules").insert(new_molecule).execute()
                self.add_test_result("DB-4.1", "Data Validation - Duplicate InChIKey", "Failed", "Able to insert molecule with duplicate inchikey")
            except Exception as e:
                if "duplicate key value violates unique constraint" in str(e):
                    self.add_test_result("DB-4.1", "Data Validation - Duplicate InChIKey", "Passed", "Properly rejected duplicate inchikey")
                else:
                    self.add_test_result("DB-4.1", "Data Validation - Duplicate InChIKey", "Failed", f"Unexpected error: {str(e)}")
            
            # Test that we can't insert a prediction with both molecule_id and mixture_id
            try:
                # Get a molecule and mixture
                molecule_response = self.supabase.table("molecules").select("id").limit(1).execute()
                mixture_response = self.supabase.table("mixtures").select("id").limit(1).execute()
                property_type_response = self.supabase.table("property_types").select("id").limit(1).execute()
                calculation_method_response = self.supabase.table("calculation_methods").select("id").limit(1).execute()
                
                if not molecule_response.data or not mixture_response.data or not property_type_response.data or not calculation_method_response.data:
                    self.add_test_result("DB-4.2", "Data Validation - Prediction Constraints", "Skipped", "Missing required data for test")
                    return
                
                molecule_id = molecule_response.data[0]["id"]
                mixture_id = mixture_response.data[0]["id"]
                property_type_id = property_type_response.data[0]["id"]
                calculation_method_id = calculation_method_response.data[0]["id"]
                
                # Try to insert a prediction with both molecule_id and mixture_id
                new_prediction = {
                    "id": str(uuid.uuid4()),
                    "molecule_id": molecule_id,
                    "mixture_id": mixture_id,
                    "property_type_id": property_type_id,
                    "calculation_method_id": calculation_method_id,
                    "numeric_value": 42.0,
                    "created_by": existing_molecule["created_by"]
                }
                self.supabase.table("predictions").insert(new_prediction).execute()
                self.add_test_result("DB-4.2", "Data Validation - Prediction Constraints", "Failed", "Able to insert prediction with both molecule_id and mixture_id")
            except Exception as e:
                if "check constraint" in str(e) or "violates check constraint" in str(e):
                    self.add_test_result("DB-4.2", "Data Validation - Prediction Constraints", "Passed", "Properly rejected prediction with both molecule_id and mixture_id")
                else:
                    self.add_test_result("DB-4.2", "Data Validation - Prediction Constraints", "Failed", f"Unexpected error: {str(e)}")
        except Exception as e:
            self.add_test_result("DB-4", "Data Validation", "Failed", f"Error testing data validation: {str(e)}")

    def test_rls_enablement(self) -> None:
        """Test that RLS is enabled on all tables."""
        try:
            query = """
            SELECT
                tablename,
                rowsecurity
            FROM
                pg_tables
            WHERE
                schemaname = 'public'
            """
            response = self.supabase.rpc("exec_sql", {"query": query}).execute()
            
            tables_without_rls = []
            for row in response.data:
                if row["tablename"] in self.expected_tables and not row["rowsecurity"]:
                    tables_without_rls.append(row["tablename"])
            
            if not tables_without_rls:
                self.add_test_result("RLS-1", "RLS Enablement", "Passed", "RLS is enabled on all tables")
            else:
                self.add_test_result("RLS-1", "RLS Enablement", "Failed", f"RLS is not enabled on tables: {', '.join(tables_without_rls)}")
        except Exception as e:
            self.add_test_result("RLS-1", "RLS Enablement", "Failed", f"Error checking RLS enablement: {str(e)}")

    def save_results(self, filename: str) -> None:
        """Save the test results to a file."""
        try:
            with open(filename, 'w') as f:
                json.dump(self.test_results, f, indent=2)
            logger.info(f"Test results saved to {filename}")
        except Exception as e:
            logger.error(f"Failed to save test results: {str(e)}")

def main():
    """Main function."""
    logger.info("Starting database schema and data integrity tests")
    
    # Create and run the tests
    test = DatabaseSchemaTest()
    results = test.run_tests()
    
    # Save the results
    test.save_results("database_schema_test_results.json")
    
    # Print summary
    logger.info(f"Test Status: {results['status']}")
    logger.info(f"Total Tests: {results['total_tests']}")
    logger.info(f"Passed Tests: {results['passed_tests']}")
    logger.info(f"Failed Tests: {results['failed_tests']}")
    logger.info(f"Skipped Tests: {results['skipped_tests']}")
    
    # Exit with appropriate status code
    if results["status"] == "Passed":
        logger.info("All tests passed!")
        sys.exit(0)
    else:
        logger.error("Some tests failed. See log for details.")
        sys.exit(1)

if __name__ == "__main__":
    main()
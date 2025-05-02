#!/usr/bin/env python3
"""
CryoProtect v2 - Database Migration Verification Suite

This script performs comprehensive verification of the database migration from singular to plural table names.
It tests database structure, data integrity, RLS policies, and API functionality.

Author: CryoProtect Team
Date: 2025-04-18
"""

import os
import sys
import json
import logging
import argparse
import time
import requests
from datetime import datetime
import psycopg2
from psycopg2 import sql
from psycopg2.extras import RealDictCursor
from dotenv import load_dotenv
import colorama
from colorama import Fore, Style
from tabulate import tabulate

# Initialize colorama
colorama.init()

# Configure logging
log_filename = f"migration_verification_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(log_filename),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

# Database connection parameters
DB_URL = os.getenv("DATABASE_URL") or os.getenv("SUPABASE_DB_URL")
if not DB_URL:
    logger.error("Database URL not found in environment variables")
    sys.exit(1)

# API base URL
API_BASE_URL = os.getenv("API_BASE_URL") or "http://localhost:5000/api/v1"

# Table mappings
TABLE_MAPPINGS = [
    {"singular": "molecule", "plural": "molecules"},
    {"singular": "mixture", "plural": "mixtures"},
    {"singular": "experiment", "plural": "experiments"},
    {"singular": "prediction", "plural": "predictions"},
    {"singular": "project", "plural": "projects"}
]

# Foreign key relationships
FK_RELATIONSHIPS = [
    {
        "table": "mixture_components",
        "columns": [
            {"name": "mixture_id", "references": "mixture", "new_references": "mixtures"},
            {"name": "molecule_id", "references": "molecule", "new_references": "molecules"}
        ]
    },
    {
        "table": "molecular_properties",
        "columns": [
            {"name": "molecule_id", "references": "molecule", "new_references": "molecules"}
        ]
    },
    {
        "table": "predictions",
        "columns": [
            {"name": "molecule_id", "references": "molecule", "new_references": "molecules"},
            {"name": "mixture_id", "references": "mixture", "new_references": "mixtures"}
        ]
    },
    {
        "table": "experiments",
        "columns": [
            {"name": "molecule_id", "references": "molecule", "new_references": "molecules"},
            {"name": "mixture_id", "references": "mixture", "new_references": "mixtures"}
        ]
    },
    {
        "table": "project_membership",
        "columns": [
            {"name": "project_id", "references": "project", "new_references": "projects"}
        ]
    }
]

# Views that should have been updated
VIEWS_TO_CHECK = [
    "experiment_with_results",
    "mixture_with_components",
    "mixtures_with_components",
    "molecule_with_properties",
    "molecules_with_properties"
]

# API endpoints to test
API_ENDPOINTS = [
    # Singular endpoints (should still work with deprecation notice)
    {"path": "/molecule", "method": "GET", "expected_status": 200},
    {"path": "/molecule/{id}", "method": "GET", "expected_status": 200},
    {"path": "/mixture", "method": "GET", "expected_status": 200},
    {"path": "/mixture/{id}", "method": "GET", "expected_status": 200},
    {"path": "/experiment", "method": "GET", "expected_status": 200},
    {"path": "/experiment/{id}", "method": "GET", "expected_status": 200},
    {"path": "/prediction", "method": "GET", "expected_status": 200},
    {"path": "/prediction/{id}", "method": "GET", "expected_status": 200},
    
    # Plural endpoints (new standard)
    {"path": "/molecules", "method": "GET", "expected_status": 200},
    {"path": "/molecules/{id}", "method": "GET", "expected_status": 200},
    {"path": "/mixtures", "method": "GET", "expected_status": 200},
    {"path": "/mixtures/{id}", "method": "GET", "expected_status": 200},
    {"path": "/experiments", "method": "GET", "expected_status": 200},
    {"path": "/experiments/{id}", "method": "GET", "expected_status": 200},
    {"path": "/predictions", "method": "GET", "expected_status": 200},
    {"path": "/predictions/{id}", "method": "GET", "expected_status": 200}
]

def connect_to_db():
    """
    Establish a connection to the database.
    
    Returns:
        tuple: (connection, cursor) or (None, None) if connection fails
    """
    try:
        conn = psycopg2.connect(DB_URL)
        conn.autocommit = True  # Auto-commit for verification queries
        cursor = conn.cursor(cursor_factory=RealDictCursor)
        logger.info("Successfully connected to the database")
        return conn, cursor
    except Exception as e:
        logger.error(f"Error connecting to the database: {str(e)}")
        return None, None

def print_section_header(title):
    """Print a formatted section header."""
    print("\n" + "=" * 80)
    print(f"{Fore.CYAN}{Style.BRIGHT}{title}{Style.RESET_ALL}")
    print("=" * 80)

def print_test_result(test_name, result, details=None):
    """Print a formatted test result."""
    if result:
        status = f"{Fore.GREEN}PASS{Style.RESET_ALL}"
    else:
        status = f"{Fore.RED}FAIL{Style.RESET_ALL}"
    
    print(f"{status} - {test_name}")
    
    if details and not result:
        print(f"       {Fore.YELLOW}Details: {details}{Style.RESET_ALL}")

def verify_table_existence(cursor, schema_name='public'):
    """
    Verify that all plural tables exist.
    
    Args:
        cursor: Database cursor
        schema_name (str): Database schema name
        
    Returns:
        dict: Verification results
    """
    print_section_header("Verifying Table Existence")
    
    results = {
        "passed": True,
        "tables_checked": [],
        "missing_tables": []
    }
    
    for mapping in TABLE_MAPPINGS:
        singular = mapping["singular"]
        plural = mapping["plural"]
        
        # Check if plural table exists
        cursor.execute("""
            SELECT EXISTS (
                SELECT FROM information_schema.tables 
                WHERE table_schema = %s AND table_name = %s
            ) as exists
        """, (schema_name, plural))
        
        plural_exists = cursor.fetchone()['exists']
        
        # Check if singular table still exists (it should)
        cursor.execute("""
            SELECT EXISTS (
                SELECT FROM information_schema.tables 
                WHERE table_schema = %s AND table_name = %s
            ) as exists
        """, (schema_name, singular))
        
        singular_exists = cursor.fetchone()['exists']
        
        results["tables_checked"].append({
            "singular": singular,
            "plural": plural,
            "singular_exists": singular_exists,
            "plural_exists": plural_exists
        })
        
        if not plural_exists:
            results["passed"] = False
            results["missing_tables"].append(plural)
            print_test_result(f"Table {plural} exists", False, "Table not found")
        else:
            print_test_result(f"Table {plural} exists", True)
    
    return results

def verify_table_structure(cursor, schema_name='public'):
    """
    Verify that plural tables have the correct structure.
    
    Args:
        cursor: Database cursor
        schema_name (str): Database schema name
        
    Returns:
        dict: Verification results
    """
    print_section_header("Verifying Table Structure")
    
    results = {
        "passed": True,
        "tables_checked": [],
        "structure_issues": []
    }
    
    for mapping in TABLE_MAPPINGS:
        singular = mapping["singular"]
        plural = mapping["plural"]
        
        # Check if both tables exist
        cursor.execute("""
            SELECT EXISTS (
                SELECT FROM information_schema.tables 
                WHERE table_schema = %s AND table_name = %s
            ) as exists
        """, (schema_name, singular))
        
        singular_exists = cursor.fetchone()['exists']
        
        cursor.execute("""
            SELECT EXISTS (
                SELECT FROM information_schema.tables 
                WHERE table_schema = %s AND table_name = %s
            ) as exists
        """, (schema_name, plural))
        
        plural_exists = cursor.fetchone()['exists']
        
        if not singular_exists or not plural_exists:
            if not singular_exists:
                logger.warning(f"Singular table {singular} does not exist, skipping structure check")
            if not plural_exists:
                logger.warning(f"Plural table {plural} does not exist, skipping structure check")
            continue
        
        # Get column information for both tables
        cursor.execute("""
            SELECT column_name, data_type, is_nullable, column_default
            FROM information_schema.columns
            WHERE table_schema = %s AND table_name = %s
            ORDER BY ordinal_position
        """, (schema_name, singular))
        
        singular_columns = cursor.fetchall()
        
        cursor.execute("""
            SELECT column_name, data_type, is_nullable, column_default
            FROM information_schema.columns
            WHERE table_schema = %s AND table_name = %s
            ORDER BY ordinal_position
        """, (schema_name, plural))
        
        plural_columns = cursor.fetchall()
        
        # Compare column count
        if len(singular_columns) != len(plural_columns):
            results["passed"] = False
            issue = {
                "table": plural,
                "issue": "Column count mismatch",
                "singular_count": len(singular_columns),
                "plural_count": len(plural_columns)
            }
            results["structure_issues"].append(issue)
            print_test_result(f"Table {plural} column count matches", False, 
                             f"Singular: {len(singular_columns)}, Plural: {len(plural_columns)}")
            continue
        else:
            print_test_result(f"Table {plural} column count matches", True)
        
        # Compare column definitions
        column_issues = []
        for i, singular_col in enumerate(singular_columns):
            plural_col = plural_columns[i]
            
            if (singular_col['column_name'] != plural_col['column_name'] or
                singular_col['data_type'] != plural_col['data_type'] or
                singular_col['is_nullable'] != plural_col['is_nullable']):
                
                column_issues.append({
                    "singular_column": singular_col['column_name'],
                    "plural_column": plural_col['column_name'],
                    "mismatch": "Column definition differs"
                })
        
        if column_issues:
            results["passed"] = False
            issue = {
                "table": plural,
                "issue": "Column definition mismatch",
                "details": column_issues
            }
            results["structure_issues"].append(issue)
            print_test_result(f"Table {plural} column definitions match", False, 
                             f"{len(column_issues)} columns have definition mismatches")
        else:
            print_test_result(f"Table {plural} column definitions match", True)
        
        # Check primary key
        cursor.execute("""
            SELECT c.column_name
            FROM information_schema.table_constraints tc
            JOIN information_schema.constraint_column_usage AS ccu USING (constraint_schema, constraint_name)
            JOIN information_schema.columns AS c ON c.table_schema = tc.constraint_schema
                AND tc.table_name = c.table_name AND ccu.column_name = c.column_name
            WHERE tc.constraint_type = 'PRIMARY KEY' AND tc.table_schema = %s AND tc.table_name = %s
        """, (schema_name, singular))
        
        singular_pk = [row['column_name'] for row in cursor.fetchall()]
        
        cursor.execute("""
            SELECT c.column_name
            FROM information_schema.table_constraints tc
            JOIN information_schema.constraint_column_usage AS ccu USING (constraint_schema, constraint_name)
            JOIN information_schema.columns AS c ON c.table_schema = tc.constraint_schema
                AND tc.table_name = c.table_name AND ccu.column_name = c.column_name
            WHERE tc.constraint_type = 'PRIMARY KEY' AND tc.table_schema = %s AND tc.table_name = %s
        """, (schema_name, plural))
        
        plural_pk = [row['column_name'] for row in cursor.fetchall()]
        
        if set(singular_pk) != set(plural_pk):
            results["passed"] = False
            issue = {
                "table": plural,
                "issue": "Primary key mismatch",
                "singular_pk": singular_pk,
                "plural_pk": plural_pk
            }
            results["structure_issues"].append(issue)
            print_test_result(f"Table {plural} primary key matches", False, 
                             f"Singular: {singular_pk}, Plural: {plural_pk}")
        else:
            print_test_result(f"Table {plural} primary key matches", True)
        
        results["tables_checked"].append({
            "singular": singular,
            "plural": plural,
            "columns_match": len(column_issues) == 0,
            "pk_matches": set(singular_pk) == set(plural_pk)
        })
    
    return results

def verify_foreign_keys(cursor, schema_name='public'):
    """
    Verify that foreign key constraints are correctly defined.
    
    Args:
        cursor: Database cursor
        schema_name (str): Database schema name
        
    Returns:
        dict: Verification results
    """
    print_section_header("Verifying Foreign Key Constraints")
    
    results = {
        "passed": True,
        "fk_checked": [],
        "fk_issues": []
    }
    
    for rel in FK_RELATIONSHIPS:
        table = rel["table"]
        
        for column in rel["columns"]:
            col_name = column["name"]
            new_references = column["new_references"]
            
            # Check if the foreign key exists and references the plural table
            cursor.execute("""
                SELECT
                    tc.constraint_name,
                    ccu.table_name AS foreign_table_name
                FROM information_schema.table_constraints AS tc
                JOIN information_schema.key_column_usage AS kcu
                    ON tc.constraint_name = kcu.constraint_name
                    AND tc.table_schema = kcu.table_schema
                JOIN information_schema.constraint_column_usage AS ccu
                    ON ccu.constraint_name = tc.constraint_name
                    AND ccu.table_schema = tc.table_schema
                WHERE tc.constraint_type = 'FOREIGN KEY'
                    AND tc.table_schema = %s
                    AND tc.table_name = %s
                    AND kcu.column_name = %s
            """, (schema_name, table, col_name))
            
            fk_info = cursor.fetchone()
            
            if not fk_info:
                results["passed"] = False
                issue = {
                    "table": table,
                    "column": col_name,
                    "issue": "Foreign key constraint not found"
                }
                results["fk_issues"].append(issue)
                print_test_result(f"Foreign key {table}.{col_name} exists", False, "Constraint not found")
                continue
            
            if fk_info['foreign_table_name'] != new_references:
                results["passed"] = False
                issue = {
                    "table": table,
                    "column": col_name,
                    "issue": "Foreign key references wrong table",
                    "expected": new_references,
                    "actual": fk_info['foreign_table_name']
                }
                results["fk_issues"].append(issue)
                print_test_result(f"Foreign key {table}.{col_name} references {new_references}", False, 
                                 f"References {fk_info['foreign_table_name']} instead")
            else:
                print_test_result(f"Foreign key {table}.{col_name} references {new_references}", True)
            
            results["fk_checked"].append({
                "table": table,
                "column": col_name,
                "expected_reference": new_references,
                "actual_reference": fk_info['foreign_table_name'] if fk_info else None,
                "constraint_name": fk_info['constraint_name'] if fk_info else None
            })
    
    return results

def verify_data_integrity(cursor, schema_name='public'):
    """
    Verify that all data was correctly migrated from singular to plural tables.
    
    Args:
        cursor: Database cursor
        schema_name (str): Database schema name
        
    Returns:
        dict: Verification results
    """
    print_section_header("Verifying Data Integrity")
    
    results = {
        "passed": True,
        "tables_checked": [],
        "data_issues": []
    }
    
    for mapping in TABLE_MAPPINGS:
        singular = mapping["singular"]
        plural = mapping["plural"]
        
        # Check if both tables exist
        cursor.execute("""
            SELECT EXISTS (
                SELECT FROM information_schema.tables 
                WHERE table_schema = %s AND table_name = %s
            ) as exists
        """, (schema_name, singular))
        
        singular_exists = cursor.fetchone()['exists']
        
        cursor.execute("""
            SELECT EXISTS (
                SELECT FROM information_schema.tables 
                WHERE table_schema = %s AND table_name = %s
            ) as exists
        """, (schema_name, plural))
        
        plural_exists = cursor.fetchone()['exists']
        
        if not singular_exists or not plural_exists:
            if not singular_exists:
                logger.warning(f"Singular table {singular} does not exist, skipping data integrity check")
            if not plural_exists:
                logger.warning(f"Plural table {plural} does not exist, skipping data integrity check")
            continue
        
        # Count rows in both tables
        cursor.execute(
            sql.SQL("SELECT COUNT(*) as count FROM {}").format(
                sql.Identifier(singular)
            )
        )
        singular_count = cursor.fetchone()['count']
        
        cursor.execute(
            sql.SQL("SELECT COUNT(*) as count FROM {}").format(
                sql.Identifier(plural)
            )
        )
        plural_count = cursor.fetchone()['count']
        
        # Compare row counts
        if singular_count != plural_count:
            results["passed"] = False
            issue = {
                "table": plural,
                "issue": "Row count mismatch",
                "singular_count": singular_count,
                "plural_count": plural_count
            }
            results["data_issues"].append(issue)
            print_test_result(f"Table {plural} row count matches", False, 
                             f"Singular: {singular_count}, Plural: {plural_count}")
        else:
            print_test_result(f"Table {plural} row count matches", True)
        
        # Get primary key column
        cursor.execute("""
            SELECT c.column_name
            FROM information_schema.table_constraints tc
            JOIN information_schema.constraint_column_usage AS ccu USING (constraint_schema, constraint_name)
            JOIN information_schema.columns AS c ON c.table_schema = tc.constraint_schema
                AND tc.table_name = c.table_name AND ccu.column_name = c.column_name
            WHERE tc.constraint_type = 'PRIMARY KEY' AND tc.table_schema = %s AND tc.table_name = %s
        """, (schema_name, singular))
        
        pk_columns = [row['column_name'] for row in cursor.fetchall()]
        
        if not pk_columns:
            logger.warning(f"No primary key found for {singular}, using all columns for data integrity check")
            # Get all columns
            cursor.execute("""
                SELECT column_name
                FROM information_schema.columns
                WHERE table_schema = %s AND table_name = %s
            """, (schema_name, singular))
            columns = [row['column_name'] for row in cursor.fetchall()]
        else:
            columns = pk_columns
        
        # Sample data from singular table
        column_identifiers = [sql.Identifier(col) for col in columns]
        
        cursor.execute(
            sql.SQL("SELECT {} FROM {} LIMIT 100").format(
                sql.SQL(", ").join(column_identifiers),
                sql.Identifier(singular)
            )
        )
        sample_rows = cursor.fetchall()
        
        # Check if each sample row exists in plural table
        missing_rows = []
        for row in sample_rows:
            where_conditions = []
            for col in columns:
                where_conditions.append(
                    sql.SQL("{} = %s").format(sql.Identifier(col))
                )
            
            cursor.execute(
                sql.SQL("SELECT COUNT(*) as count FROM {} WHERE {}").format(
                    sql.Identifier(plural),
                    sql.SQL(" AND ").join(where_conditions)
                ),
                [row[col] for col in columns]
            )
            
            match_count = cursor.fetchone()['count']
            if match_count != 1:
                missing_rows.append({
                    "primary_key": {col: str(row[col]) for col in columns},
                    "match_count": match_count
                })
        
        if missing_rows:
            results["passed"] = False
            issue = {
                "table": plural,
                "issue": "Missing or duplicate rows",
                "details": missing_rows
            }
            results["data_issues"].append(issue)
            print_test_result(f"Table {plural} data integrity check", False, 
                             f"{len(missing_rows)} rows have issues")
        else:
            print_test_result(f"Table {plural} data integrity check", True)
        
        results["tables_checked"].append({
            "singular": singular,
            "plural": plural,
            "singular_row_count": singular_count,
            "plural_row_count": plural_count,
            "sample_size": len(sample_rows),
            "missing_rows": missing_rows
        })
    
    return results

def verify_rls_policies(cursor, schema_name='public'):
    """
    Verify that RLS is enabled on all plural tables and policies are correctly defined.
    
    Args:
        cursor: Database cursor
        schema_name (str): Database schema name
        
    Returns:
        dict: Verification results
    """
    print_section_header("Verifying RLS Policies")
    
    results = {
        "passed": True,
        "tables_checked": [],
        "rls_issues": []
    }
    
    for mapping in TABLE_MAPPINGS:
        singular = mapping["singular"]
        plural = mapping["plural"]
        
        # Check if RLS is enabled on plural table
        cursor.execute("""
            SELECT relrowsecurity
            FROM pg_class c
            JOIN pg_namespace n ON n.oid = c.relnamespace
            WHERE n.nspname = %s AND c.relkind = 'r' AND c.relname = %s
        """, (schema_name, plural))
        
        result = cursor.fetchone()
        if not result:
            logger.warning(f"Table {plural} not found, skipping RLS check")
            continue
        
        rls_enabled = result['relrowsecurity']
        
        if not rls_enabled:
            results["passed"] = False
            issue = {
                "table": plural,
                "issue": "RLS not enabled"
            }
            results["rls_issues"].append(issue)
            print_test_result(f"Table {plural} has RLS enabled", False)
        else:
            print_test_result(f"Table {plural} has RLS enabled", True)
        
        # Get policies for singular table
        cursor.execute("""
            SELECT 
                policyname, 
                permissive,
                roles,
                cmd, 
                qual, 
                with_check
            FROM pg_policies
            WHERE schemaname = %s AND tablename = %s
        """, (schema_name, singular))
        
        singular_policies = cursor.fetchall()
        
        # Get policies for plural table
        cursor.execute("""
            SELECT 
                policyname, 
                permissive,
                roles,
                cmd, 
                qual, 
                with_check
            FROM pg_policies
            WHERE schemaname = %s AND tablename = %s
        """, (schema_name, plural))
        
        plural_policies = cursor.fetchall()
        
        # Compare policy count
        if len(singular_policies) != len(plural_policies):
            results["passed"] = False
            issue = {
                "table": plural,
                "issue": "Policy count mismatch",
                "singular_count": len(singular_policies),
                "plural_count": len(plural_policies)
            }
            results["rls_issues"].append(issue)
            print_test_result(f"Table {plural} policy count matches", False, 
                             f"Singular: {len(singular_policies)}, Plural: {len(plural_policies)}")
        else:
            print_test_result(f"Table {plural} policy count matches", True)
        
        # Check that each policy on singular table has a corresponding policy on plural table
        missing_policies = []
        for singular_policy in singular_policies:
            policy_name = singular_policy['policyname']
            
            # Find corresponding policy on plural table
            matching_policies = [p for p in plural_policies if p['policyname'] == policy_name]
            
            if not matching_policies:
                missing_policies.append({
                    "policy_name": policy_name,
                    "cmd": singular_policy['cmd'],
                    "roles": singular_policy['roles']
                })
                continue
            
            # Check if policy expressions are correctly updated
            plural_policy = matching_policies[0]
            
            # Check if singular table names in expressions were replaced with plural
            if singular_policy['qual']:
                for mapping_item in TABLE_MAPPINGS:
                    sing = mapping_item["singular"]
                    plur = mapping_item["plural"]
                    
                    if sing in singular_policy['qual'] and sing in plural_policy['qual']:
                        # Bad, the singular name was not replaced
                        results["passed"] = False
                        issue = {
                            "table": plural,
                            "policy": policy_name,
                            "issue": "Policy expression not updated",
                            "singular_expression": singular_policy['qual'],
                            "plural_expression": plural_policy['qual']
                        }
                        results["rls_issues"].append(issue)
                        print_test_result(f"Policy {policy_name} on {plural} has updated expressions", False, 
                                         f"Expression still contains singular table name '{sing}'")
                        break
                else:
                    print_test_result(f"Policy {policy_name} on {plural} has updated expressions", True)
        
        if missing_policies:
            results["passed"] = False
            issue = {
                "table": plural,
                "issue": "Missing policies",
                "details": missing_policies
            }
            results["rls_issues"].append(issue)
            print_test_result(f"Table {plural} has all required policies", False, 
                             f"Missing {len(missing_policies)} policies")
        else:
            print_test_result(f"Table {plural} has all required policies", True)
        
        results["tables_checked"].append({
            "singular": singular,
            "plural": plural,
            "rls_enabled": rls_enabled,
            "singular_policy_count": len(singular_policies),
            "plural_policy_count": len(plural_policies),
            "missing_policies": missing_policies
        })
    
    return results

def verify_views(cursor, schema_name='public'):
    """
    Verify that views have been updated to reference plural tables.
    
    Args:
        cursor: Database cursor
        schema_name (str): Database schema name
        
    Returns:
        dict: Verification results
    """
    print_section_header("Verifying Views")
    
    results = {
        "passed": True,
        "views_checked": [],
        "view_issues": []
    }
    
    for view_name in VIEWS_TO_CHECK:
        # Check if view exists
        cursor.execute("""
            SELECT EXISTS (
                SELECT FROM information_schema.views 
                WHERE table_schema = %s AND table_name = %s
            ) as exists
        """, (schema_name, view_name))
        
        view_exists = cursor.fetchone()['exists']
        
        if not view_exists:
            logger.warning(f"View {view_name} does not exist, skipping check")
            continue
        
        # Get view definition
        cursor.execute("""
            SELECT pg_get_viewdef(c.oid, true) as view_definition
            FROM pg_class c
            JOIN pg_namespace n ON n.oid = c.relnamespace
            WHERE c.relkind = 'v'
            AND n.nspname = %s
            AND c.relname = %s
        """, (schema_name, view_name))
        
        result = cursor.fetchone()
        if not result:
            logger.warning(f"Could not get definition for view {view_name}, skipping check")
            continue
        
        view_definition = result['view_definition']
        
        # Check if view definition references singular tables
        references_singular = False
        for mapping in TABLE_MAPPINGS:
            singular = mapping["singular"]
            if f" {singular} " in view_definition or f" {singular}." in view_definition:
                references_singular = True
                results["passed"] = False
                issue = {
                    "view": view_name,
                    "issue": "References singular table",
                    "singular_table": singular,
                    "view_definition": view_definition
                }
                results["view_issues"].append(issue)
                print_test_result(f"View {view_name} references plural tables", False, 
                                 f"References singular table '{singular}'")
                break
        
        if not references_singular:
            print_test_result(f"View {view_name} references plural tables", True)
        
        results["views_checked"].append({
            "view": view_name,
            "exists": view_exists,
            "references_singular": references_singular
        })
    
    return results

def verify_api_compatibility(api_base_url=API_BASE_URL):
    """
    Verify that the API compatibility layer is working correctly.
    
    Args:
        api_base_url (str): Base URL for the API
        
    Returns:
        dict: Verification results
    """
    print_section_header("Verifying API Compatibility")
    
    results = {
        "passed": True,
        "endpoints_checked": [],
        "endpoint_issues": []
    }
    
    # Get a sample ID for each resource type
    sample_ids = {}
    conn, cursor = connect_to_db()
    if not conn or not cursor:
        logger.error("Could not connect to database to get sample IDs")
        return {
            "passed": False,
            "error": "Database connection failed"
        }
    
    for mapping in TABLE_MAPPINGS:
        plural = mapping["plural"]
        try:
            cursor.execute(
                sql.SQL("SELECT id FROM {} LIMIT 1").format(
                    sql.Identifier(plural)
                )
            )
            result = cursor.fetchone()
            if result:
                sample_ids[plural] = result['id']
        except Exception as e:
            logger.warning(f"Could not get sample ID for {plural}: {str(e)}")
    
    conn.close()
    
    # Test each endpoint
    for endpoint in API_ENDPOINTS:
        path = endpoint["path"]
        method = endpoint["method"]
        expected_status = endpoint["expected_status"]
        
        # Replace {id} placeholders with actual IDs
        if "{id}" in path:
            resource_type = path.split("/")[1]  # e.g., "molecule" or "molecules"
            
            # Find the corresponding plural form
            plural_form = None
            for mapping in TABLE_MAPPINGS:
                if mapping["singular"] == resource_type or mapping["plural"] == resource_type:
                    plural_form = mapping["plural"]
                    break
            
            if not plural_form or plural_form not in sample_ids:
                logger.warning(f"No sample ID available for {resource_type}, skipping endpoint {path}")
                continue
            
            path = path.replace("{id}", str(sample_ids[plural_form]))
        
        # Make the request
        url = f"{api_base_url}{path}"
        try:
            if method == "GET":
                response = requests.get(url)
            elif method == "POST":
                response = requests.post(url, json={})
            elif method == "PUT":
                response = requests.put(url, json={})
            elif method == "DELETE":
                response = requests.delete(url)
            else:
                logger.warning(f"Unsupported method {method} for endpoint {path}")
                continue
            
            # Check status code
            status_matches = response.status_code == expected_status
            
            # For singular endpoints, check for deprecation notice
            is_singular = any(mapping["singular"] in path.split("/")[1] for mapping in TABLE_MAPPINGS)
            has_deprecation_notice = False
            
            if is_singular and response.status_code == expected_status:
                try:
                    response_json = response.json()
                    has_deprecation_notice = "_deprecation" in response_json
                except:
                    has_deprecation_notice = False
            
            endpoint_result = {
                "endpoint": path,
                "method": method,
                "expected_status": expected_status,
                "actual_status": response.status_code,
                "status_matches": status_matches
            }
            
            if is_singular:
                endpoint_result["is_singular"] = True
                endpoint_result["has_deprecation_notice"] = has_deprecation_notice
                
                if not has_deprecation_notice:
                    results["passed"] = False
                    issue = {
                        "endpoint": path,
                        "issue": "Missing deprecation notice",
                        "method": method,
                        "status": response.status_code
                    }
                    results["endpoint_issues"].append(issue)
                    print_test_result(f"Endpoint {method} {path} has deprecation notice", False)
                else:
                    print_test_result(f"Endpoint {method} {path} has deprecation notice", True)
            
            if not status_matches:
                results["passed"] = False
                issue = {
                    "endpoint": path,
                    "issue": "Status code mismatch",
                    "method": method,
                    "expected_status": expected_status,
                    "actual_status": response.status_code
                }
                results["endpoint_issues"].append(issue)
                print_test_result(f"Endpoint {method} {path} returns status {expected_status}", False, 
                                 f"Returned status {response.status_code}")
            else:
                print_test_result(f"Endpoint {method} {path} returns status {expected_status}", True)
            
            results["endpoints_checked"].append(endpoint_result)
            
        except Exception as e:
            results["passed"] = False
            issue = {
                "endpoint": path,
                "issue": "Request failed",
                "method": method,
                "error": str(e)
            }
            results["endpoint_issues"].append(issue)
            print_test_result(f"Endpoint {method} {path} is accessible", False, str(e))
    
    return results

def run_verification_tests(schema_name='public', api_base_url=API_BASE_URL, categories=None):
    """
    Run all verification tests.
    
    Args:
        schema_name (str): Database schema name
        api_base_url (str): Base URL for the API
        categories (list): List of test categories to run, or None for all
        
    Returns:
        dict: Verification results
    """
    conn, cursor = connect_to_db()
    if not conn or not cursor:
        return {
            "passed": False,
            "error": "Database connection failed"
        }
    
    results = {
        "timestamp": datetime.now().isoformat(),
        "passed": True,
        "categories": {}
    }
    
    try:
        # Define test categories
        test_categories = {
            "table_existence": lambda: verify_table_existence(cursor, schema_name),
            "table_structure": lambda: verify_table_structure(cursor, schema_name),
            "foreign_keys": lambda: verify_foreign_keys(cursor, schema_name),
            "data_integrity": lambda: verify_data_integrity(cursor, schema_name),
            "rls_policies": lambda: verify_rls_policies(cursor, schema_name),
            "views": lambda: verify_views(cursor, schema_name),
            "api_compatibility": lambda: verify_api_compatibility(api_base_url)
        }
        
        # Run selected tests or all tests
        if categories:
            selected_categories = {k: v for k, v in test_categories.items() if k in categories}
        else:
            selected_categories = test_categories
        
        # Run tests
        for category, test_func in selected_categories.items():
            print(f"\nRunning {category} tests...")
            category_results = test_func()
            results["categories"][category] = category_results
            
            if not category_results.get("passed", True):
                results["passed"] = False
        
    finally:
        conn.close()
    
    return results

def generate_report(results, output_file=None):
    """
    Generate a detailed report of verification results.
    
    Args:
        results (dict): Verification results
        output_file (str): Path to output file, or None for stdout
        
    Returns:
        str: Report content
    """
    report = []
    
    # Add header
    report.append("# CryoProtect v2 Database Migration Verification Report")
    report.append("")
    report.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    report.append("")
    
    # Add overall status
    if results.get("passed", False):
        report.append("## Overall Status: PASSED")
        report.append("")
        report.append("All verification tests passed successfully.")
    else:
        report.append("## Overall Status: FAILED")
        report.append("")
        report.append("Some verification tests failed. See details below.")
    
    report.append("")
    
    # Add summary table
    report.append("## Summary")
    report.append("")
    
    summary_rows = []
    for category, category_results in results.get("categories", {}).items():
        status = "PASS" if category_results.get("passed", True) else "FAIL"
        summary_rows.append([category, status])
    
    summary_table = tabulate(summary_rows, headers=["Category", "Status"], tablefmt="pipe")
    report.append(summary_table)
    report.append("")
    
    # Add detailed results for each category
    report.append("## Detailed Results")
    report.append("")
    
    for category, category_results in results.get("categories", {}).items():
        report.append(f"### {category.replace('_', ' ').title()}")
        report.append("")
        
        if category == "table_existence":
            # Table existence details
            tables_checked = category_results.get("tables_checked", [])
            missing_tables = category_results.get("missing_tables", [])
            
            report.append(f"Tables checked: {len(tables_checked)}")
            report.append(f"Missing tables: {len(missing_tables)}")
            report.append("")
            
            if missing_tables:
                report.append("#### Missing Tables")
                report.append("")
                for table in missing_tables:
                    report.append(f"- {table}")
                report.append("")
        
        elif category == "table_structure":
            # Table structure details
            tables_checked = category_results.get("tables_checked", [])
            structure_issues = category_results.get("structure_issues", [])
            
            report.append(f"Tables checked: {len(tables_checked)}")
            report.append(f"Tables with structure issues: {len(structure_issues)}")
            report.append("")
            
            if structure_issues:
                report.append("#### Structure Issues")
                report.append("")
                for issue in structure_issues:
                    report.append(f"- Table: {issue.get('table')}")
                    report.append(f"  - Issue: {issue.get('issue')}")
                    if 'details' in issue:
                        report.append(f"  - Details: {issue.get('details')}")
                    report.append("")
        
        elif category == "foreign_keys":
            # Foreign key details
            fk_checked = category_results.get("fk_checked", [])
            fk_issues = category_results.get("fk_issues", [])
            
            report.append(f"Foreign keys checked: {len(fk_checked)}")
            report.append(f"Foreign keys with issues: {len(fk_issues)}")
            report.append("")
            
            if fk_issues:
                report.append("#### Foreign Key Issues")
                report.append("")
                for issue in fk_issues:
                    report.append(f"- Table: {issue.get('table')}, Column: {issue.get('column')}")
                    report.append(f"  - Issue: {issue.get('issue')}")
                    if 'expected' in issue and 'actual' in issue:
                        report.append(f"  - Expected: {issue.get('expected')}, Actual: {issue.get('actual')}")
                    report.append("")
        
        elif category == "data_integrity":
            # Data integrity details
            tables_checked = category_results.get("tables_checked", [])
            data_issues = category_results.get("data_issues", [])
            
            report.append(f"Tables checked: {len(tables_checked)}")
            report.append(f"Tables with data issues: {len(data_issues)}")
            report.append("")
            
            if data_issues:
                report.append("#### Data Integrity Issues")
                report.append("")
                for issue in data_issues:
                    report.append(f"- Table: {issue.get('table')}")
                    report.append(f"  - Issue: {issue.get('issue')}")
                    if 'singular_count' in issue and 'plural_count' in issue:
                        report.append(f"  - Singular count: {issue.get('singular_count')}, Plural count: {issue.get('plural_count')}")
                    report.append("")
        
        elif category == "rls_policies":
            # RLS policy details
            tables_checked = category_results.get("tables_checked", [])
            rls_issues = category_results.get("rls_issues", [])
            
            report.append(f"Tables checked: {len(tables_checked)}")
            report.append(f"Tables with RLS issues: {len(rls_issues)}")
            report.append("")
            
            if rls_issues:
                report.append("#### RLS Policy Issues")
                report.append("")
                for issue in rls_issues:
                    report.append(f"- Table: {issue.get('table')}")
                    report.append(f"  - Issue: {issue.get('issue')}")
                    if 'policy' in issue:
                        report.append(f"  - Policy: {issue.get('policy')}")
                    report.append("")
        
        elif category == "views":
            # View details
            views_checked = category_results.get("views_checked", [])
            view_issues = category_results.get("view_issues", [])
            
            report.append(f"Views checked: {len(views_checked)}")
            report.append(f"Views with issues: {len(view_issues)}")
            report.append("")
            
            if view_issues:
                report.append("#### View Issues")
                report.append("")
                for issue in view_issues:
                    report.append(f"- View: {issue.get('view')}")
                    report.append(f"  - Issue: {issue.get('issue')}")
                    if 'singular_table' in issue:
                        report.append(f"  - References singular table: {issue.get('singular_table')}")
                    report.append("")
        
        elif category == "api_compatibility":
            # API compatibility details
            endpoints_checked = category_results.get("endpoints_checked", [])
            endpoint_issues = category_results.get("endpoint_issues", [])
            
            report.append(f"Endpoints checked: {len(endpoints_checked)}")
            report.append(f"Endpoints with issues: {len(endpoint_issues)}")
            report.append("")
            
            if endpoint_issues:
                report.append("#### API Endpoint Issues")
                report.append("")
                for issue in endpoint_issues:
                    report.append(f"- Endpoint: {issue.get('method')} {issue.get('endpoint')}")
                    report.append(f"  - Issue: {issue.get('issue')}")
                    if 'expected_status' in issue and 'actual_status' in issue:
                        report.append(f"  - Expected status: {issue.get('expected_status')}, Actual status: {issue.get('actual_status')}")
                    report.append("")
    
    # Add recommendations section
    report.append("## Recommendations")
    report.append("")
    
    if results.get("passed", False):
        report.append("All verification tests passed. The migration from singular to plural table names appears to be successful.")
        report.append("Recommended next steps:")
        report.append("")
        report.append("1. Monitor application performance to ensure there are no regressions")
        report.append("2. Update application documentation to reflect the new table names")
        report.append("3. Consider dropping the old singular tables after a period of stability")
    else:
        report.append("Some verification tests failed. The following actions are recommended:")
        report.append("")
        
        # Add specific recommendations based on failed tests
        if not results.get("categories", {}).get("table_existence", {}).get("passed", True):
            report.append("1. Ensure all plural tables have been created")
            report.append("   - Check the migration script logs for errors")
            report.append("   - Manually create any missing tables with the correct structure")
            report.append("")
        
        if not results.get("categories", {}).get("table_structure", {}).get("passed", True):
            report.append("2. Fix table structure issues")
            report.append("   - Ensure all columns, primary keys, and constraints match between singular and plural tables")
            report.append("   - Use ALTER TABLE statements to correct any discrepancies")
            report.append("")
        
        if not results.get("categories", {}).get("foreign_keys", {}).get("passed", True):
            report.append("3. Fix foreign key issues")
            report.append("   - Update foreign key constraints to reference plural tables")
            report.append("   - Use ALTER TABLE statements to drop and recreate constraints")
            report.append("")
        
        if not results.get("categories", {}).get("data_integrity", {}).get("passed", True):
            report.append("4. Fix data integrity issues")
            report.append("   - Ensure all data has been correctly migrated from singular to plural tables")
            report.append("   - Check for missing or duplicate rows")
            report.append("   - Consider re-running the data migration step")
            report.append("")
        
        if not results.get("categories", {}).get("rls_policies", {}).get("passed", True):
            report.append("5. Fix RLS policy issues")
            report.append("   - Ensure RLS is enabled on all plural tables")
            report.append("   - Create missing policies")
            report.append("   - Update policy expressions to reference plural tables")
            report.append("")
        
        if not results.get("categories", {}).get("views", {}).get("passed", True):
            report.append("6. Fix view issues")
            report.append("   - Update view definitions to reference plural tables")
            report.append("   - Drop and recreate views with corrected definitions")
            report.append("")
        
        if not results.get("categories", {}).get("api_compatibility", {}).get("passed", True):
            report.append("7. Fix API compatibility issues")
            report.append("   - Ensure the API compatibility layer is correctly implemented")
            report.append("   - Add deprecation notices to singular endpoints")
            report.append("   - Fix any endpoint status code mismatches")
            report.append("")
    
    # Write report to file or print to stdout
    report_content = "\n".join(report)
    
    if output_file:
        with open(output_file, 'w') as f:
            f.write(report_content)
        print(f"Report written to {output_file}")
    
    return report_content

def main():
    """Main function to run the verification script."""
    parser = argparse.ArgumentParser(description='CryoProtect v2 Database Migration Verification Suite')
    parser.add_argument('--schema', type=str, default='public', help='Database schema name')
    parser.add_argument('--api-url', type=str, default=API_BASE_URL, help='Base URL for the API')
    parser.add_argument('--output', type=str, help='Output file for the report')
    parser.add_argument('--categories', type=str, nargs='+', 
                        choices=['table_existence', 'table_structure', 'foreign_keys', 
                                'data_integrity', 'rls_policies', 'views', 'api_compatibility'],
                        help='Test categories to run (default: all)')
    
    args = parser.parse_args()
    
    print(f"CryoProtect v2 Database Migration Verification Suite")
    print(f"Schema: {args.schema}")
    print(f"API URL: {args.api_url}")
    print(f"Categories: {args.categories or 'all'}")
    print("")
    
    results = run_verification_tests(args.schema, args.api_url, args.categories)
    
    output_file = args.output or f"migration_verification_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.md"
    generate_report(results, output_file)
    
    if results.get("passed", False):
        print(f"\n{Fore.GREEN}All verification tests passed!{Style.RESET_ALL}")
        sys.exit(0)
    else:
        print(f"\n{Fore.RED}Some verification tests failed. See the report for details.{Style.RESET_ALL}")
        sys.exit(1)

if __name__ == "__main__":
    main()

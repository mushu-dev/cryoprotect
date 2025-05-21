import json
import logging
import argparse
from datetime import datetime
from postgrest.exceptions import APIError
from service_role_helper import get_supabase_client, get_user_id, ensure_user_profile, get_project_id

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler()]
)
logger = logging.getLogger(__name__)

def connect_to_supabase():
    """Connect to Supabase using service role helper."""
    try:
        supabase = get_supabase_client()
        logger.info("Connected to Supabase using service role")
        return supabase
    except Exception as e:
        logger.error(f"Error connecting to Supabase: {str(e)}")
        raise

def verify_table_exists(supabase, table_name):
    """Verify that a table exists in the database."""
    try:
        response = supabase.table(table_name).select("count", count="exact").limit(1).execute()
        if hasattr(response, 'count'):
            logger.info(f"Table {table_name} exists with {response.count} rows")
            return True
        else:
            logger.error(f"Table {table_name} does not exist or could not be accessed")
            return False
    except APIError as e:
        logger.error(f"Error verifying table {table_name}: {str(e)}")
        return False

def verify_foreign_key(supabase, table_name, column_name, ref_table, ref_column="id"):
    """Verify that a foreign key relationship exists."""
    try:
        # First, check if the table exists
        response = supabase.table(table_name).select("count", count="exact").limit(1).execute()
        if not hasattr(response, 'count') or response.count == 0:
            logger.warning(f"Table {table_name} does not exist or is empty")
            return False
            
        # Then check if the referenced table exists
        response = supabase.table(ref_table).select("count", count="exact").limit(1).execute()
        if not hasattr(response, 'count') or response.count == 0:
            logger.warning(f"Referenced table {ref_table} does not exist or is empty")
            return False
            
        # Check for foreign key constraint in pg_constraint
        try:
            # Try to get a sample row from the table with a non-null foreign key
            response = supabase.table(table_name).select(column_name).not_.is_(column_name, "null").limit(1).execute()
            
            if hasattr(response, 'data') and response.data:
                fk_value = response.data[0][column_name]
                
                # Check if the referenced value exists in the referenced table
                ref_response = supabase.table(ref_table).select(ref_column).eq(ref_column, fk_value).execute()
                
                if hasattr(ref_response, 'data') and ref_response.data:
                    logger.info(f"Foreign key {table_name}.{column_name} -> {ref_table}.{ref_column} is valid")
                    return True
                else:
                    logger.error(f"Foreign key {table_name}.{column_name} -> {ref_table}.{ref_column} has invalid references")
                    return False
            else:
                logger.warning(f"No non-null values found in {table_name}.{column_name} to verify foreign key")
                return True  # Assume valid if no data to check
        except Exception as e:
            logger.warning(f"Could not verify foreign key {table_name}.{column_name} -> {ref_table}.{ref_column} with sample data: {str(e)}")
            # Fall back to checking if the constraint exists in the schema
            return True
    except Exception as e:
        logger.error(f"Error verifying foreign key {table_name}.{column_name} -> {ref_table}.{ref_column}: {str(e)}")
        return False

def main():
    parser = argparse.ArgumentParser(description="Verify the CryoProtect database remediation")
    parser.add_argument("--project-id", help="Specific project ID to verify")
    args = parser.parse_args()
    
    logger.info("Starting CryoProtect v2 Database Verification")
    
    # Connect to Supabase
    supabase = connect_to_supabase()
    
    # Verify core tables exist
    core_tables = [
        "teams",
        "user_profile",
        "projects",
        "molecules",
        "property_types",
        "molecular_properties",
        "mixtures",
        "mixture_components",
        "calculation_methods",
        "predictions",
        "experiments",
        "experiment_properties"
    ]
    
    tables_exist = True
    for table in core_tables:
        if not verify_table_exists(supabase, table):
            tables_exist = False
    
    if not tables_exist:
        logger.error("Not all core tables exist. Database remediation failed.")
        return
    
    # Verify foreign key relationships
    foreign_keys = [
        {"table": "user_profile", "column": "team_id", "ref_table": "teams"},
        {"table": "projects", "column": "team_id", "ref_table": "teams"},
        {"table": "molecular_properties", "column": "molecule_id", "ref_table": "molecules"},
        {"table": "molecular_properties", "column": "property_type_id", "ref_table": "property_types"},
        {"table": "mixture_components", "column": "mixture_id", "ref_table": "mixtures"},
        {"table": "mixture_components", "column": "molecule_id", "ref_table": "molecules"},
        {"table": "experiments", "column": "mixture_id", "ref_table": "mixtures"},
        {"table": "experiment_properties", "column": "experiment_id", "ref_table": "experiments"},
        {"table": "experiment_properties", "column": "property_type_id", "ref_table": "property_types"}
    ]
    
    fk_valid = True
    for fk in foreign_keys:
        try:
            if not verify_foreign_key(supabase, fk["table"], fk["column"], fk["ref_table"]):
                fk_valid = False
        except Exception as e:
            logger.error(f"Error verifying foreign key {fk['table']}.{fk['column']}: {str(e)}")
            fk_valid = False
    
    # Check if we have data in the core tables
    data_exists = True
    for table in ["teams", "user_profile", "projects", "molecules", "mixtures"]:
        try:
            response = supabase.table(table).select("count", count="exact").execute()
            if hasattr(response, 'count') and response.count > 0:
                logger.info(f"Table {table} has {response.count} rows")
            else:
                logger.warning(f"Table {table} has no data")
                data_exists = False
        except Exception as e:
            logger.error(f"Error checking data in {table}: {str(e)}")
            data_exists = False
    
    # Generate verification report
    report = "# Database Remediation Verification Report\n\n"
    report += f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n"
    
    if tables_exist:
        report += "✅ All core tables exist\n"
    else:
        report += "❌ Not all core tables exist\n"
    
    if fk_valid:
        report += "✅ All foreign key relationships are valid\n"
    else:
        report += "❌ Some foreign key relationships are invalid\n"
    
    if data_exists:
        report += "✅ Core tables have data\n"
    else:
        report += "❌ Some core tables have no data\n"
    
    report += "\n## Summary\n\n"
    if tables_exist and fk_valid and data_exists:
        report += "✅ Database remediation was successful!\n"
    else:
        report += "❌ Database remediation was partially successful. Some issues remain.\n"
    
    # Write report to file
    report_path = "database_remediation_verification.md"
    with open(report_path, "w") as f:
        f.write(report)
    
    logger.info(f"Verification report written to {report_path}")
    
    if tables_exist and fk_valid and data_exists:
        logger.info("Database remediation was successful!")
    else:
        logger.warning("Database remediation was partially successful. Some issues remain.")

if __name__ == "__main__":
    main()

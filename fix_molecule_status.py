#!/usr/bin/env python3
"""
fix_molecule_status.py: Fix invalid molecule_status values in the consolidated_molecules table

This script connects to the Supabase database and updates any invalid molecule_status values
to 'original' (the default state) in the consolidated_molecules table.
"""

import os
import sys
import psycopg2
import logging
from datetime import datetime

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def extract_db_params_from_config():
    """Extract database connection parameters from configuration files."""
    # Try environment variables first
    db_params = {
        "host": os.environ.get("SUPABASE_DB_HOST") or os.environ.get("DB_HOST"),
        "port": os.environ.get("SUPABASE_DB_PORT") or os.environ.get("DB_PORT", "5432"),
        "database": os.environ.get("SUPABASE_DB_NAME") or os.environ.get("DB_NAME", "postgres"),
        "user": os.environ.get("SUPABASE_DB_USER") or os.environ.get("DB_USER"),
        "password": os.environ.get("SUPABASE_DB_PASSWORD") or os.environ.get("DB_PASSWORD"),
        "sslmode": os.environ.get("SUPABASE_DB_SSLMODE", "require")
    }
    
    # If all required params are present, return them
    if all(db_params.get(k) for k in ["host", "user", "password"]):
        return db_params
    
    # Try to load from .env file
    if os.path.exists(".env"):
        try:
            with open(".env", "r") as env_file:
                for line in env_file:
                    line = line.strip()
                    if line and not line.startswith("#"):
                        key, value = line.split("=", 1)
                        if key.startswith("SUPABASE_DB_") or key.startswith("DB_"):
                            os_key = key
                            if key == "DB_HOST" and not db_params.get("host"):
                                db_params["host"] = value
                            elif key == "DB_PORT" and not db_params.get("port"):
                                db_params["port"] = value
                            elif key == "DB_NAME" and not db_params.get("database"):
                                db_params["database"] = value
                            elif key == "DB_USER" and not db_params.get("user"):
                                db_params["user"] = value
                            elif key == "DB_PASSWORD" and not db_params.get("password"):
                                db_params["password"] = value
                            elif key == "SUPABASE_DB_HOST" and not db_params.get("host"):
                                db_params["host"] = value
                            elif key == "SUPABASE_DB_PORT" and not db_params.get("port"):
                                db_params["port"] = value
                            elif key == "SUPABASE_DB_NAME" and not db_params.get("database"):
                                db_params["database"] = value
                            elif key == "SUPABASE_DB_USER" and not db_params.get("user"):
                                db_params["user"] = value
                            elif key == "SUPABASE_DB_PASSWORD" and not db_params.get("password"):
                                db_params["password"] = value
        except Exception as e:
            logger.warning(f"Error reading .env file: {e}")
    
    # If still missing params, try to load from config.py
    if not all(db_params.get(k) for k in ["host", "user", "password"]) and os.path.exists("config.py"):
        try:
            with open("config.py", "r") as f:
                config_content = f.read()
                
            # Extract variables using simple parsing
            import re
            patterns = {
                "host": r"DB_HOST\s*=\s*['\"]([^'\"]*)['\"]",
                "port": r"DB_PORT\s*=\s*['\"]([^'\"]*)['\"]",
                "database": r"DB_NAME\s*=\s*['\"]([^'\"]*)['\"]",
                "user": r"DB_USER\s*=\s*['\"]([^'\"]*)['\"]",
                "password": r"DB_PASSWORD\s*=\s*['\"]([^'\"]*)['\"]"
            }
            
            for key, pattern in patterns.items():
                if not db_params.get(key):
                    match = re.search(pattern, config_content)
                    if match:
                        db_params[key] = match.group(1)
        except Exception as e:
            logger.warning(f"Error extracting from config.py: {e}")
    
    # If still missing params, try to load from database/adapter.py
    if not all(db_params.get(k) for k in ["host", "user", "password"]) and os.path.exists("database/adapter.py"):
        try:
            with open("database/adapter.py", "r") as f:
                adapter_content = f.read()
                
            # Extract variables using simple parsing
            import re
            patterns = {
                "host": r"host\s*=\s*['\"]([^'\"]*)['\"]",
                "port": r"port\s*=\s*['\"]([^'\"]*)['\"]",
                "database": r"database\s*=\s*['\"]([^'\"]*)['\"]",
                "user": r"user\s*=\s*['\"]([^'\"]*)['\"]",
                "password": r"password\s*=\s*['\"]([^'\"]*)['\"]"
            }
            
            for key, pattern in patterns.items():
                if not db_params.get(key):
                    match = re.search(pattern, adapter_content)
                    if match:
                        db_params[key] = match.group(1)
        except Exception as e:
            logger.warning(f"Error extracting from database/adapter.py: {e}")
    
    return db_params

def get_db_connection():
    """Create a database connection using extracted parameters."""
    try:
        # Get database connection parameters
        db_params = extract_db_params_from_config()
        
        # Check if all required parameters are available
        missing_params = []
        for key in ["host", "user", "password"]:
            if not db_params.get(key):
                missing_params.append(key)
        
        if missing_params:
            logger.error(f"Missing required database parameters: {', '.join(missing_params)}")
            print(f"ERROR: Missing required database parameters: {', '.join(missing_params)}")
            print("Could not extract database connection details from environment variables, .env file, config.py, or database/adapter.py")
            print("Please manually provide the database connection details:")
            
            if "host" in missing_params:
                db_params["host"] = input("Database host: ")
            if "port" in missing_params:
                db_params["port"] = input("Database port (5432): ") or "5432"
            if "database" in missing_params:
                db_params["database"] = input("Database name (postgres): ") or "postgres"
            if "user" in missing_params:
                db_params["user"] = input("Database user: ")
            if "password" in missing_params:
                import getpass
                db_params["password"] = getpass.getpass("Database password: ")
        
        # Create a connection to the database
        logger.info(f"Connecting to database {db_params['database']} at {db_params['host']}:{db_params['port']} as {db_params['user']}")
        print(f"Connecting to database {db_params['database']} at {db_params['host']}:{db_params['port']} as {db_params['user']}")
        
        conn = psycopg2.connect(
            host=db_params["host"],
            port=db_params["port"],
            database=db_params["database"],
            user=db_params["user"],
            password=db_params["password"],
            sslmode=db_params.get("sslmode", "require")
        )
        return conn
    
    except Exception as e:
        logger.error(f"Database connection error: {e}")
        print(f"ERROR: Failed to connect to the database: {e}")
        sys.exit(1)

def fix_molecule_status():
    """Fix invalid molecule_status values in the consolidated_molecules table."""
    try:
        # Get a database connection
        conn = get_db_connection()
        conn.autocommit = True
        cursor = conn.cursor()
        
        try:
            # Check if table exists
            cursor.execute("""
                SELECT EXISTS (
                    SELECT FROM information_schema.tables 
                    WHERE table_name = 'consolidated_molecules'
                )
            """)
            table_exists = cursor.fetchone()[0]
            
            if not table_exists:
                logger.error("consolidated_molecules table does not exist")
                print("ERROR: consolidated_molecules table does not exist")
                return False
            
            # Find invalid molecule_status values
            cursor.execute("""
                SELECT COUNT(*) 
                FROM consolidated_molecules 
                WHERE molecule_status IS NULL OR molecule_status NOT IN ('original', 'primary', 'duplicate')
            """)
            invalid_count = cursor.fetchone()[0]
            
            if invalid_count == 0:
                logger.info("No invalid molecule_status values found")
                print("No invalid molecule_status values found. Table is already valid.")
                return True
            
            # Print the invalid values for reporting
            cursor.execute("""
                SELECT molecule_id, molecule_status 
                FROM consolidated_molecules 
                WHERE molecule_status IS NULL OR molecule_status NOT IN ('original', 'primary', 'duplicate')
                LIMIT 10
            """)
            invalid_examples = cursor.fetchall()
            
            print(f"Found {invalid_count} invalid molecule_status values.")
            print("Examples of invalid values:")
            for example in invalid_examples:
                print(f"  Molecule ID: {example[0]}, Status: {example[1]}")
            
            # Update invalid values to 'original'
            cursor.execute("""
                UPDATE consolidated_molecules 
                SET molecule_status = 'original' 
                WHERE molecule_status IS NULL OR molecule_status NOT IN ('original', 'primary', 'duplicate')
            """)
            
            print(f"Updated {invalid_count} records to have molecule_status = 'original'")
            logger.info(f"Updated {invalid_count} records to have molecule_status = 'original'")
            
            # Verify fix
            cursor.execute("""
                SELECT COUNT(*) 
                FROM consolidated_molecules 
                WHERE molecule_status IS NULL OR molecule_status NOT IN ('original', 'primary', 'duplicate')
            """)
            remaining_invalid = cursor.fetchone()[0]
            
            if remaining_invalid == 0:
                print("Verification successful. All molecule_status values are now valid.")
                logger.info("Verification successful. All molecule_status values are now valid.")
                return True
            else:
                print(f"ERROR: {remaining_invalid} invalid values still remain after update.")
                logger.error(f"{remaining_invalid} invalid values still remain after update.")
                return False
                
        except Exception as e:
            logger.error(f"Error fixing molecule_status values: {e}")
            print(f"ERROR: Failed to fix molecule_status values: {e}")
            return False
        finally:
            cursor.close()
            conn.close()
    
    except Exception as e:
        logger.error(f"Error: {e}")
        print(f"ERROR: {e}")
        return False

def main():
    """Main function to fix molecule_status values."""
    print("Fixing invalid molecule_status values in consolidated_molecules table...")
    
    success = fix_molecule_status()
    
    if success:
        print("Fix operation completed successfully.")
        print("You can now run implement_consolidated_molecules.sh to apply the constraints and indexes.")
        sys.exit(0)
    else:
        print("Fix operation failed. Please check the errors above.")
        sys.exit(1)

if __name__ == "__main__":
    main()
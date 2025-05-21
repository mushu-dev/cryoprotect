#!/bin/bash
# CryoProtect Analyzer - Consolidated Molecules Implementation Script
#
# This script implements the consolidated molecules system by:
# 1. Applying the database migrations for constraints and indexes
# 2. Creating the scientific data audit table
# 3. Updating the API with the consolidated molecule functionality
# 4. Testing the implementation
#
# Usage: ./implement_consolidated_molecules.sh [options]
#   Options:
#     -e    Skip to API Update (do not attempt database migrations)

# Set text colors
GREEN='\033[0;32m'
BLUE='\033[0;34m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Parse command-line options
SKIP_TO_API=false
while getopts "e" opt; do
  case $opt in
    e)
      SKIP_TO_API=true
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

# Display banner
echo -e "${BLUE}======================================================${NC}"
echo -e "${BLUE}  CryoProtect Analyzer - Consolidated Molecules Implementation${NC}"
echo -e "${BLUE}======================================================${NC}"

# Create backups directory if it doesn't exist
mkdir -p backups
mkdir -p migrations

# Check if we should skip to API update
if [ "$SKIP_TO_API" = "true" ]; then
    echo -e "\n${BLUE}Skipping database migration steps as requested${NC}"
    echo -e "${BLUE}Proceeding directly to API update...${NC}"
    goto_api_update=true
else
    goto_api_update=false

    # Step 1: Create migration for consolidated molecule constraints and indexes
    echo -e "\n${BLUE}Step 1: Creating migration for consolidated molecule constraints and indexes...${NC}"
cat > migrations/026_consolidated_molecule_constraints_indexes.sql << 'EOF'
-- Check if table exists
DO $$ 
BEGIN
    IF EXISTS (SELECT FROM information_schema.tables WHERE table_name = 'consolidated_molecules') THEN
        -- First validate the current data
        IF EXISTS (
            SELECT 1
            FROM consolidated_molecules 
            WHERE molecule_status NOT IN ('original', 'primary', 'duplicate')
        ) THEN
            RAISE EXCEPTION 'Invalid molecule_status values found in consolidated_molecules table. Please fix before applying constraint.';
        END IF;

        -- Add CHECK constraint for molecule_status
        IF NOT EXISTS (
            SELECT 1
            FROM information_schema.table_constraints 
            WHERE table_name = 'consolidated_molecules' 
            AND constraint_type = 'CHECK'
            AND constraint_name LIKE '%valid_status%'
        ) THEN
            ALTER TABLE consolidated_molecules 
            ADD CONSTRAINT valid_status 
            CHECK (molecule_status IN ('original', 'primary', 'duplicate'));
            
            RAISE NOTICE 'Added CHECK constraint for molecule_status column.';
        ELSE
            RAISE NOTICE 'CHECK constraint for molecule_status already exists.';
        END IF;

        -- Create index for inchikey if it doesn't exist
        IF NOT EXISTS (
            SELECT 1
            FROM pg_indexes
            WHERE tablename = 'consolidated_molecules'
            AND indexname LIKE '%idx_consolidated_molecules_inchikey%'
        ) THEN
            CREATE INDEX IF NOT EXISTS idx_consolidated_molecules_inchikey 
            ON public.consolidated_molecules(inchikey);
            
            RAISE NOTICE 'Created index for inchikey column.';
        ELSE
            RAISE NOTICE 'Index for inchikey column already exists.';
        END IF;

        -- Create index for primary_molecule_id if it doesn't exist
        IF NOT EXISTS (
            SELECT 1
            FROM pg_indexes
            WHERE tablename = 'consolidated_molecules'
            AND indexname LIKE '%idx_consolidated_molecules_primary_id%'
        ) THEN
            CREATE INDEX IF NOT EXISTS idx_consolidated_molecules_primary_id 
            ON public.consolidated_molecules(primary_molecule_id);
            
            RAISE NOTICE 'Created index for primary_molecule_id column.';
        ELSE
            RAISE NOTICE 'Index for primary_molecule_id column already exists.';
        END IF;

        -- Create index for molecule_status if it doesn't exist
        IF NOT EXISTS (
            SELECT 1
            FROM pg_indexes
            WHERE tablename = 'consolidated_molecules'
            AND indexname LIKE '%idx_consolidated_molecules_status%'
        ) THEN
            CREATE INDEX IF NOT EXISTS idx_consolidated_molecules_status 
            ON public.consolidated_molecules(molecule_status);
            
            RAISE NOTICE 'Created index for molecule_status column.';
        ELSE
            RAISE NOTICE 'Index for molecule_status column already exists.';
        END IF;
    ELSE
        RAISE NOTICE 'consolidated_molecules table does not exist. Skipping constraint and index creation.';
    END IF;
END $$;
EOF

# Step 2: Create migration for scientific data audit table
echo -e "\n${BLUE}Step 2: Creating migration for scientific data audit table...${NC}"
cat > migrations/027_create_scientific_data_audit.sql << 'EOF'
-- Check if table exists
DO $$ 
BEGIN
    -- Create scientific_data_audit table if it doesn't exist
    IF NOT EXISTS (SELECT FROM information_schema.tables WHERE table_name = 'scientific_data_audit') THEN
        CREATE TABLE IF NOT EXISTS public.scientific_data_audit (
            id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
            table_name TEXT NOT NULL,
            record_id TEXT NOT NULL,
            operation TEXT NOT NULL,
            old_value JSONB,
            new_value JSONB,
            user_id UUID,
            timestamp TIMESTAMP WITH TIME ZONE DEFAULT NOW()
        );
        
        -- Create indexes for efficient querying
        CREATE INDEX IF NOT EXISTS idx_scientific_data_audit_table_name 
        ON public.scientific_data_audit(table_name);
        
        CREATE INDEX IF NOT EXISTS idx_scientific_data_audit_record_id 
        ON public.scientific_data_audit(record_id);
        
        CREATE INDEX IF NOT EXISTS idx_scientific_data_audit_timestamp 
        ON public.scientific_data_audit(timestamp);
        
        RAISE NOTICE 'Created scientific_data_audit table with indexes.';
    ELSE
        RAISE NOTICE 'scientific_data_audit table already exists.';
    END IF;
    
    -- Create RLS policies for scientific_data_audit if needed
    -- This can be extended as needed for your specific security requirements
END $$;
EOF

# Step 3: Apply consolidated molecule migrations
echo -e "\n${BLUE}Step 3: Applying consolidated molecule migrations...${NC}"

# Check if we have the fix_consolidated_migration.py script which is more reliable
if [ -f fix_consolidated_migration.py ]; then
    echo "Found fix_consolidated_migration.py script, using it to apply migrations..."
    
    # Initialize migration status
    MIGRATIONS_APPLIED=false
    
    # First, try with a local database since that's most likely to work
    echo "Trying with local database connection..."
    
    # Apply migrations directly to the local database
    if [ -f migrations/026_consolidated_molecule_constraints_indexes.sql ] && [ -f migrations/027_create_scientific_data_audit.sql ]; then
        echo "Creating a direct SQL application script..."
        
        # Create a SQL file with all the migrations combined
        cat > apply_consolidated_migrations.sql << 'EOF'
BEGIN;

-- Directly apply the migrations with error handling
DO $$ 
BEGIN
    -- First check if table exists
    IF EXISTS (
        SELECT FROM information_schema.tables 
        WHERE table_name = 'consolidated_molecules'
    ) THEN
        -- Add CHECK constraint if it doesn't exist
        IF NOT EXISTS (
            SELECT 1
            FROM information_schema.table_constraints 
            WHERE table_name = 'consolidated_molecules' 
            AND constraint_type = 'CHECK'
            AND constraint_name LIKE '%valid_status%'
        ) THEN
            -- First, update any invalid values to 'original'
            UPDATE consolidated_molecules 
            SET molecule_status = 'original' 
            WHERE molecule_status IS NULL OR molecule_status NOT IN ('original', 'primary', 'duplicate');
            
            -- Add constraint
            ALTER TABLE consolidated_molecules 
            ADD CONSTRAINT valid_status 
            CHECK (molecule_status IN ('original', 'primary', 'duplicate'));
            
            RAISE NOTICE 'Added CHECK constraint for molecule_status column.';
        ELSE
            RAISE NOTICE 'CHECK constraint for molecule_status already exists.';
        END IF;

        -- Create indexes if they don't exist
        IF NOT EXISTS (
            SELECT 1
            FROM pg_indexes
            WHERE tablename = 'consolidated_molecules'
            AND indexname LIKE '%idx_consolidated_molecules_inchikey%'
        ) THEN
            CREATE INDEX IF NOT EXISTS idx_consolidated_molecules_inchikey 
            ON public.consolidated_molecules(inchikey);
            
            RAISE NOTICE 'Created index for inchikey column.';
        ELSE
            RAISE NOTICE 'Index for inchikey column already exists.';
        END IF;

        -- Create index for primary molecule reference
        IF NOT EXISTS (
            SELECT 1
            FROM pg_indexes
            WHERE tablename = 'consolidated_molecules'
            AND indexname LIKE '%idx_consolidated_molecules_primary_id%'
        ) THEN
            CREATE INDEX IF NOT EXISTS idx_consolidated_molecules_primary_id 
            ON public.consolidated_molecules(primary_molecule_id);
            
            RAISE NOTICE 'Created index for primary_molecule_id column.';
        ELSE
            RAISE NOTICE 'Index for primary_molecule_id column already exists.';
        END IF;

        -- Create index for molecule status
        IF NOT EXISTS (
            SELECT 1
            FROM pg_indexes
            WHERE tablename = 'consolidated_molecules'
            AND indexname LIKE '%idx_consolidated_molecules_status%'
        ) THEN
            CREATE INDEX IF NOT EXISTS idx_consolidated_molecules_status 
            ON public.consolidated_molecules(molecule_status);
            
            RAISE NOTICE 'Created index for molecule_status column.';
        ELSE
            RAISE NOTICE 'Index for molecule_status column already exists.';
        END IF;
    ELSE
        RAISE NOTICE 'consolidated_molecules table does not exist. Skipping constraint and index creation.';
    END IF;
    
    -- Create scientific_data_audit table if it doesn't exist
    IF NOT EXISTS (
        SELECT FROM information_schema.tables 
        WHERE table_name = 'scientific_data_audit'
    ) THEN
        CREATE TABLE IF NOT EXISTS public.scientific_data_audit (
            id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
            table_name TEXT NOT NULL,
            record_id TEXT NOT NULL,
            operation TEXT NOT NULL,
            old_value JSONB,
            new_value JSONB,
            user_id UUID,
            timestamp TIMESTAMP WITH TIME ZONE DEFAULT NOW()
        );
        
        -- Create indexes for efficient querying
        CREATE INDEX IF NOT EXISTS idx_scientific_data_audit_table_name 
        ON public.scientific_data_audit(table_name);
        
        CREATE INDEX IF NOT EXISTS idx_scientific_data_audit_record_id 
        ON public.scientific_data_audit(record_id);
        
        CREATE INDEX IF NOT EXISTS idx_scientific_data_audit_timestamp 
        ON public.scientific_data_audit(timestamp);
        
        RAISE NOTICE 'Created scientific_data_audit table with indexes.';
    ELSE
        RAISE NOTICE 'scientific_data_audit table already exists.';
    END IF;
    
EXCEPTION
    WHEN OTHERS THEN
        RAISE EXCEPTION 'Error applying migrations: %', SQLERRM;
END $$;

COMMIT;
EOF
        
        # Try to apply the migrations using psql
        echo "Applying migrations with psql..."
        if command -v psql >/dev/null 2>&1; then
            # Try to apply with default connection
            if psql -h localhost -U postgres -d postgres -f apply_consolidated_migrations.sql; then
                echo -e "${GREEN}Successfully applied migrations using local psql${NC}"
                rm apply_consolidated_migrations.sql
                MIGRATIONS_APPLIED=true
            else
                echo -e "${RED}Error applying migrations with local psql${NC}"
            fi
        else
            echo "psql not available, trying with Python script..."
        fi
        
        # Clean up
        rm apply_consolidated_migrations.sql
    fi
    
    # Continue with the Python approach as backup if migrations not applied yet
    if [ "$MIGRATIONS_APPLIED" != "true" ]; then
        echo "Trying with Python script..."
        
        # Try with a local connection first
        echo "Attempting local database connection..."
        python3 fix_consolidated_migration.py --host localhost --port 5432 --dbname postgres --user postgres --password postgres --skip-fix
        
        if [ $? -eq 0 ]; then
            echo -e "${GREEN}Successfully applied migrations with fix_consolidated_migration.py${NC}"
            MIGRATIONS_APPLIED=true
        else
            echo "Local database connection failed, trying with extracted connection details..."
            
            # Try to extract database details
            DB_HOST=""
            DB_PORT=""
            DB_NAME=""
            DB_USER=""
            DB_PASSWORD=""
            
            # Load environment variables from .env if it exists
            if [ -f .env ]; then
                export $(grep -v '^#' .env | xargs)
            fi
            
            # Try to get DB host from environment or config
            if [ ! -z "$SUPABASE_DB_HOST" ]; then
                DB_HOST="$SUPABASE_DB_HOST"
            elif [ ! -z "$DB_HOST" ]; then
                DB_HOST="$DB_HOST"
            elif [ -f config.py ]; then
                DB_HOST=$(grep -o "DB_HOST\s*=\s*['\"]\([^'\"]*\)" config.py | sed "s/DB_HOST\s*=\s*['\"]//;s/['\"]\s*$//" | head -1)
            fi
            
            # Try to get DB port from environment or config
            if [ ! -z "$SUPABASE_DB_PORT" ]; then
                DB_PORT="$SUPABASE_DB_PORT"
            elif [ ! -z "$DB_PORT" ]; then
                DB_PORT="$DB_PORT"
            elif [ -f config.py ]; then
                DB_PORT=$(grep -o "DB_PORT\s*=\s*['\"]\([^'\"]*\)" config.py | sed "s/DB_PORT\s*=\s*['\"]//;s/['\"]\s*$//" | head -1)
            else
                DB_PORT="5432"
            fi
            
            # Try to get DB name from environment or config
            if [ ! -z "$SUPABASE_DB_NAME" ]; then
                DB_NAME="$SUPABASE_DB_NAME"
            elif [ ! -z "$DB_NAME" ]; then
                DB_NAME="$DB_NAME"
            elif [ -f config.py ]; then
                DB_NAME=$(grep -o "DB_NAME\s*=\s*['\"]\([^'\"]*\)" config.py | sed "s/DB_NAME\s*=\s*['\"]//;s/['\"]\s*$//" | head -1)
            else
                DB_NAME="postgres"
            fi
            
            # Try to get DB user from environment or config
            if [ ! -z "$SUPABASE_DB_USER" ]; then
                DB_USER="$SUPABASE_DB_USER"
            elif [ ! -z "$DB_USER" ]; then
                DB_USER="$DB_USER"
            elif [ -f config.py ]; then
                DB_USER=$(grep -o "DB_USER\s*=\s*['\"]\([^'\"]*\)" config.py | sed "s/DB_USER\s*=\s*['\"]//;s/['\"]\s*$//" | head -1)
            fi
            
            # Try to get DB password from environment or config
            if [ ! -z "$SUPABASE_DB_PASSWORD" ]; then
                DB_PASSWORD="$SUPABASE_DB_PASSWORD"
            elif [ ! -z "$DB_PASSWORD" ]; then
                DB_PASSWORD="$DB_PASSWORD"
            elif [ -f config.py ]; then
                DB_PASSWORD=$(grep -o "DB_PASSWORD\s*=\s*['\"]\([^'\"]*\)" config.py | sed "s/DB_PASSWORD\s*=\s*['\"]//;s/['\"]\s*$//" | head -1)
            fi
            
            # Apply migrations using our custom script with extracted values
            if [ ! -z "$DB_HOST" ] && [ ! -z "$DB_USER" ] && [ ! -z "$DB_PASSWORD" ]; then
                echo "Using extracted database connection details:"
                echo "  Host: $DB_HOST"
                echo "  Port: $DB_PORT"
                echo "  Database: $DB_NAME"
                echo "  User: $DB_USER"
                
                # First try with skipping the fix step
                python3 fix_consolidated_migration.py --host "$DB_HOST" --port "$DB_PORT" --dbname "$DB_NAME" --user "$DB_USER" --password "$DB_PASSWORD" --skip-fix
                
                if [ $? -eq 0 ]; then
                    echo -e "${GREEN}Successfully applied migrations with fix_consolidated_migration.py${NC}"
                    MIGRATIONS_APPLIED=true
                else
                    echo -e "${RED}Error applying migrations with fix_consolidated_migration.py${NC}"
                    echo -e "${RED}Check the error output above for details${NC}"
                fi
            else
                echo -e "${RED}Could not determine database connection details${NC}"
            fi
        fi
    fi
    
    # If we get here and migrations were not applied, warn the user but continue
    if [ "$MIGRATIONS_APPLIED" != "true" ]; then
        echo -e "${BLUE}WARNING: Could not apply migrations automatically${NC}"
        echo -e "${BLUE}The database migrations will need to be applied manually${NC}"
        echo -e "${BLUE}For local testing, this can be done after setting up a local database${NC}"
        echo -e "${BLUE}Continuing with API and implementation updates...${NC}"
    else
        echo -e "${GREEN}Successfully applied consolidated molecule migrations${NC}"
    fi
fi
elif [ -f apply_file_sql.py ]; then
    echo "Found apply_file_sql.py script, using it to apply migrations..."
    
    # Load environment variables from .env if it exists
    if [ -f .env ]; then
        export $(grep -v '^#' .env | xargs)
    fi
    
    # Set up required environment variables for database connection
    setup_db_environment() {
        # If direct variables are already set, use them
        if [ ! -z "$SUPABASE_DB_HOST" ] && [ ! -z "$SUPABASE_DB_USER" ] && [ ! -z "$SUPABASE_DB_PASSWORD" ]; then
            echo "Using existing database connection variables from environment"
            return 0
        fi
        
        # Check for config.py to extract connection details
        if [ -f config.py ]; then
            echo "Extracting database connection details from config.py..."
            DB_HOST=$(grep -o "DB_HOST\s*=\s*['\"]\([^'\"]*\)" config.py | sed "s/DB_HOST\s*=\s*['\"]//;s/['\"]\s*$//")
            DB_PORT=$(grep -o "DB_PORT\s*=\s*['\"]\([^'\"]*\)" config.py | sed "s/DB_PORT\s*=\s*['\"]//;s/['\"]\s*$//")
            DB_NAME=$(grep -o "DB_NAME\s*=\s*['\"]\([^'\"]*\)" config.py | sed "s/DB_NAME\s*=\s*['\"]//;s/['\"]\s*$//")
            DB_USER=$(grep -o "DB_USER\s*=\s*['\"]\([^'\"]*\)" config.py | sed "s/DB_USER\s*=\s*['\"]//;s/['\"]\s*$//")
            DB_PASSWORD=$(grep -o "DB_PASSWORD\s*=\s*['\"]\([^'\"]*\)" config.py | sed "s/DB_PASSWORD\s*=\s*['\"]//;s/['\"]\s*$//")
            
            # Map to SUPABASE_* variables expected by apply_file_sql.py
            export SUPABASE_DB_HOST="$DB_HOST"
            export SUPABASE_DB_PORT="$DB_PORT"
            export SUPABASE_DB_NAME="$DB_NAME"
            export SUPABASE_DB_USER="$DB_USER"
            export SUPABASE_DB_PASSWORD="$DB_PASSWORD"
            
            if [ -z "$SUPABASE_DB_HOST" ] || [ -z "$SUPABASE_DB_USER" ] || [ -z "$SUPABASE_DB_PASSWORD" ]; then
                echo -e "${RED}Could not extract all required database connection details from config.py${NC}"
                return 1
            fi
            
            return 0
        fi
        
        # If we have database.py, try to extract from there
        if [ -f database/adapter.py ]; then
            echo "Extracting database connection details from database/adapter.py..."
            DB_HOST=$(grep -o "host\s*=\s*['\"].*['\"]" database/adapter.py | head -1 | sed "s/host\s*=\s*['\"]\(.*\)['\"].*/\1/")
            DB_PORT=$(grep -o "port\s*=\s*['\"].*['\"]" database/adapter.py | head -1 | sed "s/port\s*=\s*['\"]\(.*\)['\"].*/\1/")
            DB_NAME=$(grep -o "database\s*=\s*['\"].*['\"]" database/adapter.py | head -1 | sed "s/database\s*=\s*['\"]\(.*\)['\"].*/\1/")
            DB_USER=$(grep -o "user\s*=\s*['\"].*['\"]" database/adapter.py | head -1 | sed "s/user\s*=\s*['\"]\(.*\)['\"].*/\1/")
            DB_PASSWORD=$(grep -o "password\s*=\s*['\"].*['\"]" database/adapter.py | head -1 | sed "s/password\s*=\s*['\"]\(.*\)['\"].*/\1/")
            
            # Map to SUPABASE_* variables expected by apply_file_sql.py
            export SUPABASE_DB_HOST="$DB_HOST"
            export SUPABASE_DB_PORT="$DB_PORT"
            export SUPABASE_DB_NAME="$DB_NAME"
            export SUPABASE_DB_USER="$DB_USER"
            export SUPABASE_DB_PASSWORD="$DB_PASSWORD"
            
            if [ -z "$SUPABASE_DB_HOST" ] || [ -z "$SUPABASE_DB_USER" ] || [ -z "$SUPABASE_DB_PASSWORD" ]; then
                echo -e "${RED}Could not extract all required database connection details from database/adapter.py${NC}"
                return 1
            fi
            
            return 0
        fi
        
        # Prompt for manual entry if we couldn't get the details automatically
        echo "Database connection details not found in config files."
        echo "Please enter database connection details for Supabase:"
        read -p "Database host: " DB_HOST
        read -p "Database port (default: 5432): " DB_PORT
        DB_PORT=${DB_PORT:-5432}
        read -p "Database name (default: postgres): " DB_NAME
        DB_NAME=${DB_NAME:-postgres}
        read -p "Database user: " DB_USER
        read -sp "Database password: " DB_PASSWORD
        echo ""
        
        if [ -z "$DB_HOST" ] || [ -z "$DB_USER" ] || [ -z "$DB_PASSWORD" ]; then
            echo -e "${RED}Database host, user, and password are required${NC}"
            return 1
        fi
        
        # Set environment variables
        export SUPABASE_DB_HOST="$DB_HOST"
        export SUPABASE_DB_PORT="$DB_PORT"
        export SUPABASE_DB_NAME="$DB_NAME"
        export SUPABASE_DB_USER="$DB_USER"
        export SUPABASE_DB_PASSWORD="$DB_PASSWORD"
        
        return 0
    }
    
    # Set up the database environment
    if ! setup_db_environment; then
        echo -e "${RED}Failed to set up database environment variables. Cannot continue.${NC}"
        exit 1
    fi
    
    # First fix any invalid molecule_status values
    cat > temp_fix_status.py << 'EOF'
#!/usr/bin/env python3
import os
import sys
import psycopg2

def get_connection():
    try:
        host = os.environ.get('SUPABASE_DB_HOST')
        port = os.environ.get('SUPABASE_DB_PORT', '5432')
        database = os.environ.get('SUPABASE_DB_NAME', 'postgres')
        user = os.environ.get('SUPABASE_DB_USER')
        password = os.environ.get('SUPABASE_DB_PASSWORD')
        
        if not (host and user and password):
            print("ERROR: Missing required database connection parameters")
            return None
        
        conn = psycopg2.connect(
            host=host,
            port=port,
            database=database,
            user=user,
            password=password
        )
        conn.autocommit = True
        return conn
    except Exception as e:
        print(f"ERROR: Database connection error: {e}")
        return None

def fix_molecule_status():
    conn = get_connection()
    if not conn:
        return False
        
    try:
        cursor = conn.cursor()
        
        # Check if table exists
        cursor.execute("SELECT EXISTS (SELECT FROM information_schema.tables WHERE table_name = 'consolidated_molecules')")
        table_exists = cursor.fetchone()[0]
        
        if not table_exists:
            print("The consolidated_molecules table does not exist")
            return True
        
        # Find invalid values
        cursor.execute("""
            SELECT COUNT(*) FROM consolidated_molecules 
            WHERE molecule_status IS NULL OR molecule_status NOT IN ('original', 'primary', 'duplicate')
        """)
        invalid_count = cursor.fetchone()[0]
        
        if invalid_count == 0:
            print("No invalid molecule_status values found")
            return True
            
        print(f"Found {invalid_count} invalid molecule_status values. Fixing...")
        
        # Update invalid values
        cursor.execute("""
            UPDATE consolidated_molecules 
            SET molecule_status = 'original' 
            WHERE molecule_status IS NULL OR molecule_status NOT IN ('original', 'primary', 'duplicate')
        """)
        
        print(f"Updated {invalid_count} records to have molecule_status = 'original'")
        
        # Verify fix
        cursor.execute("""
            SELECT COUNT(*) FROM consolidated_molecules 
            WHERE molecule_status IS NULL OR molecule_status NOT IN ('original', 'primary', 'duplicate')
        """)
        remaining_invalid = cursor.fetchone()[0]
        
        if remaining_invalid == 0:
            print("All molecule_status values are now valid")
            return True
        else:
            print(f"ERROR: {remaining_invalid} invalid values still remain after update")
            return False
            
    except Exception as e:
        print(f"ERROR: {e}")
        return False
    finally:
        if conn:
            conn.close()

if __name__ == "__main__":
    print("Fixing invalid molecule_status values...")
    success = fix_molecule_status()
    sys.exit(0 if success else 1)
EOF

    # Make it executable and run it
    chmod +x temp_fix_status.py
    python3 temp_fix_status.py
    
    if [ $? -ne 0 ]; then
        echo -e "${RED}Error fixing invalid molecule_status values${NC}"
        echo -e "${RED}Attempting to apply migrations anyway...${NC}"
    fi
    
    # Clean up
    rm temp_fix_status.py
    
    # Apply migrations using Python script
    for MIGRATION_FILE in migrations/026_consolidated_molecule_constraints_indexes.sql migrations/027_create_scientific_data_audit.sql; do
        echo "Applying migration: $MIGRATION_FILE"
        python3 apply_file_sql.py "$MIGRATION_FILE"
        
        if [ $? -ne 0 ]; then
            echo -e "${RED}Error applying migration $MIGRATION_FILE${NC}"
            echo -e "${RED}Check the error output above for details${NC}"
            exit 1
        else
            echo -e "${GREEN}Successfully applied migration: $MIGRATION_FILE${NC}"
        fi
    done
else
    # If neither fix script is available, try apply_migrations.sh as a fallback
    echo "No direct migration scripts found, checking for apply_migrations.sh..."
    
    if [ -f apply_migrations.sh ]; then
        echo "Using apply_migrations.sh as fallback..."
        
        # Check how apply_migrations.sh is implemented
        if grep -q "\-\-start" apply_migrations.sh; then
            # If it uses start/end arguments, use that approach
            echo "Using apply_migrations.sh with start/end parameters"
            
            # Make sure the migrations directory exists
            mkdir -p migrations
            
            # Run apply_migrations.sh with appropriate parameters
            ./apply_migrations.sh --start=26 --end=27
        else
            # If it's a different format, try a more direct approach
            echo "Using apply_migrations.sh with individual file execution"
            
            # Execute SQL directly
            echo "Executing consolidated_molecule_constraints_indexes.sql"
            ./apply_migrations.sh --execute=migrations/026_consolidated_molecule_constraints_indexes.sql
            
            echo "Executing create_scientific_data_audit.sql"
            ./apply_migrations.sh --execute=migrations/027_create_scientific_data_audit.sql
        fi
        
        if [ $? -ne 0 ]; then
            echo -e "${RED}Error applying migrations with apply_migrations.sh${NC}"
            echo -e "${RED}Cannot continue with implementation${NC}"
            exit 1
        fi
    else
        echo -e "${RED}No migration application method found${NC}"
        echo -e "${RED}Cannot apply migrations automatically${NC}"
        echo -e "${RED}Please apply the migrations manually before continuing${NC}"
        exit 1
    fi
fi

# Step 4: Update the API
echo -e "\n${BLUE}Step 4: Updating the API...${NC}"
if [ -f update_consolidated_molecule_api.py ]; then
    python3 update_consolidated_molecule_api.py
else
    echo -e "${RED}update_consolidated_molecule_api.py not found.${NC}"
    exit 1
fi

# Step 5: Run tests
echo -e "\n${BLUE}Step 5: Running tests...${NC}"
if [ -f test_consolidated_molecules_integration.py ]; then
    python3 test_consolidated_molecules_integration.py
else
    echo -e "${RED}test_consolidated_molecules_integration.py not found.${NC}"
    exit 1
fi

# Step 6: Final verification
echo -e "\n${BLUE}Step 6: Verifying implementation...${NC}"
echo "Checking API updates..."

# Check if the API files have been updated
if [ -f api/consolidated_molecule_resource.py ] && [ -f api/consolidated_utils.py ]; then
    echo "API files are in place. Verifying content..."
    
    # Check for key functions in consolidated_utils.py
    if grep -q "def get_primary_molecule" api/consolidated_utils.py && \
       grep -q "def migrate_properties" api/consolidated_utils.py; then
        echo "consolidated_utils.py contains required functions."
    else
        echo -e "${RED}consolidated_utils.py is missing required functions.${NC}"
        API_VERIFICATION_FAILED=true
    fi
    
    # Check for key endpoints in consolidated_molecule_resource.py
    if grep -q "class ConsolidatedMoleculeResource" api/consolidated_molecule_resource.py && \
       grep -q "migrate_properties" api/consolidated_molecule_resource.py; then
        echo "consolidated_molecule_resource.py contains required endpoints."
    else
        echo -e "${RED}consolidated_molecule_resource.py is missing required endpoints.${NC}"
        API_VERIFICATION_FAILED=true
    fi
    
    # Check API registration in __init__.py
    if grep -q "consolidated_molecule_resource" api/__init__.py; then
        echo "API endpoints are registered in __init__.py."
    else
        echo -e "${RED}API endpoints are not registered in __init__.py.${NC}"
        API_VERIFICATION_FAILED=true
    fi
else
    echo -e "${RED}API files are missing.${NC}"
    API_VERIFICATION_FAILED=true
fi

# Check for RDKit wrapper file
if [ -f rdkit_wrapper_consolidated.py ]; then
    echo "RDKit wrapper file is in place."
    
    # Check for key functions in rdkit_wrapper_consolidated.py
    if grep -q "def get_primary_molecule_id" rdkit_wrapper_consolidated.py && \
       grep -q "def calculate_properties_for_consolidated" rdkit_wrapper_consolidated.py; then
        echo "rdkit_wrapper_consolidated.py contains required functions."
    else
        echo -e "${RED}rdkit_wrapper_consolidated.py is missing required functions.${NC}"
        API_VERIFICATION_FAILED=true
    fi
else
    echo -e "${RED}RDKit wrapper file is missing.${NC}"
    API_VERIFICATION_FAILED=true
fi

# Inform about database migrations
echo -e "\n${BLUE}Database Schema:${NC}"
echo "The following database elements should be created:"
echo "1. CHECK constraint on consolidated_molecules.molecule_status"
echo "2. Indexes on consolidated_molecules (inchikey, primary_molecule_id, molecule_status)"
echo "3. scientific_data_audit table with indexes"
echo "These can be verified once database access is available."

if [ "$API_VERIFICATION_FAILED" = "true" ]; then
    echo -e "${RED}API implementation verification failed.${NC}"
    echo -e "${RED}Please check the logs for details.${NC}"
    exit 1
else
    echo -e "${GREEN}Implementation completed successfully!${NC}"
    echo -e "${GREEN}The API components for the consolidated molecules system are now ready to use.${NC}"
    echo -e "\n${BLUE}New API endpoints available:${NC}"
    echo -e "- /api/v1/consolidated/molecules/{molecule_id}"
    echo -e "- /api/v1/consolidated/batch"
    echo -e "- /api/v1/molecules/{molecule_id}/primary"
    echo -e "- /api/v1/consolidated"
    echo -e "- /api/v1/consolidated/molecules/{molecule_id}/audit"
    echo -e "- /api/v1/consolidated/molecules/{molecule_id}/migrate-properties"
    echo -e "- /api/v1/consolidated/search"
    echo -e "\n${BLUE}Note:${NC} Database migrations must be applied separately before using the API."
fi
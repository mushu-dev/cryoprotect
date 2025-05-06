#!/bin/bash
# CryoProtect v2 - Linux Environment Setup Helper
# This script helps set up the Linux environment for CryoProtect migration

echo "Creating Linux environment for CryoProtect..."

# 1. Ensure .env file exists
if [ ! -f .env ]; then
    echo "Creating .env file from template..."
    cp .env.template .env
    echo "✅ Created .env file from template. Please edit it with your credentials."
else
    echo "✅ .env file already exists."
fi

# 2. Update Docker secrets path for Linux
DOCKER_SECRETS_PATH="/run/secrets"
echo "Setting Docker secrets path to $DOCKER_SECRETS_PATH"
# Check if directory exists, if not create it for testing
if [ ! -d "$DOCKER_SECRETS_PATH" ]; then
    echo "Creating Docker secrets directory for testing..."
    sudo mkdir -p "$DOCKER_SECRETS_PATH"
    sudo chmod 755 "$DOCKER_SECRETS_PATH"
fi

# 3. Fix file permissions for all scripts
echo "Fixing permissions for all shell scripts..."
find . -name "*.sh" -exec chmod +x {} \;
echo "✅ Made all shell scripts executable."

# 4. Create local database if not exists
echo "Checking PostgreSQL installation..."
if command -v psql >/dev/null 2>&1; then
    echo "✅ PostgreSQL is installed."
    
    # Try to connect to PostgreSQL
    if psql -U postgres -c '\l' >/dev/null 2>&1; then
        echo "✅ Connected to PostgreSQL server."
        
        # Check if cryoprotect database exists
        if ! psql -U postgres -lqt | cut -d \| -f 1 | grep -qw cryoprotect; then
            echo "Creating cryoprotect database..."
            psql -U postgres -c "CREATE DATABASE cryoprotect;"
            echo "✅ Created cryoprotect database."
        else
            echo "✅ Cryoprotect database already exists."
        fi
    else
        echo "⚠️ Could not connect to PostgreSQL server. Please check your PostgreSQL installation."
    fi
else
    echo "⚠️ PostgreSQL is not installed. Please install PostgreSQL."
    echo "  You can install it with: sudo apt install postgresql postgresql-contrib"
fi

# 5. Update environment setup instructions
echo "
✅ Linux environment setup complete!

Next steps:
1. Edit the .env file with your Supabase credentials
2. Set up conda environment with: ./setup_environment.sh
3. Apply database migrations with: node migrations/apply_migration.js
4. Run the application with: ./run_app_with_fix.sh

For troubleshooting:
- Check PostgreSQL status with: systemctl status postgresql
- Check RDKit installation with: conda list rdkit
- Check database connection with: python check_supabase_connection.py
"

# Make this script executable
chmod +x "$(realpath "$0")"
#!/bin/bash

echo "==============================================================="
echo "CryoProtect v2 - Setup Database Remediation Environment"
echo "==============================================================="
echo

# Check if Python is installed
if ! command -v python3 &> /dev/null; then
    echo "Error: Python is not installed or not in the PATH."
    echo "Please install Python and try again."
    exit 1
fi

echo "This script will set up the environment for the database remediation process."
echo "It will:"
echo "1. Install required Python packages"
echo "2. Create a .env file for your Supabase credentials"
echo

read -p "Do you want to continue? (y/n): " CONFIRM
if [[ $CONFIRM != "y" && $CONFIRM != "Y" ]]; then
    echo "Operation cancelled."
    exit 0
fi

echo
echo "Installing required Python packages..."
pip3 install supabase python-dotenv

if [ $? -ne 0 ]; then
    echo "Error: Failed to install required packages."
    echo "Please check your internet connection and try again."
    exit 1
fi

echo
echo "Creating .env file..."

if [ -f .env ]; then
    echo ".env file already exists."
    read -p "Do you want to overwrite it? (y/n): " OVERWRITE
    if [[ $OVERWRITE != "y" && $OVERWRITE != "Y" ]]; then
        echo "Keeping existing .env file."
        SKIP_ENV=true
    fi
fi

if [ "$SKIP_ENV" != "true" ]; then
    echo "# CryoProtect v2 - Database Remediation Environment Variables" > .env
    echo >> .env

    echo "Please enter your Supabase credentials:"
    echo

    read -p "Supabase URL (e.g., https://abcdefghijklm.supabase.co): " SUPABASE_URL
    read -p "Supabase Service Role Key: " SUPABASE_KEY
    read -p "Supabase Project ID (optional): " SUPABASE_PROJECT_ID

    echo >> .env
    echo "# Supabase Project URL" >> .env
    echo "SUPABASE_URL=$SUPABASE_URL" >> .env
    echo >> .env
    echo "# Supabase Service Role Key" >> .env
    echo "SUPABASE_KEY=$SUPABASE_KEY" >> .env

    if [ ! -z "$SUPABASE_PROJECT_ID" ]; then
        echo >> .env
        echo "# Supabase Project ID" >> .env
        echo "SUPABASE_PROJECT_ID=$SUPABASE_PROJECT_ID" >> .env
    fi

    echo ".env file created successfully."
fi

echo
echo "Environment setup completed successfully."
echo "You can now run the database remediation process using:"
echo "  ./run_database_remediation.sh"
echo

read -p "Press Enter to exit..."

# Make the remediation scripts executable
chmod +x run_database_remediation.sh
chmod +x verify_database_remediation.sh
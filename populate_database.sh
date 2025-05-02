#!/bin/bash

echo "CryoProtect v2 - Database Population Script"
echo "=========================================="
echo

# Check if Python is installed
if ! command -v python3 &> /dev/null; then
    echo "Error: Python 3 is not installed or not in PATH."
    echo "Please install Python 3.8 or higher and try again."
    exit 1
fi

# Check if .env file exists
if [ ! -f .env ]; then
    echo "Warning: .env file not found."
    echo "Creating a template .env file. Please edit it with your Supabase credentials."
    cat > .env << EOL
SUPABASE_URL=your_supabase_url
SUPABASE_KEY=your_supabase_key
SUPABASE_USER=your_supabase_user_email
SUPABASE_PASSWORD=your_supabase_user_password
EOL
    echo
    echo "Template .env file created. Please edit it and run this script again."
    exit 1
fi

# Check if required packages are installed
echo "Checking required packages..."
if ! python3 -c "import dotenv, supabase" &> /dev/null; then
    echo "Installing required packages..."
    pip3 install python-dotenv supabase
fi

# Make the script executable
chmod +x populate_database_supabase.py
chmod +x verify_database_population.py

# Run the population script
echo
echo "Running database population script..."
echo
python3 ./populate_database_supabase.py "$@"

# Check if the script ran successfully
if [ $? -ne 0 ]; then
    echo
    echo "Error: Database population failed."
    exit 1
fi

echo
echo "Database population completed successfully."
echo
echo "Running verification script..."
echo
python3 ./verify_database_population.py

echo
echo "Process completed."
echo "See database_population.log and database_verification.log for details."
echo
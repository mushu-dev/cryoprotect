#!/bin/bash

echo "CryoProtect v2 - Create exec_sql Function"
echo "========================================"
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

# Create a temporary Python script to execute the SQL
cat > create_function_temp.py << EOL
import os
from dotenv import load_dotenv
from supabase import create_client

load_dotenv()

SUPABASE_URL = os.getenv("SUPABASE_URL")
SUPABASE_KEY = os.getenv("SUPABASE_KEY")

if not SUPABASE_URL or not SUPABASE_KEY:
    raise ValueError("SUPABASE_URL and SUPABASE_KEY must be set in .env file")

supabase = create_client(SUPABASE_URL, SUPABASE_KEY)

with open("create_exec_sql_function.sql", "r") as f:
    sql = f.read()

print("Executing SQL to create exec_sql function...")
response = supabase.rpc("exec_sql", {"query": sql}).execute()

print("Function created successfully!")
print("You can now use the exec_sql function in your verification scripts.")
EOL

# Make the script executable
chmod +x create_function_temp.py

# Run the temporary Python script
echo
echo "Creating exec_sql function in Supabase..."
python3 ./create_function_temp.py

# Check if the script ran successfully
if [ $? -ne 0 ]; then
    echo
    echo "Error: Failed to create exec_sql function."
    echo "Please check your Supabase credentials and try again."
    echo "You can also create the function manually using the SQL in create_exec_sql_function.sql."
    rm create_function_temp.py
    exit 1
fi

# Clean up
rm create_function_temp.py

echo
echo "exec_sql function created successfully!"
echo "You can now run the verification script."
echo
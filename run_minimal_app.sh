#\!/bin/bash

# Activate virtual environment
source quick_env/bin/activate

# Set Supabase environment variables from .env file if it exists
if [ -f .env ]; then
    echo "Loading environment variables from .env file"
    export $(grep -v '^#' .env  < /dev/null |  xargs)
fi

# Run the minimal app
echo "Starting minimal app..."
python minimal_app.py

# Deactivate virtual environment on exit
deactivate

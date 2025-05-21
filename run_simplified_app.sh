#\!/bin/bash

# Activate virtual environment
source quick_env/bin/activate

# Set environment variables from .env file if it exists
if [ -f .env ]; then
    echo "Loading environment variables from .env file"
    export $(grep -v '^#' .env  < /dev/null |  xargs)
fi

# Run the simplified app
echo "Starting simplified app on port 5000..."
export PORT=5000
python simplified_app.py

# Deactivate virtual environment on exit
deactivate

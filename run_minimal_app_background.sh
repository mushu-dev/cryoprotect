#\!/bin/bash

# Activate virtual environment
source quick_env/bin/activate

# Set Supabase environment variables from .env file if it exists
if [ -f .env ]; then
    echo "Loading environment variables from .env file"
    export $(grep -v '^#' .env  < /dev/null |  xargs)
fi

# Kill any existing instances of the app
pkill -f "python minimal_app.py" || echo "No existing instances found"

# Run the minimal app in the background
echo "Starting minimal app in background..."
nohup python minimal_app.py > minimal_app.log 2>&1 &
APP_PID=$!

echo "App started with PID $APP_PID"
echo "Logs are being written to minimal_app.log"
echo "Wait 2 seconds for app to start..."
sleep 2

# Check if the app is still running
if ps -p $APP_PID > /dev/null; then
    echo "App is running. You can access it at http://localhost:5000"
    echo "You can stop it later with: kill $APP_PID"
else
    echo "App failed to start. Check minimal_app.log for details."
    exit 1
fi

# Deactivate virtual environment
deactivate

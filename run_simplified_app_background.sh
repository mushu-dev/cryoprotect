#\!/bin/bash

# Activate virtual environment
cd /home/mushu/Projects/CryoProtect
source quick_env/bin/activate

# Set environment variables from .env file if it exists
if [ -f .env ]; then
    echo "Loading environment variables from .env file"
    export $(grep -v '^#' .env  < /dev/null |  xargs)
fi

# Kill any existing instances
pkill -f "python simplified_app.py" || echo "No existing instances found"

# Run the simplified app in the background
echo "Starting simplified app in background..."
nohup python simplified_app.py > simplified_app.log 2>&1 &
APP_PID=$!

echo "App started with PID $APP_PID"
echo "Logs are being written to simplified_app.log"
echo "Wait 2 seconds for app to start..."
sleep 2

# Check if the app is running
if ps -p $APP_PID > /dev/null; then
    echo "App is running successfully!"
else
    echo "App failed to start. Check simplified_app.log for details."
    exit 1
fi

# Display log
echo "First few lines of log:"
head -10 simplified_app.log

# Deactivate virtual environment
deactivate


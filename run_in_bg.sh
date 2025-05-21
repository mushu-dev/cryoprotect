#\!/bin/bash

# Activate virtual environment
cd /home/mushu/Projects/CryoProtect
source quick_env/bin/activate

# Kill any existing app
pkill -f "python minimal_app.py" || echo "No existing instances found"

# Start in background
nohup python minimal_app.py > minimal_app.log 2>&1 &
APP_PID=$\!
echo "Started app with PID $APP_PID"

# Wait for app to start
sleep 2

# Check dependencies
curl http://localhost:5000/dependencies

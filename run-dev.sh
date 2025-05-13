#!/bin/bash
# Script to run both frontend and backend in development mode

# Set some colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

echo -e "${YELLOW}CryoProtect Development Server${NC}"

# Check for running processes
check_port() {
    local port=$1
    if lsof -Pi :$port -sTCP:LISTEN -t >/dev/null ; then
        return 0  # Port is in use
    else
        return 1  # Port is free
    fi
}

# Check if ports are already in use
if check_port 3000; then
    echo -e "${RED}Error: Port 3000 is already in use by another process${NC}"
    echo "Please stop that process first or use a different port"
    exit 1
fi

if check_port 5000; then
    echo -e "${RED}Error: Port 5000 is already in use by another process${NC}"
    echo "Please stop that process first or use a different port"
    exit 1
fi

# Start backend in a separate terminal/process
start_backend() {
    echo -e "${GREEN}Starting backend server on port 5000...${NC}"
    cd "$PWD"
    
    # Use the appropriate command to start the Flask app
    if [ -f "run_app.sh" ]; then
        ./run_app.sh &
    else
        # Default Flask command
        FLASK_APP=app.py FLASK_ENV=development flask run --host=0.0.0.0 --port=5000 &
    fi
    
    # Store the PID
    BACKEND_PID=$!
    echo "Backend server started with PID: $BACKEND_PID"
}

# Start frontend in a separate terminal/process
start_frontend() {
    echo -e "${GREEN}Starting frontend server on port 3000...${NC}"
    cd "$PWD/frontend"
    
    # Start the development server
    npm run dev &
    
    # Store the PID
    FRONTEND_PID=$!
    echo "Frontend server started with PID: $FRONTEND_PID"
}

# Setup signal handling
cleanup() {
    echo -e "${YELLOW}Stopping servers...${NC}"
    
    # Kill the frontend and backend processes if they exist
    if [ ! -z "$FRONTEND_PID" ]; then
        echo "Stopping frontend server (PID: $FRONTEND_PID)"
        kill -TERM $FRONTEND_PID 2>/dev/null || true
    fi
    
    if [ ! -z "$BACKEND_PID" ]; then
        echo "Stopping backend server (PID: $BACKEND_PID)"
        kill -TERM $BACKEND_PID 2>/dev/null || true
    fi
    
    echo -e "${GREEN}All servers stopped${NC}"
    exit 0
}

# Setup signal handlers
trap cleanup SIGINT SIGTERM

# Start the servers
start_backend
start_frontend

# Show URLs
echo
echo -e "${GREEN}Development servers started:${NC}"
echo -e "Frontend: ${YELLOW}http://localhost:3000${NC}"
echo -e "Backend API: ${YELLOW}http://localhost:5000${NC}"
echo
echo -e "Press ${YELLOW}Ctrl+C${NC} to stop both servers"

# Keep the script running to handle signals
while true; do
    sleep 1
done
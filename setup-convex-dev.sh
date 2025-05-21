#!/bin/bash
# setup-convex-dev.sh
#
# Setup script for CryoProtect development with Convex integration.
# This script configures the environment variables and starts the 
# appropriate development server.

# Display usage information
usage() {
  echo "Usage: $0 [main|minimal]"
  echo "  main     - Setup and run the main frontend with Convex"
  echo "  minimal  - Setup and run the minimal frontend with Convex"
  exit 1
}

# Check for environment argument
if [ $# -ne 1 ]; then
  usage
fi

# Configure based on argument
case "$1" in
  "main")
    # Configure for main frontend
    export USE_CONVEX=true
    export CONVEX_URL=https://hallowed-malamute-424.convex.cloud
    export NEXT_PUBLIC_USE_CONVEX=true
    export NEXT_PUBLIC_CONVEX_URL=https://hallowed-malamute-424.convex.cloud
    
    echo "Setting up main frontend with Convex integration"
    echo "NEXT_PUBLIC_USE_CONVEX=$NEXT_PUBLIC_USE_CONVEX"
    echo "NEXT_PUBLIC_CONVEX_URL=$NEXT_PUBLIC_CONVEX_URL"
    
    # Change to frontend directory
    cd frontend
    ;;
    
  "minimal")
    # Configure for minimal frontend
    export USE_CONVEX=true
    export CONVEX_URL=https://hallowed-malamute-424.convex.cloud
    export NEXT_PUBLIC_USE_CONVEX=true
    export NEXT_PUBLIC_CONVEX_URL=https://hallowed-malamute-424.convex.cloud
    
    echo "Setting up minimal frontend with Convex integration"
    echo "NEXT_PUBLIC_USE_CONVEX=$NEXT_PUBLIC_USE_CONVEX"
    echo "NEXT_PUBLIC_CONVEX_URL=$NEXT_PUBLIC_CONVEX_URL"
    
    # Change to minimal-frontend directory
    cd minimal-frontend
    ;;
    
  *)
    # Invalid argument
    usage
    ;;
esac

# Check if we're in the right directory
if [ ! -f "package.json" ]; then
  echo "Error: Could not find package.json in $(pwd)"
  exit 1
fi

# Run development server
echo "Starting development server..."
npm run dev
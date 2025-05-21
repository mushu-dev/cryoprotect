#!/bin/bash
# Simple Netlify build script that doesn't rely on submodules
set -e

echo "Running Netlify build script..."

# Check if we're in a Netlify environment
if [ -n "$NETLIFY" ]; then
  echo "Running in Netlify environment"
  
  # Ensure we ignore any submodule errors
  if [ -f ".gitmodules" ]; then
    echo "Found .gitmodules file, but we're configured to ignore submodules"
  fi
  
  # If the frontend directory exists, navigate to it
  if [ -d "frontend" ]; then
    cd frontend
    echo "Switched to frontend directory"
    
    # Install dependencies
    echo "Installing dependencies..."
    npm ci || npm install
    
    # Build the project
    echo "Building project..."
    npm run build
    
    echo "Build completed successfully!"
  else
    echo "Frontend directory not found, using pre-built files"
    
    # If using pre-built files, no build needed
    echo "Using pre-built files, no build needed"
  fi
else
  echo "Running in local environment"
  
  # Local build process (optional)
  if [ -d "frontend" ]; then
    cd frontend
    npm install
    npm run build
  fi
fi

echo "Netlify build script completed."
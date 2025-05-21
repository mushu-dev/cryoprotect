#!/bin/bash

# Test backend connectivity for CryoProtect services
# This script runs the test-connection.js Node.js script

# Check if Node.js is installed
if ! command -v node &> /dev/null; then
    echo "Node.js not found. Please install it first."
    exit 1
fi

# Check if axios is installed
if [ ! -d "node_modules/axios" ]; then
    echo "Installing axios..."
    npm install axios --no-save
fi

# Get custom service URLs from arguments if provided
FRONTEND_URL=${1:-"https://cryoprotect.netlify.app"}
API_URL=${2:-"https://cryoprotect-8030e4025428.herokuapp.com"}
RDKIT_URL=${3:-"https://rdkit.cryoprotect.app"}
CONVEX_URL=${4:-"https://dynamic-mink-63.convex.cloud"}

# Show information about the test
echo "Running backend connectivity tests for CryoProtect services:"
echo " Frontend URL: $FRONTEND_URL"
echo " API URL: $API_URL"
echo " RDKit URL: $RDKIT_URL"
echo " Convex URL: $CONVEX_URL"
echo

# Export the URLs as environment variables for the Node.js script
export FRONTEND_URL
export API_URL
export RDKIT_URL
export CONVEX_URL

# Run the test script
node test-connection.js

# Check exit code
if [ $? -eq 0 ]; then
    echo "Connection test completed successfully."
else
    echo "Connection test failed."
fi
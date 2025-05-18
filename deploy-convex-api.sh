#!/bin/bash
# Script to deploy Convex API functions

set -e  # Exit on any error

echo "Deploying Convex API functions..."

# Check if npx is available
if ! command -v npx &> /dev/null; then
    echo "Error: npx could not be found. Please install Node.js and npm."
    exit 1
fi

# Check if we're in the right directory (should have a convex folder)
if [ ! -d "convex" ]; then
    echo "Error: convex directory not found. Please run this script from the project root."
    exit 1
fi

# Deploy to Convex
echo "Deploying to Convex..."
npx convex deploy

echo "Deployment complete!"
echo "Testing Convex API endpoints..."

# Test the API endpoints
echo "Testing query endpoint..."
curl -X POST "https://dynamic-mink-63.convex.cloud/api/query" \
  -H "Content-Type: application/json" \
  --data '{"table":"molecules","limit":1}' \
  --fail || echo "Query endpoint test failed"

echo "Testing insert endpoint..."
curl -X POST "https://dynamic-mink-63.convex.cloud/api/insert" \
  -H "Content-Type: application/json" \
  --data '{"table":"test","data":{"name":"test"}}' \
  --fail || echo "Insert endpoint test failed"

echo "Testing update endpoint..."
curl -X POST "https://dynamic-mink-63.convex.cloud/api/update" \
  -H "Content-Type: application/json" \
  --data '{"table":"test","data":{"name":"updated"},"filters":{"name":"test"}}' \
  --fail || echo "Update endpoint test failed"

echo "Testing delete endpoint..."
curl -X POST "https://dynamic-mink-63.convex.cloud/api/delete" \
  -H "Content-Type: application/json" \
  --data '{"table":"test","filters":{"name":"updated"}}' \
  --fail || echo "Delete endpoint test failed"

echo "All tests complete!"
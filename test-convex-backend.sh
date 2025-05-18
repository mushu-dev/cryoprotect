#!/bin/bash
# Test the backend integration with Convex

echo "Testing backend integration with Convex..."

# Test the API root endpoint
echo "Testing API root endpoint..."
curl -s https://cryoprotect-8030e4025428.herokuapp.com/

# Test fetching molecules
echo
echo "Testing molecules endpoint..."
curl -s https://cryoprotect-8030e4025428.herokuapp.com/api/molecules?limit=2

# Get the database type
echo
echo "Testing database configuration..."
curl -s https://cryoprotect-8030e4025428.herokuapp.com/api/config | grep database_type

echo
echo "Backend integration test complete!"
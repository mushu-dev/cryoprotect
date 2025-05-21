#!/bin/bash

# API connectivity test script
BACKEND_URL="https://cryoprotect-8030e4025428.herokuapp.com"
FRONTEND_ORIGIN="https://frontend-cryoprotect.vercel.app"

echo "Testing API connectivity with backend: $BACKEND_URL"
echo "Using frontend origin: $FRONTEND_ORIGIN"
echo

# Test health endpoint
echo "Testing health endpoint..."
curl -s "$BACKEND_URL/health" | jq .
echo

# Test API connect endpoint
echo "Testing API connect endpoint..."
curl -s -H "Origin: $FRONTEND_ORIGIN" "$BACKEND_URL/api/connect" | jq .
echo

# Test direct molecule endpoint 
echo "Testing molecule endpoint..."
curl -s -H "Origin: $FRONTEND_ORIGIN" "$BACKEND_URL/api/molecules?limit=1" | jq .
echo

echo "All tests completed."
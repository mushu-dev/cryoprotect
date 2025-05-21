#!/bin/bash
# Script to test connectivity between Heroku backend and Vercel frontend

# Default values
HEROKU_APP_NAME=${HEROKU_APP_NAME:-cryoprotect}
VERCEL_FRONTEND_URL=${VERCEL_FRONTEND_URL:-https://frontend-cryoprotect.vercel.app}

# Construct the backend API URL
BACKEND_URL="https://cryoprotect-8030e4025428.herokuapp.com"

echo "Testing connectivity between:"
echo "- Backend: $BACKEND_URL"
echo "- Frontend: $VERCEL_FRONTEND_URL"
echo ""

# Test 1: Backend Health Check
echo "Test 1: Backend Health Check"
if curl -s "$BACKEND_URL/health" | grep -q "status"; then
    echo "✅ Backend health check successful"
else
    echo "❌ Backend health check failed"
fi

# Test 2: API Connectivity Endpoint
echo ""
echo "Test 2: API Connectivity Endpoint"
if curl -s "$BACKEND_URL/api/v1/health/connectivity" | grep -q "connected"; then
    echo "✅ API connectivity endpoint successful"
else
    echo "❌ API connectivity endpoint failed"
fi

# Test 3: CORS Headers Check
echo ""
echo "Test 3: CORS Headers Check"
CORS_HEADERS=$(curl -s -I -X OPTIONS "$BACKEND_URL/api/v1/health/connectivity" -H "Origin: $VERCEL_FRONTEND_URL")
if echo "$CORS_HEADERS" | grep -q "Access-Control-Allow-Origin"; then
    echo "✅ CORS headers are properly configured"
    echo "CORS Headers:"
    echo "$CORS_HEADERS" | grep "Access-Control-"
else
    echo "❌ CORS headers are missing"
fi

# Test 4: Frontend to Backend API Request
echo ""
echo "Test 4: Frontend to Backend API Request"
echo "This test requires the frontend to be deployed."
echo "Please visit $VERCEL_FRONTEND_URL in your browser and check the network tab"
echo "for requests to $BACKEND_URL to verify connectivity."

echo ""
echo "------------------------------------------------------------"
echo "If any tests failed, check the following:"
echo "1. Both applications are deployed and running"
echo "2. CORS settings are configured correctly"
echo "3. Environment variables are set properly"
echo "4. Network rules allow communication between services"
echo "------------------------------------------------------------"
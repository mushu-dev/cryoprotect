#!/bin/bash
# Script to deploy to Vercel with protection bypass but without Vercel Authentication

# Generate a random NEXTAUTH_SECRET if not provided
NEXTAUTH_SECRET=${NEXTAUTH_SECRET:-$(openssl rand -base64 32)}

# Get Heroku app name from environment or default
HEROKU_APP_NAME=${HEROKU_APP_NAME:-cryoprotect}

# Construct the backend API URL
API_URL="https://cryoprotect-8030e4025428.herokuapp.com/v1"

# Protection Bypass Token
PROTECTION_BYPASS=${PROTECTION_BYPASS:-TAt23KbtFE8dkZobJU3hpgTP4L5ja07V}

echo "Deploying to Vercel with the following configuration:"
echo "Backend API URL: $API_URL"
echo "Protection Bypass: Configured"

# Deploy to Vercel with updated environment variables
# Use --no-clipboard to prevent copying the deployment URL to clipboard
vercel deploy --prod --yes \
  -e NEXT_PUBLIC_API_URL=$API_URL \
  -e NEXT_PUBLIC_USE_MOCK_DATA=false \
  -e NEXT_PUBLIC_ENABLE_API_LOGGING=true \
  -e NEXT_PUBLIC_ENVIRONMENT=production \
  -e NEXT_PUBLIC_FRONTEND_PROTECTION_BYPASS="$PROTECTION_BYPASS" \
  -e NEXTAUTH_SECRET="$NEXTAUTH_SECRET" \
  -e PROTECTION_BYPASS="$PROTECTION_BYPASS"

# Store the deployment URL for testing
DEPLOYMENT_URL=$(vercel ls | head -n 1)
echo "Deployment URL: $DEPLOYMENT_URL"
echo "Testing protection bypass..."

echo "1. Testing without bypass token (should fail):"
curl -s "$DEPLOYMENT_URL" -o /dev/null -w "Status: %{http_code}\n"

echo "2. Testing with bypass token in header (should succeed):"
curl -s -H "x-protection-bypass: $PROTECTION_BYPASS" "$DEPLOYMENT_URL" -o /dev/null -w "Status: %{http_code}\n"

echo "3. Testing with bypass token in query parameter (should succeed):"
curl -s "$DEPLOYMENT_URL?bypass=$PROTECTION_BYPASS" -o /dev/null -w "Status: %{http_code}\n"

echo "4. Testing API health endpoint with bypass token:"
curl -s -H "x-protection-bypass: $PROTECTION_BYPASS" "$DEPLOYMENT_URL/api/v1/health"
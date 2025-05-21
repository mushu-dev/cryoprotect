#!/bin/bash
# Script to deploy the frontend to Vercel with environment variables

# Generate a random NEXTAUTH_SECRET if not provided
NEXTAUTH_SECRET=${NEXTAUTH_SECRET:-$(openssl rand -base64 32)}

# Get Heroku app name from environment or default
HEROKU_APP_NAME=${HEROKU_APP_NAME:-cryoprotect}

# Construct the backend API URL
API_URL="https://cryoprotect-8030e4025428.herokuapp.com/v1"

# Get Vercel frontend URL from environment or default
VERCEL_URL=${VERCEL_URL:-frontend-cryoprotect.vercel.app}
FRONTEND_URL="https://${VERCEL_URL}"

# Get protection bypass token from environment or use default
PROTECTION_BYPASS=${PROTECTION_BYPASS:-TAt23KbtFE8dkZobJU3hpgTP4L5ja07V}

echo "Deploying frontend to Vercel with the following configuration:"
echo "Backend API URL: $API_URL"
echo "Frontend URL: $FRONTEND_URL"
echo "Protection Bypass: Configured"

# Deploy to Vercel with updated environment variables
vercel deploy --archive=tgz --prod --yes \
  -e NEXT_PUBLIC_API_URL=$API_URL \
  -e NEXT_PUBLIC_USE_MOCK_DATA=false \
  -e NEXT_PUBLIC_ENABLE_API_LOGGING=true \
  -e NEXT_PUBLIC_ENVIRONMENT=production \
  -e NEXT_PUBLIC_FRONTEND_PROTECTION_BYPASS="$PROTECTION_BYPASS" \
  -e NEXTAUTH_URL=$FRONTEND_URL \
  -e NEXTAUTH_SECRET="$NEXTAUTH_SECRET" \
  -e PROTECTION_BYPASS="$PROTECTION_BYPASS" \
  -e VERCEL_ANALYTICS_ID=true \
  -e VERCEL_SPEED_INSIGHTS=true

# Once deployment is complete, notify the Heroku backend about the frontend URL
# This helps with CORS configuration
if [ $? -eq 0 ]; then
  echo "Vercel deployment successful!"
  echo "Notifying Heroku backend about frontend URL..."
  
  # Only run this if HEROKU_API_KEY is set
  if [ -n "$HEROKU_API_KEY" ]; then
    curl -X PATCH \
      -H "Content-Type: application/json" \
      -H "Accept: application/vnd.heroku+json; version=3" \
      -H "Authorization: Bearer $HEROKU_API_KEY" \
      -d "{\"config\": {\"VERCEL_FRONTEND_URL\": \"$FRONTEND_URL\"}}" \
      "https://api.heroku.com/apps/${HEROKU_APP_NAME}/config-vars"
      
    echo "Heroku configuration updated with frontend URL."
  else
    echo "HEROKU_API_KEY not set. Skipping Heroku configuration update."
    echo "Please manually set VERCEL_FRONTEND_URL config var on Heroku to: $FRONTEND_URL"
  fi
else
  echo "Vercel deployment failed."
fi
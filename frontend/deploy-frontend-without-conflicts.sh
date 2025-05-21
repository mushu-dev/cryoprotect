#!/bin/bash
# Script to deploy only the frontend with explicit API exclusion

# Change to the frontend directory
cd "$(dirname "$0")"

# Generate a random NEXTAUTH_SECRET if not provided
NEXTAUTH_SECRET=${NEXTAUTH_SECRET:-$(openssl rand -base64 32)}

# Get environment variables
HEROKU_APP_NAME=${HEROKU_APP_NAME:-cryoprotect}
API_URL="https://cryoprotect-8030e4025428.herokuapp.com/v1"
VERCEL_URL=${VERCEL_URL:-frontend-cryoprotect.vercel.app}
FRONTEND_URL="https://${VERCEL_URL}"
PROTECTION_BYPASS=${PROTECTION_BYPASS:-TAt23KbtFE8dkZobJU3hpgTP4L5ja07V}

echo "Creating temp directory for clean deployment..."
TEMP_DIR=$(mktemp -d)
echo "Temporary directory: $TEMP_DIR"

# Copy frontend files only to temp directory
echo "Copying frontend files to temp directory..."
cp -r ./* $TEMP_DIR/
cp .vercelignore $TEMP_DIR/
cp vercel.json $TEMP_DIR/

# Move into the temp directory for deployment
cd $TEMP_DIR

echo "Deploying frontend to Vercel with the following configuration:"
echo "Backend API URL: $API_URL"
echo "Frontend URL: $FRONTEND_URL"
echo "Protection Bypass: Configured"

# Deploy frontend to Vercel with updated environment variables
vercel --archive=tgz --prod --yes \
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

# Save the result
DEPLOY_RESULT=$?

# Clean up the temp directory
echo "Cleaning up temporary directory..."
cd -
rm -rf $TEMP_DIR

# Check if the deployment was successful
if [ $DEPLOY_RESULT -eq 0 ]; then
  echo "Vercel deployment successful!"
  echo "Your site should be live at: $FRONTEND_URL"
  
  # Notify Heroku backend about the frontend URL if HEROKU_API_KEY is set
  if [ -n "$HEROKU_API_KEY" ]; then
    echo "Notifying Heroku backend about frontend URL..."
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
  echo "Vercel deployment failed. Please check the error messages above."
  exit 1
fi
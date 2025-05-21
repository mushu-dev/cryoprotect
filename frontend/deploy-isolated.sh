#!/bin/bash
# Script to deploy frontend to Vercel from a completely isolated directory
# This avoids all path conflicts by creating a clean deployment separate from the repository

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

# Create temp directory outside the repository
TEMP_DIR=$(mktemp -d)
echo "Created temporary directory: $TEMP_DIR"

# Create a fresh vercel project in the temp directory
mkdir -p $TEMP_DIR/pages $TEMP_DIR/public $TEMP_DIR/styles

# Copy only essential Next.js files
echo "Copying essential Next.js files to temp directory..."
cp -r ./src $TEMP_DIR/
cp -r ./public $TEMP_DIR/
cp ./next.config.js $TEMP_DIR/
cp ./package.json $TEMP_DIR/
cp ./package-lock.json $TEMP_DIR/
cp ./postcss.config.js $TEMP_DIR/
cp ./tailwind.config.js $TEMP_DIR/
cp ./tsconfig.json $TEMP_DIR/

# Create a minimal vercel.json
cat > $TEMP_DIR/vercel.json << EOL
{
  "version": 2,
  "buildCommand": "npm install && npm run build",
  "outputDirectory": ".next",
  "framework": "nextjs",
  "regions": ["iad1"],
  "env": {
    "NEXT_PUBLIC_API_URL": "$API_URL",
    "NEXT_PUBLIC_USE_MOCK_DATA": "false",
    "NEXT_PUBLIC_ENABLE_API_LOGGING": "true",
    "NEXT_PUBLIC_ENVIRONMENT": "production",
    "NEXT_PUBLIC_FRONTEND_PROTECTION_BYPASS": "$PROTECTION_BYPASS",
    "NEXTAUTH_URL": "$FRONTEND_URL",
    "NEXTAUTH_SECRET": "$NEXTAUTH_SECRET",
    "PROTECTION_BYPASS": "$PROTECTION_BYPASS",
    "VERCEL_ANALYTICS_ID": "true",
    "VERCEL_SPEED_INSIGHTS": "true"
  }
}
EOL

echo "Moving to temporary directory to deploy from..."
cd $TEMP_DIR

# Initialize vercel in the temp directory
echo "Initializing Vercel project..."
echo 'y' | vercel link --yes

echo "Deploying isolated frontend to Vercel"
echo "Backend API URL: $API_URL"
echo "Frontend URL: $FRONTEND_URL"

# Deploy to Vercel
vercel --prod --yes

# Save the result
DEPLOY_RESULT=$?

# Clean up the temp directory
echo "Cleaning up temporary directory..."
cd -  # Return to original directory
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
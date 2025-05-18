#!/bin/bash
# Script to deploy Convex integration to Heroku

set -e  # Exit on any error

HEROKU_APP_NAME="cryoprotect"

echo "Deploying Convex integration to Heroku..."

# Check if Heroku CLI is installed
if ! command -v heroku &> /dev/null; then
    echo "Error: Heroku CLI could not be found. Please install it first."
    exit 1
fi

# Check if logged in to Heroku
heroku whoami &> /dev/null || { 
    echo "Not logged in to Heroku. Please run 'heroku login' first."; 
    exit 1; 
}

# First, deploy the Convex API functions
echo "Deploying Convex API functions..."
bash ./deploy-convex-api.sh

# Set Heroku environment variables for Convex
echo "Setting Heroku environment variables for Convex..."
heroku config:set USE_CONVEX=true \
    CONVEX_URL="https://dynamic-mink-63.convex.cloud" \
    CONVEX_DEPLOYMENT_KEY="" \
    --app $HEROKU_APP_NAME

# Deploy the Convex adapter and database factory to Heroku
echo "Deploying to Heroku..."
git add database/convex_adapter.py database/db_factory.py
git commit -m "Add Convex adapter and database factory" || echo "No changes to commit"
git push heroku convex-implementation:master

# Test the integration
echo "Testing the Convex integration..."
python test-convex-adapter.py

echo "Deployment complete!"
echo "The backend has been configured to use Convex instead of Supabase."
echo "API requests will now be routed through the Convex adapter."
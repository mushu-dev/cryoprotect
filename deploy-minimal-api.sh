#!/bin/bash
# Script to deploy minimal Convex API functions and update Heroku

set -e

HEROKU_APP_NAME="cryoprotect"

echo "=== Deploying Minimal Convex API Functions ==="
echo

# Step 1: Deploy Convex API Functions
echo "Step 1: Deploying Convex API Functions..."
cd convex-minimal
./deploy.sh
cd ..

# Step 2: Update Heroku with Convex configuration
echo
echo "Step 2: Updating Heroku with Convex configuration..."
heroku config:set USE_CONVEX=true \
    CONVEX_URL="https://upbeat-parrot-866.convex.cloud" \
    CONVEX_DEPLOYMENT_KEY="" \
    --app $HEROKU_APP_NAME

# Step 3: Deploy Heroku (add the Convex adapter and database factory)
echo
echo "Step 3: Deploying Convex adapter to Heroku..."
git add database/convex_adapter.py database/db_factory.py
git commit -m "Add Convex adapter and database factory" || echo "No changes to commit"
git push heroku convex-implementation:master

# Step 4: Test the Heroku deployment
echo
echo "Step 4: Testing the Heroku deployment..."
curl -s https://$HEROKU_APP_NAME.herokuapp.com/v1/health | grep "status"

echo
echo "=== Deployment Complete ==="
echo
echo "Your Heroku app is now configured to use Convex."
echo "The minimal API functions have been deployed to Convex."
echo "You can now test the integration with:"
echo "  python test-convex-adapter.py"
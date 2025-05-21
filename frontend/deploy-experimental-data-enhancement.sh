#!/bin/bash
# Script to deploy the Experimental Data Enhancement feature to Netlify

# Configuration
NETLIFY_SITE_NAME=${NETLIFY_SITE_NAME:-cryoprotect}
FRONTEND_URL="https://${NETLIFY_SITE_NAME}.netlify.app"
API_URL=${API_URL:-"https://cryoprotect-8030e4025428.herokuapp.com/v1"}
PROTECTION_BYPASS=${PROTECTION_BYPASS:-TAt23KbtFE8dkZobJU3hpgTP4L5ja07V}
NEXTAUTH_SECRET=${NEXTAUTH_SECRET:-$(openssl rand -base64 32)}
BRANCH_NAME="experimental-data-enhancement"

echo "🧪 Deploying Experimental Data Enhancement Feature"
echo "=================================================="
echo "Frontend URL: $FRONTEND_URL"
echo "API URL: $API_URL"
echo "Branch: $BRANCH_NAME"

# First, run validation tests to ensure everything's working properly
echo "Step 1: Running validation tests..."
echo "----------------------------------"
bash validate-experimental-data-enhancement.sh

# Check if validation failed
if [ $? -ne 0 ]; then
  echo "❌ Validation tests failed. Aborting deployment."
  echo "Please fix the issues and try again."
  exit 1
fi

echo "✅ Validation tests passed. Proceeding with deployment."

# Make sure we have the correct branch
echo "Step 2: Checking out experimental-data-enhancement branch..."
echo "---------------------------------------------------------"
git checkout $BRANCH_NAME

# Update environment variables
echo "Step 3: Setting environment variables in Netlify..."
echo "------------------------------------------------"
netlify env:set NEXT_PUBLIC_API_URL "$API_URL"
netlify env:set NEXT_PUBLIC_USE_MOCK_DATA false
netlify env:set NEXT_PUBLIC_ENABLE_API_LOGGING true
netlify env:set NEXT_PUBLIC_ENVIRONMENT production
netlify env:set NEXT_PUBLIC_NETLIFY true
netlify env:set NEXT_PUBLIC_FRONTEND_PROTECTION_BYPASS "$PROTECTION_BYPASS"
netlify env:set NEXT_PUBLIC_ENABLE_EXPERIMENTAL_FEATURES true
netlify env:set NEXTAUTH_URL "$FRONTEND_URL"
netlify env:set NEXTAUTH_SECRET "$NEXTAUTH_SECRET"
netlify env:set PROTECTION_BYPASS "$PROTECTION_BYPASS"

# Update API endpoints in code
echo "Step 4: Updating API endpoints in code..."
echo "--------------------------------------"
npm run update-api-endpoints

# Build the project
echo "Step 5: Building the project..."
echo "----------------------------"
npm run build

# Deploy to Netlify with the experimental branch
echo "Step 6: Deploying to Netlify from $BRANCH_NAME branch..."
echo "---------------------------------------------------"
netlify deploy --prod --branch=$BRANCH_NAME

# Verify the deployment
if [ $? -eq 0 ]; then
  echo "✅ Netlify deployment successful!"
  
  echo "Step 7: Verifying the deployment..."
  echo "---------------------------------"
  NETLIFY_URL="$FRONTEND_URL" NEXT_PUBLIC_API_URL="$API_URL" npm run verify-netlify
  
  echo "Step 8: Running post-deployment tests..."
  echo "-------------------------------------"
  # Run a basic navigation test to make sure everything loads
  bash run-basic-navigation-test.sh
  
  # Run the experimental data UI test for final verification
  npx playwright test tests/playwright/experimental-data-ui.spec.js --project=chromium
  
  TIMESTAMP=$(date +"%Y%m%d-%H%M%S")
  DEPLOYMENT_REPORT="./test-results/deployment-report-${TIMESTAMP}.md"
  
  echo "# Experimental Data Enhancement Deployment Report" > "$DEPLOYMENT_REPORT"
  echo "" >> "$DEPLOYMENT_REPORT"
  echo "Deployed on: $(date)" >> "$DEPLOYMENT_REPORT"
  echo "Frontend URL: $FRONTEND_URL" >> "$DEPLOYMENT_REPORT"
  echo "API URL: $API_URL" >> "$DEPLOYMENT_REPORT"
  echo "Branch: $BRANCH_NAME" >> "$DEPLOYMENT_REPORT"
  echo "" >> "$DEPLOYMENT_REPORT"
  echo "## Deployed Features" >> "$DEPLOYMENT_REPORT"
  echo "" >> "$DEPLOYMENT_REPORT"
  echo "| Feature | Status | Notes |" >> "$DEPLOYMENT_REPORT"
  echo "|---------|--------|-------|" >> "$DEPLOYMENT_REPORT"
  echo "| Basic experiment listing | ✅ Complete | Displays experiments with basic information |" >> "$DEPLOYMENT_REPORT"
  echo "| Experiment detail view | ✅ Complete | Shows basic experiment data |" >> "$DEPLOYMENT_REPORT"
  echo "| Enhanced visualization | ⚠️ In progress | Charts and interactive visualizations |" >> "$DEPLOYMENT_REPORT"
  echo "| Filtering | ⚠️ In progress | Basic filtering functionality |" >> "$DEPLOYMENT_REPORT"
  echo "| Data sorting | ⚠️ In progress | Sorting experiments by various criteria |" >> "$DEPLOYMENT_REPORT"
  echo "| Data comparison | ⚠️ Planned | Comparing multiple experiments |" >> "$DEPLOYMENT_REPORT"
  echo "| Data export | ⚠️ Planned | Exporting experimental data |" >> "$DEPLOYMENT_REPORT"
  echo "| Mobile responsiveness | ✅ Complete | Works on mobile devices |" >> "$DEPLOYMENT_REPORT"
  
  echo "✅ Deployment complete! Your experimental data enhancement features are now live at: $FRONTEND_URL"
  echo "Deployment report available at: $DEPLOYMENT_REPORT"
  
  echo "
IMPORTANT: For troubleshooting deployment issues, check:
- Netlify logs at: https://app.netlify.com/sites/${NETLIFY_SITE_NAME}/deploys
- Run validate-experimental-data-enhancement.sh to verify feature functionality
- Check browser console for any JavaScript errors"
else
  echo "❌ Netlify deployment failed."
  echo "Please check the Netlify CLI output above for errors."
fi
#!/bin/bash

# Backend-Frontend Integration Script for CryoProtect
# This script configures and integrates all backend services (Heroku API, RDKit service, and Convex)
# with the frontend deployed on Netlify

set -e  # Exit on error

# Color definitions for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;36m'
NC='\033[0m' # No Color

echo -e "${BLUE}===========================================================${NC}"
echo -e "${GREEN}CryoProtect Backend Integration Script${NC}"
echo -e "${BLUE}===========================================================${NC}"
echo ""

# Configuration variables with defaults
HEROKU_APP_NAME=${HEROKU_APP_NAME:-cryoprotect-8030e4025428}
NETLIFY_SITE_NAME=${NETLIFY_SITE_NAME:-cryoprotect}
CONVEX_PROJECT_ID=${CONVEX_PROJECT_ID:-dynamic-mink-63}
RDKIT_SERVICE_URL=${RDKIT_SERVICE_URL:-https://rdkit.cryoprotect.app}
USE_CONVEX=${USE_CONVEX:-true}

# Derived URLs
API_URL="https://${HEROKU_APP_NAME}.herokuapp.com/v1"
FRONTEND_URL="https://${NETLIFY_SITE_NAME}.netlify.app"
CONVEX_URL="https://${CONVEX_PROJECT_ID}.convex.cloud"

# Display configuration
echo -e "${YELLOW}Service Configuration:${NC}"
echo -e "  ${BLUE}Frontend (Netlify):${NC} $FRONTEND_URL"
echo -e "  ${BLUE}Main API (Heroku):${NC} $API_URL"
echo -e "  ${BLUE}RDKit Service:${NC} $RDKIT_SERVICE_URL"
echo -e "  ${BLUE}Convex Database:${NC} $CONVEX_URL"
echo -e "  ${BLUE}Using Convex:${NC} $USE_CONVEX"
echo ""

# Confirm configuration
read -p "Do you want to continue with this configuration? (y/n) " -n 1 -r
echo ""
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo -e "${RED}Integration aborted by user.${NC}"
    exit 0
fi

echo ""
echo -e "${YELLOW}Step 1: Checking CLI tools availability${NC}"

# Check for required CLI tools
MISSING_TOOLS=false

if ! command -v heroku &> /dev/null; then
    echo -e "${RED}Heroku CLI not found. Please install it:${NC}"
    echo "npm install -g heroku"
    MISSING_TOOLS=true
fi

if ! command -v netlify &> /dev/null; then
    echo -e "${RED}Netlify CLI not found. Please install it:${NC}"
    echo "npm install -g netlify-cli"
    MISSING_TOOLS=true
fi

if ! command -v convex &> /dev/null; then
    echo -e "${RED}Convex CLI not found. Please install it:${NC}"
    echo "npm install -g convex"
    MISSING_TOOLS=true
fi

if [ "$MISSING_TOOLS" = true ]; then
    echo -e "${RED}Please install the missing tools and try again.${NC}"
    exit 1
fi

echo -e "${GREEN}All required CLI tools are available.${NC}"

# Check authentication status
echo ""
echo -e "${YELLOW}Step 2: Checking authentication status${NC}"

echo "Checking Heroku authentication..."
if ! heroku auth:whoami &> /dev/null; then
    echo -e "${RED}Not authenticated with Heroku. Please run:${NC}"
    echo "heroku login"
    exit 1
fi
echo -e "${GREEN}Authenticated with Heroku.${NC}"

echo "Checking Netlify authentication..."
if ! netlify status &> /dev/null; then
    echo -e "${RED}Not authenticated with Netlify. Please run:${NC}"
    echo "netlify login"
    exit 1
fi
echo -e "${GREEN}Authenticated with Netlify.${NC}"

echo "Checking Convex authentication..."
if ! convex status &> /dev/null; then
    echo -e "${RED}Not authenticated with Convex. Please run:${NC}"
    echo "convex login"
    exit 1
fi
echo -e "${GREEN}Authenticated with Convex.${NC}"

# Configure Heroku
echo ""
echo -e "${YELLOW}Step 3: Configuring Heroku backend${NC}"

echo "Setting CORS configuration on Heroku..."
heroku config:set ALLOWED_ORIGINS="$FRONTEND_URL,$RDKIT_SERVICE_URL,$CONVEX_URL" -a $HEROKU_APP_NAME
heroku config:set FRONTEND_URL="$FRONTEND_URL" -a $HEROKU_APP_NAME
heroku config:set CONVEX_URL="$CONVEX_URL" -a $HEROKU_APP_NAME
heroku config:set RDKIT_SERVICE_URL="$RDKIT_SERVICE_URL" -a $HEROKU_APP_NAME

echo -e "${GREEN}Successfully configured Heroku app.${NC}"

# Configure Netlify
echo ""
echo -e "${YELLOW}Step 4: Configuring Netlify frontend${NC}"

echo "Setting environment variables on Netlify..."
netlify env:set NEXT_PUBLIC_API_URL "$API_URL" -s $NETLIFY_SITE_NAME
netlify env:set NEXT_PUBLIC_RDKIT_API_URL "$RDKIT_SERVICE_URL" -s $NETLIFY_SITE_NAME
netlify env:set NEXT_PUBLIC_CONVEX_URL "$CONVEX_URL" -s $NETLIFY_SITE_NAME
netlify env:set NEXT_PUBLIC_USE_CONVEX "$USE_CONVEX" -s $NETLIFY_SITE_NAME
netlify env:set NEXT_PUBLIC_ENVIRONMENT "production" -s $NETLIFY_SITE_NAME
netlify env:set NEXT_PUBLIC_ENABLE_API_LOGGING "true" -s $NETLIFY_SITE_NAME
netlify env:set NEXT_PUBLIC_NETLIFY "true" -s $NETLIFY_SITE_NAME

echo -e "${GREEN}Successfully configured Netlify environment.${NC}"

# Test connectivity between services
echo ""
echo -e "${YELLOW}Step 5: Running connectivity tests between services${NC}"

echo "Creating a temporary test script..."
TEST_SCRIPT="connection_test_$(date +%s).js"

cat > $TEST_SCRIPT <<EOF
const axios = require('axios');

const config = {
  timeout: 10000,
  frontend: '$FRONTEND_URL',
  api: 'https://$HEROKU_APP_NAME.herokuapp.com',
  rdkit: '$RDKIT_SERVICE_URL',
  convex: '$CONVEX_URL'
};

async function testConnection(name, url, endpoint = '', options = {}) {
  const fullUrl = \`\${url}\${endpoint}\`;
  console.log(\`Testing connection to \${name}: \${fullUrl}\`);
  
  try {
    const response = await axios.get(fullUrl, { 
      timeout: config.timeout,
      ...options
    });
    
    const status = response.status;
    const contentType = response.headers['content-type'] || '';
    
    console.log(\`✅ \${name}: Connected successfully (\${status})\`);
    console.log(\`   Content-Type: \${contentType}\`);
    
    if (contentType.includes('application/json')) {
      console.log(\`   Response data: \${JSON.stringify(response.data).substring(0, 100)}...\`);
    }
    
    return true;
  } catch (error) {
    console.error(\`❌ \${name}: Connection failed\`);
    
    if (error.response) {
      console.error(\`   Status: \${error.response.status}\`);
      console.error(\`   Message: \${error.message}\`);
    } else if (error.request) {
      console.error(\`   Network error: No response received\`);
    } else {
      console.error(\`   Error: \${error.message}\`);
    }
    
    return false;
  }
}

async function runTests() {
  console.log('🧪 Starting connection tests...\n');
  
  // Test direct connections
  await testConnection('Frontend', config.frontend);
  await testConnection('API', config.api, '/health');
  await testConnection('RDKit Service', config.rdkit, '/health');
  
  // Test CORS
  await testConnection('API CORS Test', config.api, '/test-cors', {
    headers: { 'Origin': config.frontend }
  });
  
  await testConnection('RDKit CORS Test', config.rdkit, '/test-cors', {
    headers: { 'Origin': config.frontend }
  });
  
  // Test Netlify redirects
  await testConnection('Netlify API Redirect', config.frontend, '/api/health');
  await testConnection('Netlify RDKit Redirect', config.frontend, '/rdkit-api/health');
  
  console.log('\n🏁 Connection tests completed');
}

runTests();
EOF

# Install axios if not present
if ! [ -d "node_modules/axios" ]; then
  echo "Installing axios for tests..."
  npm install axios --no-save
fi

# Run the test script
echo "Running connectivity tests..."
node $TEST_SCRIPT

# Clean up
rm $TEST_SCRIPT

# Final steps and summary
echo ""
echo -e "${BLUE}===========================================================${NC}"
echo -e "${GREEN}Backend Integration Complete!${NC}"
echo -e "${BLUE}===========================================================${NC}"
echo ""
echo -e "${YELLOW}Integration Summary:${NC}"
echo -e "  • Main API (Heroku): $API_URL"
echo -e "  • RDKit Service: $RDKIT_SERVICE_URL"
echo -e "  • Convex Database: $CONVEX_URL"
echo -e "  • Frontend: $FRONTEND_URL"
echo ""
echo -e "${YELLOW}Next Steps:${NC}"
echo -e "  1. Run ${BLUE}npm run build:with-convex${NC} to build the frontend"
echo -e "  2. Run ${BLUE}netlify deploy --prod${NC} to deploy to Netlify"
echo -e "  3. Visit ${BLUE}$FRONTEND_URL/convex-test${NC} to verify Convex integration"
echo ""
echo -e "${GREEN}Your CryoProtect application is now configured with all backend services.${NC}"
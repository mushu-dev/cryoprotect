#!/bin/bash
# Script to set up Heroku environment variables for connecting with Vercel frontend

# Default values
HEROKU_APP_NAME=${HEROKU_APP_NAME:-cryoprotect}
VERCEL_FRONTEND_URL=${VERCEL_FRONTEND_URL:-https://frontend-cryoprotect.vercel.app}

# Check if Heroku CLI is installed
if ! command -v heroku &> /dev/null; then
    echo "Error: Heroku CLI is not installed. Please install it first."
    exit 1
fi

# Check if logged in to Heroku
heroku_auth=$(heroku auth:whoami 2>&1)
if [[ $heroku_auth == *"not logged in"* ]]; then
    echo "You're not logged in to Heroku. Please login first with 'heroku login'."
    exit 1
fi

echo "Setting up environment variables for $HEROKU_APP_NAME..."

# Set required environment variables
heroku config:set \
  VERCEL_FRONTEND_URL="$VERCEL_FRONTEND_URL" \
  FLASK_ENV="production" \
  FLASK_APP="app.py" \
  API_VERSION="1.0.0" \
  ALLOW_CORS="true" \
  --app "$HEROKU_APP_NAME"

echo "Environment variables set successfully."
echo ""
echo "Current Heroku configuration:"
heroku config --app "$HEROKU_APP_NAME"

# Make the script executable
chmod +x setup-heroku-env.sh

echo ""
echo "------------------------------------------------------------"
echo "Next steps:"
echo "1. Deploy your backend to Heroku: heroku git:remote -a $HEROKU_APP_NAME && git push heroku master"
echo "2. Deploy your frontend to Vercel: cd frontend && ./deploy-to-vercel.sh"
echo "3. Test the connection: curl $VERCEL_FRONTEND_URL/api/connectivity-test"
echo "------------------------------------------------------------"
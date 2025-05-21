#!/bin/bash
# Script for clean deployment of just the Next.js frontend

# Base directory
BASE_DIR="/home/mushu/Projects/cryoprotect"
TEMP_DIR="$BASE_DIR/deploy_temp_nextjs"

# Create a temporary deployment directory
echo "Creating temporary deployment directory..."
mkdir -p $TEMP_DIR

# Copy essential frontend files
echo "Copying Next.js files..."
cp -r $BASE_DIR/frontend/src $TEMP_DIR/src
cp -r $BASE_DIR/frontend/public $TEMP_DIR/public
cp $BASE_DIR/frontend/package.json $TEMP_DIR/
cp $BASE_DIR/frontend/next.config.js $TEMP_DIR/
cp $BASE_DIR/frontend/postcss.config.js $TEMP_DIR/
cp $BASE_DIR/frontend/tailwind.config.js $TEMP_DIR/
cp $BASE_DIR/frontend/tsconfig.json $TEMP_DIR/

# Create a simple vercel.json for the deployment
cat > $TEMP_DIR/vercel.json << EOF
{
  "version": 2,
  "buildCommand": "npm install && npm run build",
  "framework": "nextjs",
  "outputDirectory": ".next",
  "env": {
    "NEXT_PUBLIC_API_URL": "https://api.cryoprotect.app/v1",
    "NEXT_PUBLIC_ENVIRONMENT": "production",
    "NEXT_PUBLIC_ENABLE_API_LOGGING": "true",
    "VERCEL_ANALYTICS_ID": "true",
    "VERCEL_SPEED_INSIGHTS": "true",
    "PROTECTION_BYPASS": "TAt23KbtFE8dkZobJU3hpgTP4L5ja07V",
    "NEXT_PUBLIC_FRONTEND_PROTECTION_BYPASS": "TAt23KbtFE8dkZobJU3hpgTP4L5ja07V"
  }
}
EOF

# Change to the deployment directory
cd $TEMP_DIR

echo "Starting Vercel deployment from temporary directory..."
vercel deploy --prod --yes

# Capture the deployment result
DEPLOY_RESULT=$?

# Go back to the original directory
cd $BASE_DIR

# Clean up
echo "Cleaning up temporary files..."
# Uncomment to clean up when finished testing
# rm -rf $TEMP_DIR

# Exit with the deployment result
exit $DEPLOY_RESULT
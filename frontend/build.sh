#!/bin/bash
# Script to build the frontend application for production

# Set environment variables for the build
export NODE_ENV=production
export NEXT_PUBLIC_API_URL=https://api.cryoprotect.app/v1
export NEXTAUTH_URL=https://www.cryoprotect.app

# Generate a random NEXTAUTH_SECRET if not set (for local testing only)
# In production, this should be set to a consistent value
if [ -z "$NEXTAUTH_SECRET" ]; then
  echo "NEXTAUTH_SECRET not set, generating a temporary one for testing..."
  export NEXTAUTH_SECRET=$(openssl rand -base64 32)
  echo "Generated NEXTAUTH_SECRET: $NEXTAUTH_SECRET"
  echo "NOTE: This is for testing only. In production, set a permanent NEXTAUTH_SECRET."
fi

# Clean up previous build artifacts
echo "Cleaning up previous build..."
rm -rf .next

# Install dependencies if needed
if [ "$1" == "--install" ]; then
  echo "Installing dependencies..."
  npm install
fi

# Build the application
echo "Building Next.js application..."
npm run build

# Check build status
if [ $? -eq 0 ]; then
  echo "Build completed successfully!"
  echo "To start the production server: npm run start"
else
  echo "Build failed. Check the errors above."
  exit 1
fi
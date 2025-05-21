#!/bin/bash

# Script to run the application in Convex mode

# Export environment variables
export USE_CONVEX=true
export NEXT_PUBLIC_USE_CONVEX=true

# Run the Convex development service in the background
cd convex
npx convex dev &
CONVEX_PID=$!

# Wait for Convex to initialize
echo "Starting Convex server..."
sleep 5

# Run the Next.js development server with Convex
cd ../frontend
npm run dev:with-convex

# Cleanup: Kill the Convex process when Next.js is stopped
kill $CONVEX_PID
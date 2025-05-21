#!/bin/bash

# Exit on error
set -e

echo "=== Starting CryoProtect Frontend with Convex Enabled ==="
echo ""

# Change to frontend directory
cd frontend

# Set environment variable and run development server
NEXT_PUBLIC_USE_CONVEX=true npm run dev
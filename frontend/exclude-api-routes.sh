#!/bin/bash
# Script to temporarily move API routes during build
set -e

echo "Temporarily excluding API routes for static build..."

# Create a backup directory
mkdir -p .api-routes-backup

# Move API routes to backup
if [ -d "src/app/api" ]; then
  mv src/app/api .api-routes-backup/
  echo "API routes moved to backup"
fi

echo "API routes excluded, ready for static build"
#!/bin/bash
# Script to restore API routes after build
set -e

echo "Restoring API routes after static build..."

# Restore API routes from backup
if [ -d ".api-routes-backup/api" ]; then
  mkdir -p src/app
  mv .api-routes-backup/api src/app/
  echo "API routes restored"
fi

# Remove backup directory
rm -rf .api-routes-backup

echo "API routes restoration complete"
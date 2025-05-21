#!/bin/bash
# Script to commit the Netlify configuration fixes

set -e

echo "Committing Netlify configuration fixes..."

# Navigate to the project root
cd "$(dirname "$0")"

# Add the files to git
git add frontend/netlify.toml
git add deploy-site-fix.sh

# Commit the changes
git commit -m "Fix Netlify deployment configuration

This change:
- Adds explicit base directory setting in netlify.toml to resolve build path issues
- Updates build command to ensure proper build process is used
- Maintains all redirect rules for dynamic routes
- Adds a new deployment script that works with the corrected configuration"

echo "Changes committed successfully!"
echo "To push the changes, run: git push origin experimental-data-enhancement"
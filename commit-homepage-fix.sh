#!/bin/bash
# Script to commit the homepage fix

set -e

echo "Committing homepage fix..."

# Navigate to the project root
cd "$(dirname "$0")"

# Add the files to git
git add frontend/src/pages/index.js
git add frontend/package.json
git add deploy-homepage-fix.sh
git add NETLIFY_DEPLOYMENT_UPDATE.md

# Commit the changes
git commit -m "Fix homepage layout and deployment process

This change:
- Updates the homepage to use a classic vertical layout with sections
- Maintains links to Molecules, Mixtures, Experiments, and Protocols
- Adds appropriate header and footer elements
- Adds the export script to package.json for static site generation
- Updates deployment script to use static export for Netlify
- Adds documentation about the Netlify deployment configuration"

echo "Changes committed successfully!"
echo "To push the changes, run: git push origin experimental-data-enhancement"
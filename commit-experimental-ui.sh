#!/bin/bash
# Script to commit the experimental data enhancement UI implementation

set -e

echo "Committing experimental data enhancement UI implementation..."

# Navigate to the project root
cd "$(dirname "$0")"

# Add the files to git
git add frontend/next.config.js
git add frontend/netlify.toml
git add frontend/src/pages/index.js
git add frontend/src/pages/experiments
git add frontend/src/pages/protocols
git add check-deployment.js
git add deploy-experimental-ui-fix.sh
git add EXPERIMENTAL_DATA_ENHANCEMENT_SUMMARY.md
git add EXPERIMENTAL_DATA_ENHANCEMENT_VERIFICATION.md
git add test-experimental-ui.js
git add direct-ui-test.js

# Commit the changes
git commit -m "Implement experimental data enhancement UI with dynamic routing fixes

This commit:
- Adds feature cards to the homepage for experiments and protocols
- Implements experiments list, detail, and creation pages
- Implements protocols list, detail, and creation pages
- Fixes dynamic routing issues for detail pages
- Updates Netlify configuration for proper deployment
- Adds testing scripts for UI verification
- Includes comprehensive documentation"

echo "Changes committed successfully!"
echo "To push the changes, run: git push origin experimental-data-enhancement"
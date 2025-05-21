#!/bin/bash
set -e

echo "🔒 Fixing npm vulnerabilities..."

# Run npm audit fix with force flag
npm audit fix --force

echo "🔄 Updating potentially vulnerable dependencies..."
npm install --save react@latest react-dom@latest next@latest

# Install specific versions if needed for compatibility
# npm install --save package-name@specific-version

echo "✅ Vulnerability fixes applied!"
echo "🛠️ Run 'npm audit' to verify all issues are resolved."
#!/bin/bash
# Script to prepare and test frontend dependencies locally before deployment

echo "====================================================="
echo "Preparing and testing frontend dependencies locally"
echo "====================================================="

# Navigate to frontend directory
cd frontend || { echo "Error: frontend directory not found"; exit 1; }

# Create .npmrc file with legacy peer deps setting
echo "Creating .npmrc file to handle dependency conflicts..."
cat > .npmrc << 'EOF'
legacy-peer-deps=true
strict-peer-dependencies=false
EOF

# Add engine restrictions to handle Node.js compatibility
echo "Updating package.json with Node.js engine constraints..."
if ! grep -q "\"engines\":" package.json; then
  # Create a temporary file with the updated content
  jq '. + {"engines": {"node": "^18.0.0", "npm": "^9.0.0"}}' package.json > package.json.tmp
  mv package.json.tmp package.json
fi

# Prune any existing node_modules to start clean
echo "Cleaning existing node_modules..."
rm -rf node_modules
rm -f package-lock.json

# Install dependencies using special flags
echo "Installing dependencies with legacy peer deps..."
npm install --legacy-peer-deps --no-fund --no-audit --loglevel=error

# Test if the build works locally
echo "Testing build locally..."
npm run build

# Check build result
if [ $? -eq 0 ]; then
  echo "✅ Local build successful! Dependencies are properly resolved."
  echo "You can now run your Vercel deployment."
else
  echo "❌ Local build failed. This indicates there's still a dependency issue."
  echo "Check the error messages above for more details."
fi

# Return to the original directory
cd ..

# Create a Vercel-specific .npmrc in the root
echo "Creating .npmrc in the project root for Vercel..."
cat > .npmrc << 'EOF'
legacy-peer-deps=true
strict-peer-dependencies=false
EOF

echo "====================================================="
echo "Preparation complete. Try deploying again with:"
echo "vercel deploy"
echo "====================================================="
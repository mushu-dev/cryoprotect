#!/bin/bash
# Script to update package.json with overrides for React versions

echo "Adding React version overrides to package.json..."

# Check if original package.json exists
if [ ! -f "frontend/package.json" ]; then
  echo "Error: frontend/package.json not found"
  exit 1
fi

# Make a backup of the original
cp frontend/package.json frontend/package.json.bak

# Check if overrides file exists
if [ ! -f "frontend/package.json.overrides" ]; then
  echo "Creating overrides file..."
  
  # Extract the content of the current package.json
  content=$(cat frontend/package.json)
  
  # Replace the closing brace with overrides and resolutions sections
  modified_content="${content%\}*}, 
  \"overrides\": {
    \"react\": \"18.2.0\",
    \"react-dom\": \"18.2.0\"
  },
  \"resolutions\": {
    \"react\": \"18.2.0\",
    \"react-dom\": \"18.2.0\"
  }
}"
  
  # Write the modified content to the overrides file
  echo "$modified_content" > frontend/package.json.overrides
fi

# Replace the package.json with the overrides version
cp frontend/package.json.overrides frontend/package.json

# Add .npmrc file for good measure
echo "legacy-peer-deps=true" > frontend/.npmrc

echo "Package.json updated with React version overrides"
echo "Vercel should now be able to install dependencies correctly"
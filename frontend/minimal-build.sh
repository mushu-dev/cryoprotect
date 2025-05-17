#!/bin/bash
# Minimal build script for Netlify
set -ex

# Clean install with yarn which handles conflicts better
echo "Setting up environment..."
npm install -g yarn

# Show versions
echo "Node version: $(node -v)"
echo "Yarn version: $(yarn -v)"

# Install dependencies with yarn
echo "Installing dependencies with yarn..."
yarn install

# Build Next.js
echo "Building Next.js app..."
yarn build

echo "Build completed successfully!"
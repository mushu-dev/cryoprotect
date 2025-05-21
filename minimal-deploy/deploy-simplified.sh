#!/bin/bash

# Super Minimal Deployment Script for Netlify Debug
set -e

echo "=== CryoProtect Ultra Minimal Deployment ==="

# Clean previous builds
rm -rf node_modules .next package-lock.json

# Install dependencies
npm install

# Deploy to Netlify
echo "Deploying to Netlify..."
netlify deploy --build --prod

echo "Deployment complete! Check the site URL above."
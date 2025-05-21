#!/bin/bash
set -e

echo "🔄 Setting up Next.js to use Pages Router only for Netlify deployment..."

# Backup app directory
if [ -d "src/app" ]; then
  echo "📦 Backing up App Router (src/app) directory..."
  mkdir -p .app-router-backup
  cp -r src/app/* .app-router-backup/
  
  # Remove the app directory temporarily
  echo "🗑️ Temporarily removing App Router directory..."
  rm -rf src/app
fi

echo "✅ Successfully configured for Pages Router only deployment!"
echo "🚀 You can now proceed with the build process."
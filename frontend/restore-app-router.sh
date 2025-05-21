#!/bin/bash
set -e

echo "🔄 Restoring App Router configuration..."

# Restore app directory if backup exists
if [ -d ".app-router-backup" ] && [ "$(ls -A .app-router-backup)" ]; then
  echo "📦 Restoring App Router (src/app) directory from backup..."
  mkdir -p src/app
  cp -r .app-router-backup/* src/app/
  
  echo "✅ Successfully restored App Router directory!"
else
  echo "❌ No App Router backup found. Nothing to restore."
fi
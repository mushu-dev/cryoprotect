#!/bin/bash
set -e

echo "🔄 Restoring original navigation component..."

# Restore navigation-header.tsx if backup exists
if [ -f ".typescript-backup/navigation-header.tsx" ]; then
  mv .typescript-backup/navigation-header.tsx src/components/navigation-header.tsx
  echo "✅ Original navigation component restored"
else
  echo "❌ No backup found for navigation-header.tsx"
fi
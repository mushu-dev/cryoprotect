#!/bin/bash
set -e

echo "🔄 Restoring TypeScript components after build..."

# Restore circuit breaker components
if [ -f "src/components/circuit-breaker-dashboard.tsx.disabled" ]; then
  mv src/components/circuit-breaker-dashboard.tsx.disabled src/components/circuit-breaker-dashboard.tsx
fi

if [ -f "src/components/circuit-breaker-indicator.tsx.disabled" ]; then
  mv src/components/circuit-breaker-indicator.tsx.disabled src/components/circuit-breaker-indicator.tsx
fi

if [ -f "src/components/circuit-breaker-status.tsx.disabled" ]; then
  mv src/components/circuit-breaker-status.tsx.disabled src/components/circuit-breaker-status.tsx
fi

# Restore circuit provider component
if [ -f "src/components/circuit-breaker/circuit-provider.tsx.disabled" ]; then
  mv src/components/circuit-breaker/circuit-provider.tsx.disabled src/components/circuit-breaker/circuit-provider.tsx
fi

# Restore analytics provider
if [ -f "src/components/analytics/AnalyticsProvider.tsx.disabled" ]; then
  mv src/components/analytics/AnalyticsProvider.tsx.disabled src/components/analytics/AnalyticsProvider.tsx
fi

# Remove placeholder directory if it exists
if [ -d "src/components/placeholder" ]; then
  rm -rf src/components/placeholder
fi

echo "✅ TypeScript components restored"
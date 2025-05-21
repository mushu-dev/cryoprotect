#!/bin/bash
# Script for creating a public analytics demo

# Base directory
BASE_DIR="/home/mushu/Projects/cryoprotect"
TEMP_DIR="$BASE_DIR/public_analytics_demo"

# Create a temporary deployment directory
echo "Creating public analytics demo directory..."
rm -rf $TEMP_DIR
mkdir -p $TEMP_DIR/pages
mkdir -p $TEMP_DIR/public

# Create a minimal Next.js app
cat > $TEMP_DIR/pages/index.js << 'EOF'
import { Analytics } from '@vercel/analytics/react'
import { SpeedInsights } from '@vercel/speed-insights/next'

export default function Home() {
  return (
    <div style={{ 
      display: 'flex', 
      flexDirection: 'column', 
      alignItems: 'center', 
      justifyContent: 'center', 
      height: '100vh',
      padding: '0 2rem',
      textAlign: 'center'
    }}>
      <h1 style={{ fontSize: '2.5rem', marginBottom: '1rem' }}>CryoProtect Analytics Demo</h1>
      <p style={{ fontSize: '1.2rem', marginBottom: '2rem', maxWidth: '600px' }}>
        This page demonstrates Vercel Analytics and Speed Insights integration.
        It's a minimal deployment to test these features.
      </p>
      <div style={{ 
        background: '#f4f4f4', 
        padding: '1.5rem', 
        borderRadius: '8px',
        maxWidth: '600px'
      }}>
        <h2 style={{ marginTop: 0 }}>Features Enabled:</h2>
        <ul style={{ textAlign: 'left' }}>
          <li>Vercel Analytics</li>
          <li>Vercel Speed Insights</li>
        </ul>
        <p>Visit your Vercel dashboard to see the analytics data.</p>
      </div>
      
      {/* Include Analytics and SpeedInsights */}
      <Analytics />
      <SpeedInsights />
    </div>
  )
}
EOF

# Create simplified package.json
cat > $TEMP_DIR/package.json << 'EOF'
{
  "name": "cryoprotect-analytics-demo",
  "version": "0.1.0",
  "private": true,
  "scripts": {
    "dev": "next dev",
    "build": "next build",
    "start": "next start"
  },
  "dependencies": {
    "@vercel/analytics": "^1.0.1",
    "@vercel/speed-insights": "^1.0.1",
    "next": "14.0.4",
    "react": "18.2.0",
    "react-dom": "18.2.0"
  }
}
EOF

# Create a simple next.config.js
cat > $TEMP_DIR/next.config.js << 'EOF'
/** @type {import('next').NextConfig} */
const nextConfig = {
  reactStrictMode: true,
  swcMinify: true,
}

module.exports = nextConfig
EOF

# Create a specific vercel.json file for public access
cat > $TEMP_DIR/vercel.json << EOF
{
  "version": 2,
  "public": true,
  "buildCommand": "npm install && npm run build",
  "outputDirectory": ".next",
  "framework": "nextjs",
  "env": {
    "VERCEL_ANALYTICS_ID": "true",
    "VERCEL_SPEED_INSIGHTS": "true"
  }
}
EOF

# Change to the deployment directory
cd $TEMP_DIR

echo "Creating a new Vercel project..."
vercel project rm -y analytics-demo 2>/dev/null || true
vercel project add analytics-demo

echo "Starting Vercel deployment from public analytics demo directory..."
vercel --name analytics-demo --confirm --prod

# Capture the deployment result
DEPLOY_RESULT=$?

# Go back to the original directory
cd $BASE_DIR

# Exit with the deployment result
echo "Deployment complete! Check your Vercel dashboard for the URL."
exit $DEPLOY_RESULT
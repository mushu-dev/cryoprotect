#!/bin/bash
# Minimal deployment approach for Vercel

# Environment variables
PROTECTION_BYPASS="TAt23KbtFE8dkZobJU3hpgTP4L5ja07V" 
API_URL="https://cryoprotect-8030e4025428.herokuapp.com/v1"
FRONTEND_URL="https://cryoprotect.vercel.app"
NEXTAUTH_SECRET=${NEXTAUTH_SECRET:-$(openssl rand -base64 32)}

echo "====================================================="
echo "CryoProtect Minimal Vercel Deployment"
echo "====================================================="

# Create a minimal deployment directory
DEPLOY_DIR=$(mktemp -d)
echo "Created deployment directory: $DEPLOY_DIR"

# Create a minimal Next.js app
mkdir -p $DEPLOY_DIR/pages
mkdir -p $DEPLOY_DIR/public

# Create minimal package.json with only what's needed
cat > $DEPLOY_DIR/package.json << 'EOF'
{
  "name": "cryoprotect-minimal",
  "version": "0.1.0",
  "private": true,
  "scripts": {
    "dev": "next dev",
    "build": "next build",
    "start": "next start"
  },
  "dependencies": {
    "@vercel/analytics": "^1.1.1",
    "@vercel/speed-insights": "^1.0.2",
    "next": "13.4.19",
    "react": "18.2.0",
    "react-dom": "18.2.0"
  }
}
EOF

# Create a simple _app.js with Analytics
cat > $DEPLOY_DIR/pages/_app.js << 'EOF'
import { Analytics } from '@vercel/analytics/react'
import { SpeedInsights } from '@vercel/speed-insights/next'

function MyApp({ Component, pageProps }) {
  return (
    <>
      <Component {...pageProps} />
      <Analytics />
      <SpeedInsights />
    </>
  )
}

export default MyApp
EOF

# Create a simple index.js page
cat > $DEPLOY_DIR/pages/index.js << 'EOF'
export default function Home() {
  return (
    <div style={{ 
      padding: '2rem', 
      maxWidth: '800px', 
      margin: '0 auto',
      fontFamily: 'Arial, sans-serif'
    }}>
      <h1>CryoProtect</h1>
      <p>This is a minimal deployment test page with Vercel Analytics enabled.</p>
      
      <div style={{ 
        marginTop: '2rem', 
        padding: '1rem', 
        backgroundColor: '#f5f5f5', 
        borderRadius: '8px' 
      }}>
        <h2>Testing Vercel Deployment</h2>
        <p>If you can see this page, the minimal deployment was successful!</p>
        <p>This page includes:</p>
        <ul>
          <li>Vercel Analytics</li>
          <li>Vercel Speed Insights</li>
          <li>Protection bypass token</li>
        </ul>
      </div>

      <div style={{ 
        marginTop: '2rem', 
        padding: '1rem', 
        backgroundColor: '#e6f7ff', 
        borderRadius: '8px' 
      }}>
        <h2>Next Steps</h2>
        <p>Now that we've successfully deployed a minimal app, we can iterate to add more functionality.</p>
      </div>
    </div>
  )
}
EOF

# Create next.config.js
cat > $DEPLOY_DIR/next.config.js << 'EOF'
/** @type {import('next').NextConfig} */
const nextConfig = {
  reactStrictMode: true,
}

module.exports = nextConfig
EOF

# Create .npmrc
cat > $DEPLOY_DIR/.npmrc << 'EOF'
legacy-peer-deps=true
EOF

# Deploy from our minimal directory
echo "Deploying to Vercel..."
cd $DEPLOY_DIR
vercel deploy --prod --yes \
  --name cryoprotect-minimal \
  -e NEXT_PUBLIC_API_URL=$API_URL \
  -e NEXT_PUBLIC_FRONTEND_PROTECTION_BYPASS="$PROTECTION_BYPASS" \
  -e PROTECTION_BYPASS="$PROTECTION_BYPASS" \
  -e NEXTAUTH_URL=$FRONTEND_URL \
  -e NEXTAUTH_SECRET="$NEXTAUTH_SECRET" \
  -e VERCEL_ANALYTICS_ID=true \
  -e VERCEL_SPEED_INSIGHTS=true

RESULT=$?

# Clean up
echo "Cleaning up deployment directory..."
cd /home/mushu/Projects/cryoprotect
rm -rf $DEPLOY_DIR

if [ $RESULT -eq 0 ]; then
  echo "✅ Minimal deployment successful!"
  echo "This demonstrates that Vercel Analytics is working."
  echo "The minimal deployment has a different URL from your main project."
  echo "You can now work on gradually expanding this deployment."
else
  echo "❌ Deployment failed."
  echo "Please try running with vercel --debug for more information."
fi
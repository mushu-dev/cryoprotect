#!/bin/bash
# Script to deploy a minimal version of the site to Netlify

echo "Creating and deploying minimal version of the site..."

# Change to the frontend directory
cd frontend || { echo "Frontend directory not found"; exit 1; }

# Create a simple index.js file
cat > src/pages/index.js << 'EOL'
import React from 'react';
import Head from 'next/head';

export default function Home() {
  return (
    <div style={{ 
      display: 'flex', 
      flexDirection: 'column', 
      alignItems: 'center', 
      justifyContent: 'center', 
      minHeight: '100vh',
      fontFamily: 'Arial, sans-serif',
      padding: '20px',
      textAlign: 'center'
    }}>
      <Head>
        <title>CryoProtect</title>
        <meta name="description" content="Cryoprotectant research and analysis platform" />
        <link rel="icon" href="/favicon.ico" />
      </Head>

      <h1 style={{ color: '#0070f3', marginBottom: '20px' }}>
        Welcome to CryoProtect
      </h1>

      <p style={{ fontSize: '1.2rem', maxWidth: '600px', lineHeight: 1.6 }}>
        CryoProtect is a platform for cryoprotectant research, analysis, and experimental tracking. 
        We're currently working on enhancements to our platform to better serve the scientific community.
      </p>

      <div style={{ 
        marginTop: '40px',
        padding: '20px',
        border: '1px solid #eaeaea',
        borderRadius: '10px',
        backgroundColor: '#f9f9f9',
        maxWidth: '600px'
      }}>
        <h2 style={{ color: '#0070f3', marginBottom: '10px' }}>Coming Soon</h2>
        <ul style={{ 
          textAlign: 'left', 
          paddingLeft: '20px',
          lineHeight: 1.6
        }}>
          <li>Enhanced molecular visualization</li>
          <li>Improved experimental data tracking</li>
          <li>Protocol and experiment sharing</li>
          <li>Advanced analysis tools</li>
          <li>Integration with laboratory equipment</li>
        </ul>
      </div>

      <footer style={{
        marginTop: '60px',
        color: '#666',
        fontSize: '0.9rem'
      }}>
        © {new Date().getFullYear()} CryoProtect - All rights reserved
      </footer>
    </div>
  );
}
EOL

# Create a basic _app.js file 
cat > src/pages/_app.js << 'EOL'
import React from 'react';
import Head from 'next/head';

function MyApp({ Component, pageProps }) {
  return (
    <>
      <Head>
        <meta name="viewport" content="width=device-width, initial-scale=1" />
      </Head>
      <Component {...pageProps} />
    </>
  );
}

export default MyApp;
EOL

# Update next.config.js to use extremely minimal settings
cat > next.config.js << 'EOL'
/** @type {import('next').NextConfig} */
module.exports = {
  reactStrictMode: true,
  swcMinify: true,
  output: 'export',
  images: {
    unoptimized: true
  },
  eslint: {
    ignoreDuringBuilds: true
  },
  typescript: {
    ignoreBuildErrors: true
  },
  trailingSlash: true
};
EOL

# Update netlify.toml with minimal configuration
cat > netlify.toml << 'EOL'
[build]
  base = "frontend"
  command = "npm run build && npm run export"
  publish = "out"
  
[build.environment]
  NEXT_PUBLIC_ENVIRONMENT = "production"
  NEXT_PUBLIC_NETLIFY = "true"

# These headers help secure your site
[[headers]]
  for = "/*"
  [headers.values]
    Content-Security-Policy = "default-src 'self'; script-src 'self' 'unsafe-eval' 'unsafe-inline'; style-src 'self' 'unsafe-inline'; img-src 'self' data:; font-src 'self';"
    X-Frame-Options = "DENY"
    X-Content-Type-Options = "nosniff"
    Referrer-Policy = "strict-origin-when-cross-origin"
    Permissions-Policy = "camera=(), microphone=(), geolocation=()"
EOL

# Update package.json with export script
npm pkg set scripts.export="next export"

# Run the build 
echo "Running build..."
npm run build

# Check if build was successful
if [ ! -d "out" ]; then
    echo "Build failed - 'out' directory not found"
    npm run export
    if [ ! -d "out" ]; then
        echo "Export failed - 'out' directory still not found"
        exit 1
    fi
fi

# Deploy to Netlify
echo "Deploying to Netlify..."
netlify deploy --prod

echo "Deployment completed!"
#!/bin/bash
set -e

echo "🚀 Running Netlify static export script..."

# Create a directory for static export
mkdir -p netlify-static

# Create simplified pages that don't use React hooks
echo "📝 Creating simplified static pages..."
mkdir -p netlify-static/pages

# Create a static 404 page
cat > netlify-static/pages/404.js << EOL
export default function Custom404() {
  return (
    <div style={{ padding: '40px', textAlign: 'center', maxWidth: '600px', margin: '0 auto' }}>
      <h1 style={{ fontSize: '2rem', marginBottom: '20px' }}>404 - Page Not Found</h1>
      <p style={{ marginBottom: '30px' }}>The page you are looking for does not exist.</p>
      <a 
        href="/"
        style={{ 
          display: 'inline-block', 
          padding: '10px 20px', 
          backgroundColor: '#0070f3', 
          color: 'white', 
          borderRadius: '5px', 
          textDecoration: 'none' 
        }}
      >
        Back to Home
      </a>
    </div>
  )
}
EOL

# Create a static 500 page
cat > netlify-static/pages/500.js << EOL
export default function Custom500() {
  return (
    <div style={{ padding: '40px', textAlign: 'center', maxWidth: '600px', margin: '0 auto' }}>
      <h1 style={{ fontSize: '2rem', marginBottom: '20px' }}>500 - Server Error</h1>
      <p style={{ marginBottom: '30px' }}>Something went wrong on our end. Please try again later.</p>
      <a 
        href="/"
        style={{ 
          display: 'inline-block', 
          padding: '10px 20px', 
          backgroundColor: '#0070f3', 
          color: 'white', 
          borderRadius: '5px', 
          textDecoration: 'none' 
        }}
      >
        Back to Home
      </a>
    </div>
  )
}
EOL

# Create a simplified index page
cat > netlify-static/pages/index.js << EOL
export default function Home() {
  return (
    <div style={{ padding: '40px', textAlign: 'center', maxWidth: '800px', margin: '0 auto' }}>
      <h1 style={{ fontSize: '2.5rem', marginBottom: '20px' }}>CryoProtect</h1>
      <p style={{ marginBottom: '30px' }}>Cryoprotectant Analysis Platform</p>
      <div style={{ display: 'flex', justifyContent: 'center', gap: '20px', flexWrap: 'wrap' }}>
        <a 
          href="/molecules"
          style={{ 
            display: 'inline-block', 
            padding: '10px 20px', 
            backgroundColor: '#0070f3', 
            color: 'white', 
            borderRadius: '5px', 
            textDecoration: 'none' 
          }}
        >
          Explore Molecules
        </a>
        <a 
          href="/mixtures"
          style={{ 
            display: 'inline-block', 
            padding: '10px 20px', 
            backgroundColor: '#0070f3', 
            color: 'white', 
            borderRadius: '5px', 
            textDecoration: 'none' 
          }}
        >
          View Mixtures
        </a>
      </div>
      <p style={{ marginTop: '40px', fontSize: '0.9rem', color: '#666' }}>
        This is a static fallback page. Please enable JavaScript for the full experience.
      </p>
    </div>
  )
}
EOL

# Create simplified placeholder pages for molecules and mixtures static routes
mkdir -p netlify-static/pages/molecules
mkdir -p netlify-static/pages/mixtures

cat > netlify-static/pages/molecules/[id].js << EOL
export default function MoleculePage() {
  return (
    <div style={{ padding: '40px', textAlign: 'center', maxWidth: '600px', margin: '0 auto' }}>
      <h1 style={{ fontSize: '2rem', marginBottom: '20px' }}>Molecule Details</h1>
      <p style={{ marginBottom: '30px' }}>Loading molecule information...</p>
      <p style={{ fontSize: '0.9rem', color: '#666' }}>
        Please enable JavaScript to view molecule details.
      </p>
      <a 
        href="/molecules"
        style={{ 
          display: 'inline-block', 
          marginTop: '20px',
          padding: '10px 20px', 
          backgroundColor: '#0070f3', 
          color: 'white', 
          borderRadius: '5px', 
          textDecoration: 'none' 
        }}
      >
        Back to All Molecules
      </a>
    </div>
  )
}

export function getStaticPaths() {
  return {
    paths: [
      { params: { id: 'placeholder' } },
      { params: { id: '962' } },
      { params: { id: '176' } },
      { params: { id: '6276' } },
      { params: { id: '8857' } }
    ],
    fallback: true
  }
}

export function getStaticProps() {
  return {
    props: {}
  }
}
EOL

cat > netlify-static/pages/mixtures/[id].js << EOL
export default function MixturePage() {
  return (
    <div style={{ padding: '40px', textAlign: 'center', maxWidth: '600px', margin: '0 auto' }}>
      <h1 style={{ fontSize: '2rem', marginBottom: '20px' }}>Mixture Details</h1>
      <p style={{ marginBottom: '30px' }}>Loading mixture information...</p>
      <p style={{ fontSize: '0.9rem', color: '#666' }}>
        Please enable JavaScript to view mixture details.
      </p>
      <a 
        href="/mixtures"
        style={{ 
          display: 'inline-block', 
          marginTop: '20px',
          padding: '10px 20px', 
          backgroundColor: '#0070f3', 
          color: 'white', 
          borderRadius: '5px', 
          textDecoration: 'none' 
        }}
      >
        Back to All Mixtures
      </a>
    </div>
  )
}

export function getStaticPaths() {
  return {
    paths: [
      { params: { id: 'placeholder' } },
      { params: { id: '1' } },
      { params: { id: '2' } },
      { params: { id: '3' } }
    ],
    fallback: true
  }
}

export function getStaticProps() {
  return {
    props: {}
  }
}
EOL

# Create a simplified next.config.js for static export
echo "📝 Creating simplified next.config.js for static export..."
cat > netlify-static/next.config.js << EOL
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
  }
};
EOL

# Create minimal package.json for static export
echo "📝 Creating minimal package.json for static export..."
cat > netlify-static/package.json << EOL
{
  "name": "cryoprotect-frontend-static",
  "version": "0.1.0",
  "private": true,
  "scripts": {
    "dev": "next dev",
    "build": "next build",
    "start": "next start"
  },
  "dependencies": {
    "next": "14.0.4",
    "react": "^18.2.0",
    "react-dom": "^18.2.0"
  }
}
EOL

# Install dependencies and build static export
echo "🔧 Installing dependencies and building static export..."
cd netlify-static
npm install
npm run build

# Create Netlify redirects file
echo "📝 Creating _redirects file for Netlify..."
cat > out/_redirects << EOL
# API routes should redirect to the backend
/api/*  https://cryoprotect-8030e4025428.herokuapp.com/api/:splat  200

# Handle dynamic routes 
/molecules/*  /molecules/[id].html  200
/mixtures/*  /mixtures/[id].html  200

# SPA fallback
/*  /index.html  200
EOL

# Copy static files to original out directory
echo "📦 Copying static files to the original out directory..."
cd ..
mkdir -p out
cp -r netlify-static/out/* out/

echo "✅ Static export completed successfully!"
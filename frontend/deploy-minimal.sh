#!/bin/bash
set -e

echo "🚀 Running minimal Next.js deployment script..."

# Create backup directory
mkdir -p .backup

# Backup TypeScript files
echo "📦 Backing up TypeScript files..."
find src -name "*.tsx" | while read file; do
  mkdir -p ".backup/$(dirname "$file")"
  cp "$file" ".backup/$file"
  # Rename to prevent compilation
  mv "$file" "${file}.bak"
done

# Set up environment variables
echo "📝 Creating .env.local file..."
cat > .env.local << EOL
NEXT_PUBLIC_API_URL=https://cryoprotect-8030e4025428.herokuapp.com/v1
NEXT_PUBLIC_RDKIT_API_URL=https://cryoprotect-rdkit.fly.dev
NEXT_PUBLIC_CONVEX_URL=https://upbeat-parrot-866.convex.cloud
NEXT_PUBLIC_USE_CONVEX=true
NEXT_PUBLIC_ENVIRONMENT=production
NEXT_PUBLIC_ENABLE_API_LOGGING=true
NEXT_PUBLIC_NETLIFY=true
EOL

# Hide the App Router directory
echo "🔧 Hiding App Router directory..."
if [ -d "src/app" ]; then
  mkdir -p .backup/app
  cp -r src/app/* .backup/app/
  rm -rf src/app
fi

# Create simplified next.config.js
echo "📝 Creating simplified next.config.js..."
cp next.config.js .backup/next.config.js
cat > next.config.js << EOL
/** @type {import('next').NextConfig} */
module.exports = {
  reactStrictMode: true,
  swcMinify: true,
  images: {
    unoptimized: true,
  },
  eslint: {
    ignoreDuringBuilds: true
  },
  typescript: {
    ignoreBuildErrors: true
  },
  
  // Static export configuration
  output: 'export',
  trailingSlash: true,
  distDir: 'out'
};
EOL

# Create package.json.override for build
cp package.json .backup/package.json
# Remove all TypeScript dev dependencies from package.json
sed -i '/"typescript"/d' package.json
sed -i '/"@types/react"/d' package.json
sed -i '/"@types/node"/d' package.json

# Create a simple landing page for static export
echo "📝 Creating simple landing page..."
mkdir -p out
cat > out/index.html << EOL
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>CryoProtect</title>
  <style>
    body {
      font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
      line-height: 1.6;
      color: #333;
      background-color: #f8f9fa;
      margin: 0;
      padding: 0;
    }
    header {
      background-color: #0070f3;
      color: white;
      padding: 1rem;
      text-align: center;
    }
    main {
      max-width: 1200px;
      margin: 0 auto;
      padding: 2rem;
    }
    .hero {
      text-align: center;
      margin-bottom: 2rem;
    }
    .hero h1 {
      font-size: 2.5rem;
      color: #0070f3;
    }
    .hero p {
      font-size: 1.25rem;
      color: #666;
    }
    .features {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
      grid-gap: 2rem;
      margin-top: 3rem;
    }
    .feature-card {
      background: white;
      border-radius: 8px;
      box-shadow: 0 4px 6px rgba(0,0,0,0.1);
      padding: 1.5rem;
      transition: transform 0.3s ease;
    }
    .feature-card:hover {
      transform: translateY(-5px);
    }
    .feature-card h2 {
      color: #0070f3;
      margin-top: 0;
    }
    .btn {
      display: inline-block;
      background-color: #0070f3;
      color: white;
      padding: 0.75rem 1.5rem;
      border-radius: 4px;
      text-decoration: none;
      font-weight: bold;
      margin-top: 1rem;
    }
    footer {
      background-color: #f1f1f1;
      text-align: center;
      padding: 1rem;
      margin-top: 3rem;
    }
  </style>
</head>
<body>
  <header>
    <h1>CryoProtect</h1>
  </header>
  
  <main>
    <section class="hero">
      <h1>Cryoprotectant Analysis Platform</h1>
      <p>A comprehensive tool for analyzing and managing cryoprotectant molecules</p>
    </section>
    
    <div class="features">
      <div class="feature-card">
        <h2>Molecule Database</h2>
        <p>Explore our comprehensive database of cryoprotectant molecules and their properties.</p>
        <a href="/molecules" class="btn">View Molecules</a>
      </div>
      
      <div class="feature-card">
        <h2>Mixtures</h2>
        <p>Discover optimized cryoprotectant mixtures and their performance characteristics.</p>
        <a href="/mixtures" class="btn">View Mixtures</a>
      </div>
      
      <div class="feature-card">
        <h2>Experiments</h2>
        <p>Design, track, and analyze cryopreservation experiments with detailed protocols.</p>
        <a href="/experiments" class="btn">Manage Experiments</a>
      </div>
      
      <div class="feature-card">
        <h2>Protocols</h2>
        <p>Create and manage standardized protocols for reproducible cryopreservation procedures.</p>
        <a href="/protocols" class="btn">Browse Protocols</a>
      </div>
    </div>
  </main>
  
  <footer>
    <p>&copy; 2025 CryoProtect. All rights reserved.</p>
  </footer>
</body>
</html>
EOL

# Create Netlify redirects file
echo "📝 Creating _redirects file for Netlify..."
cat > out/_redirects << EOL
# API routes should redirect to the backend
/api/*  https://cryoprotect-8030e4025428.herokuapp.com/v1/:splat  200

# Backend API
/health  https://cryoprotect-8030e4025428.herokuapp.com/health  200

# RDKit API
/rdkit-api/health  https://cryoprotect-rdkit.fly.dev/health  200

# Handle dynamic routes 
/molecules/*  /index.html  200
/mixtures/*  /index.html  200
/experiments/*  /index.html  200
/protocols/*  /index.html  200
/dashboard/*  /index.html  200

# SPA fallback for all other routes
/*  /index.html  200
EOL

# Create other static pages
mkdir -p out/molecules out/mixtures out/experiments out/protocols
cp out/index.html out/molecules/index.html
cp out/index.html out/mixtures/index.html
cp out/index.html out/experiments/index.html
cp out/index.html out/protocols/index.html

echo "✅ Static site created!"
echo "🚀 Ready to deploy to Netlify!"
echo ""
echo "To deploy to Netlify, run:"
echo "netlify deploy --prod --dir=out"

# Restore App Router
echo "🔄 Restoring App Router directory..."
if [ -d ".backup/app" ] && [ "$(ls -A .backup/app)" ]; then
  mkdir -p src/app
  cp -r .backup/app/* src/app/
fi

# Restore next.config.js
mv .backup/next.config.js next.config.js

# Restore package.json
mv .backup/package.json package.json

# Restore TypeScript files
find .backup/src -name "*.tsx.bak" | while read file; do
  target_file=$(echo "$file" | sed 's/\.backup\///' | sed 's/\.bak$//')
  mkdir -p "$(dirname "$target_file")"
  mv "$file" "$target_file"
done
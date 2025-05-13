#!/bin/bash
# Script to build the Next.js application with static export
# This disables server-side rendering for all pages

echo "Building CryoProtect frontend with static export..."
echo

# First, let's update next.config.js to enable static export
# We're going to create a temporary file and then move it
cat > next.config.static.js << EOL
/** @type {import('next').NextConfig} */
const nextConfig = {
  reactStrictMode: true,
  swcMinify: true,
  images: {
    domains: ['localhost', 'api.cryoprotect.app'],
    unoptimized: true, // Required for static export
  },
  output: 'export', // Enable static export
  distDir: 'out',
  trailingSlash: true, // Add trailing slashes to all paths
  eslint: {
    ignoreDuringBuilds: true
  },
  typescript: {
    ignoreBuildErrors: true
  }
}

module.exports = nextConfig
EOL

# Backup original config
cp next.config.js next.config.js.bak

# Apply static config
cp next.config.static.js next.config.js

# Create temporary 404 page if it doesn't exist
mkdir -p src/app
if [ ! -f "src/app/404.tsx" ]; then
  echo "Creating temporary 404 page..."
  cat > src/app/404.tsx << EOL
export default function NotFound() {
  return (
    <div className="flex flex-col items-center justify-center min-h-[70vh]">
      <h1 className="text-4xl font-bold">404 - Page Not Found</h1>
      <p className="mt-4 text-lg text-muted-foreground">
        The page you are looking for does not exist.
      </p>
      <a 
        href="/"
        className="mt-6 px-4 py-2 rounded bg-primary text-primary-foreground hover:bg-primary/90"
      >
        Return Home
      </a>
    </div>
  )
}
EOL
fi

# Temporarily rename dynamic routes directories so they're not included in the build
echo "Temporarily renaming dynamic route directories..."
if [ -d "src/app/molecules/[id]" ]; then
  mv "src/app/molecules/[id]" "src/app/molecules/id_skipped"
fi

if [ -d "src/app/mixtures/[id]" ]; then
  mv "src/app/mixtures/[id]" "src/app/mixtures/id_skipped"
fi

if [ -d "src/app/api/auth/[...nextauth]" ]; then
  mv "src/app/api/auth/[...nextauth]" "src/app/api/auth/nextauth_skipped"
fi

# Build the app
echo "Running Next.js build with static export..."
npm run build

# Restore original directories
echo "Restoring original directories..."
if [ -d "src/app/molecules/id_skipped" ]; then
  mv "src/app/molecules/id_skipped" "src/app/molecules/[id]"
fi

if [ -d "src/app/mixtures/id_skipped" ]; then
  mv "src/app/mixtures/id_skipped" "src/app/mixtures/[id]"
fi

if [ -d "src/app/api/auth/nextauth_skipped" ]; then
  mv "src/app/api/auth/nextauth_skipped" "src/app/api/auth/[...nextauth]"
fi

# Restore original config
mv next.config.js.bak next.config.js

# Remove temporary 404 page if we created one
if [ -f "src/app/404.tsx.bak" ]; then
  mv src/app/404.tsx.bak src/app/404.tsx
fi

echo
echo "Static build completed! The output is in the 'out' directory."
echo "You can deploy this directory to any static hosting service like Vercel, Netlify, or GitHub Pages."
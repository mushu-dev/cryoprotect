#!/bin/bash
# Setup script for Convex integration

echo "Setting up Convex integration for CryoProtect..."

# Copy configuration files if they don't exist
if [ ! -f .env ]; then
    echo "Creating .env from .env.convex template..."
    cp .env.convex .env
    echo "Please update the .env file with your actual Convex deployment key and JWT secret."
fi

if [ ! -f frontend/.env.local ]; then
    echo "Creating frontend/.env.local for Convex integration..."
    mkdir -p frontend
    echo "NEXT_PUBLIC_USE_CONVEX=true" > frontend/.env.local
    echo "NEXT_PUBLIC_CONVEX_URL=https://upbeat-parrot-866.convex.cloud" >> frontend/.env.local
fi

# Install dependencies
echo "Installing dependencies..."
npm install

# Go to frontend directory
cd frontend

# Install Convex dependencies
echo "Installing Convex dependencies..."
npm install convex@latest

# Create _generated directory if it doesn't exist
mkdir -p src/convex/_generated

# Generate Convex types if possible
if command -v npx &> /dev/null; then
    echo "Generating Convex types..."
    npx convex codegen || echo "Warning: Could not generate Convex types. Please run 'npx convex codegen' manually after setting up your Convex project."
else
    echo "Warning: npx not found. Please run 'npx convex codegen' manually after setting up your Convex project."
fi

echo ""
echo "Setup complete. Next steps:"
echo "1. Update the .env file with your actual Convex deployment key and JWT secret"
echo "2. Run './deploy-convex.sh' to deploy your Convex functions"
echo "3. Run 'npm run dev:with-convex' to start the development server with Convex enabled"
echo ""
echo "For production deployment:"
echo "1. Make sure Netlify environment variables are set (already updated in netlify.toml)"
echo "2. Deploy to Netlify: 'cd frontend && npm run build && netlify deploy --prod'"
echo ""
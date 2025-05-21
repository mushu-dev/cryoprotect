#!/bin/bash
# Verification script for Convex integration

echo "Verifying Convex integration..."
echo "==============================="

# Check environment variables
echo "Checking environment variables..."
if [ -f .env ]; then
    if grep -q "CONVEX_DB_ENABLED=true" .env; then
        echo "✅ CONVEX_DB_ENABLED is set to true in .env"
    else
        echo "❌ CONVEX_DB_ENABLED is not set to true in .env"
    fi

    if grep -q "CONVEX_URL" .env; then
        echo "✅ CONVEX_URL is set in .env"
    else
        echo "❌ CONVEX_URL is not set in .env"
    fi
    
    if grep -q "CONVEX_DEPLOYMENT_KEY" .env; then
        echo "✅ CONVEX_DEPLOYMENT_KEY is set in .env"
    else
        echo "❌ CONVEX_DEPLOYMENT_KEY is not set in .env"
    fi
    
    if grep -q "JWT_SECRET" .env; then
        echo "✅ JWT_SECRET is set in .env"
    else
        echo "❌ JWT_SECRET is not set in .env"
    fi
else
    echo "❌ .env file not found. Run ./setup-convex.sh to create it."
fi

# Check frontend environment variables
echo -e "\nChecking frontend environment variables..."
if [ -f frontend/.env.local ]; then
    if grep -q "NEXT_PUBLIC_USE_CONVEX=true" frontend/.env.local; then
        echo "✅ NEXT_PUBLIC_USE_CONVEX is set to true in frontend/.env.local"
    else
        echo "❌ NEXT_PUBLIC_USE_CONVEX is not set to true in frontend/.env.local"
    fi
    
    if grep -q "NEXT_PUBLIC_CONVEX_URL" frontend/.env.local; then
        echo "✅ NEXT_PUBLIC_CONVEX_URL is set in frontend/.env.local"
    else
        echo "❌ NEXT_PUBLIC_CONVEX_URL is not set in frontend/.env.local"
    fi
else
    echo "❌ frontend/.env.local file not found. Run ./setup-convex.sh to create it."
fi

# Check Netlify configuration
echo -e "\nChecking Netlify configuration..."
if [ -f frontend/netlify.toml ]; then
    if grep -q "NEXT_PUBLIC_USE_CONVEX = \"true\"" frontend/netlify.toml; then
        echo "✅ NEXT_PUBLIC_USE_CONVEX is set to true in netlify.toml"
    else
        echo "❌ NEXT_PUBLIC_USE_CONVEX is not set to true in netlify.toml"
    fi
    
    if grep -q "NEXT_PUBLIC_CONVEX_URL" frontend/netlify.toml; then
        echo "✅ NEXT_PUBLIC_CONVEX_URL is set in netlify.toml"
    else
        echo "❌ NEXT_PUBLIC_CONVEX_URL is not set in netlify.toml"
    fi
else
    echo "❌ netlify.toml file not found in frontend directory."
fi

# Check Convex dependencies
echo -e "\nChecking Convex dependencies..."
if grep -q "\"convex\":" frontend/package.json; then
    echo "✅ Convex is listed as a dependency in package.json"
else
    echo "❌ Convex is not listed as a dependency in package.json"
fi

# Check Convex API types
echo -e "\nChecking Convex API types..."
if [ -f frontend/src/convex/_generated/api.js ]; then
    echo "✅ Convex API types file exists"
else
    echo "❌ Convex API types file does not exist. Run 'npx convex codegen' to generate it."
fi

# Check Convex client configuration
echo -e "\nChecking Convex client configuration..."
if [ -f frontend/src/convex/client.ts ]; then
    echo "✅ Convex client configuration exists"
else
    echo "❌ Convex client configuration does not exist."
fi

echo -e "\nVerification complete. If any items are marked with ❌, please fix them before deploying."
echo "Run './setup-convex.sh' to fix most issues automatically."
echo "Run './deploy-convex.sh' to deploy your Convex functions."
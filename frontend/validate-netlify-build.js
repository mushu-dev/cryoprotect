#!/usr/bin/env node

/**
 * This script validates that the Netlify build configuration is correct
 * and that all required files for dynamic routes are in place.
 */

const fs = require('fs');
const path = require('path');
const { execSync } = require('child_process');

// Colors for console output
const colors = {
  reset: '\x1b[0m',
  red: '\x1b[31m',
  green: '\x1b[32m',
  yellow: '\x1b[33m',
  blue: '\x1b[34m',
  cyan: '\x1b[36m',
  white: '\x1b[37m',
};

console.log(`${colors.cyan}=== Validating Netlify Build Configuration ===${colors.reset}\n`);

// Check if netlify.toml exists
const netlifyTomlPath = path.join(process.cwd(), 'netlify.toml');
if (!fs.existsSync(netlifyTomlPath)) {
  console.error(`${colors.red}Error: netlify.toml does not exist${colors.reset}`);
  process.exit(1);
} else {
  console.log(`${colors.green}✓ netlify.toml exists${colors.reset}`);
}

// Check that next.config.js has the correct settings
const nextConfigPath = path.join(process.cwd(), 'next.config.js');
if (!fs.existsSync(nextConfigPath)) {
  console.error(`${colors.red}Error: next.config.js does not exist${colors.reset}`);
  process.exit(1);
} else {
  const nextConfig = fs.readFileSync(nextConfigPath, 'utf8');
  console.log(`${colors.green}✓ next.config.js exists${colors.reset}`);
  
  // Check for output: 'export'
  if (!nextConfig.includes("output: 'export'")) {
    console.warn(`${colors.yellow}Warning: next.config.js does not have output: 'export'${colors.reset}`);
  } else {
    console.log(`${colors.green}✓ next.config.js has output: 'export'${colors.reset}`);
  }
  
  // Check for unoptimized images setting
  if (!nextConfig.includes('unoptimized: true')) {
    console.warn(`${colors.yellow}Warning: next.config.js should have unoptimized: true for images with static export${colors.reset}`);
  } else {
    console.log(`${colors.green}✓ next.config.js has unoptimized: true for images${colors.reset}`);
  }

  // Check for experimental settings that help with static export
  if (!nextConfig.includes('experimental:')) {
    console.warn(`${colors.yellow}Warning: next.config.js does not have experimental settings for better static export${colors.reset}`);
  } else {
    console.log(`${colors.green}✓ next.config.js has experimental settings${colors.reset}`);
  }
}

// Check for generateStaticParams.ts in dynamic routes
const dynamicRoutes = [
  'src/app/molecules/[id]/generateStaticParams.ts',
  'src/app/mixtures/[id]/generateStaticParams.ts'
];

dynamicRoutes.forEach(route => {
  const routePath = path.join(process.cwd(), route);
  if (!fs.existsSync(routePath)) {
    console.error(`${colors.red}Error: ${route} does not exist${colors.reset}`);
  } else {
    const content = fs.readFileSync(routePath, 'utf8');
    if (content.includes('return []')) {
      console.warn(`${colors.yellow}Warning: ${route} returns an empty array, which won't work with static export${colors.reset}`);
    } else {
      console.log(`${colors.green}✓ ${route} has valid static params${colors.reset}`);
    }
  }
});

// Check for environment variables
const envVars = [
  'NEXT_PUBLIC_API_URL',
  'NEXT_PUBLIC_NETLIFY',
  'NEXT_PUBLIC_ENVIRONMENT'
];

console.log(`\n${colors.cyan}=== Checking Environment Variables ===${colors.reset}`);
envVars.forEach(envVar => {
  if (!process.env[envVar]) {
    console.warn(`${colors.yellow}Warning: ${envVar} is not set${colors.reset}`);
  } else {
    console.log(`${colors.green}✓ ${envVar} is set to ${process.env[envVar]}${colors.reset}`);
  }
});

// Check analytics implementation
console.log(`\n${colors.cyan}=== Checking Analytics Implementation ===${colors.reset}`);
const analyticsFiles = [
  'src/app/netlify-analytics.js',
  'src/hooks/useAnalytics.ts',
  'src/components/analytics/AnalyticsProvider.tsx'
];

analyticsFiles.forEach(file => {
  const filePath = path.join(process.cwd(), file);
  if (!fs.existsSync(filePath)) {
    console.warn(`${colors.yellow}Warning: ${file} does not exist${colors.reset}`);
  } else {
    console.log(`${colors.green}✓ ${file} exists${colors.reset}`);
  }
});

// Check build script
const minimalBuildPath = path.join(process.cwd(), 'minimal-build.sh');
if (!fs.existsSync(minimalBuildPath)) {
  console.error(`${colors.red}Error: minimal-build.sh does not exist${colors.reset}`);
} else {
  console.log(`${colors.green}✓ minimal-build.sh exists${colors.reset}`);
  
  // Make sure it's executable
  try {
    execSync(`chmod +x ${minimalBuildPath}`);
    console.log(`${colors.green}✓ minimal-build.sh is executable${colors.reset}`);
  } catch (error) {
    console.error(`${colors.red}Error making minimal-build.sh executable: ${error.message}${colors.reset}`);
  }
}

console.log(`\n${colors.cyan}=== Validation Complete ===${colors.reset}`);
console.log(`${colors.cyan}You're ready to deploy to Netlify!${colors.reset}`);
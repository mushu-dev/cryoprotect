#!/usr/bin/env node

const fs = require('fs');
const path = require('path');
const readline = require('readline');

const rl = readline.createInterface({
  input: process.stdin,
  output: process.stdout
});

// Define the environment file path
const envFilePath = path.join(process.cwd(), '.env.local');

// Check if .env.local already exists
const envFileExists = fs.existsSync(envFilePath);

let existingEnvContent = '';
if (envFileExists) {
  existingEnvContent = fs.readFileSync(envFilePath, 'utf8');
} else {
  // Create empty file if it doesn't exist
  fs.writeFileSync(envFilePath, '');
  console.log(`Created new .env.local file at ${envFilePath}`);
}

console.log('=== Convex Environment Setup ===');
console.log('This script will help you set up your Convex environment variables.');
console.log('Press Enter to keep the current value (shown in parentheses).\n');

// Extract current values
const currentConvexUrl = (existingEnvContent.match(/NEXT_PUBLIC_CONVEX_URL=(.*)/) || [])[1] || '';
const currentUseConvex = (existingEnvContent.match(/NEXT_PUBLIC_USE_CONVEX=(.*)/) || [])[1] || 'false';
const currentClerkKey = (existingEnvContent.match(/NEXT_PUBLIC_CLERK_PUBLISHABLE_KEY=(.*)/) || [])[1] || '';

// Prompt for values
rl.question(`Enter your Convex URL (${currentConvexUrl || 'https://dynamic-mink-63.convex.cloud'}): `, (convexUrl) => {
  const finalConvexUrl = convexUrl || currentConvexUrl || 'https://dynamic-mink-63.convex.cloud';
  
  rl.question(`Enable Convex integration? (${currentUseConvex}): `, (useConvex) => {
    const finalUseConvex = useConvex || currentUseConvex;
    
    rl.question(`Enter your Clerk publishable key (${currentClerkKey || 'not set'}): `, (clerkKey) => {
      const finalClerkKey = clerkKey || currentClerkKey;
      
      // Update or create .env.local file
      let updatedEnvContent = existingEnvContent;
      
      // Update NEXT_PUBLIC_CONVEX_URL
      if (existingEnvContent.includes('NEXT_PUBLIC_CONVEX_URL=')) {
        updatedEnvContent = updatedEnvContent.replace(/NEXT_PUBLIC_CONVEX_URL=.*/, `NEXT_PUBLIC_CONVEX_URL=${finalConvexUrl}`);
      } else {
        updatedEnvContent += `\nNEXT_PUBLIC_CONVEX_URL=${finalConvexUrl}`;
      }
      
      // Update NEXT_PUBLIC_USE_CONVEX
      if (existingEnvContent.includes('NEXT_PUBLIC_USE_CONVEX=')) {
        updatedEnvContent = updatedEnvContent.replace(/NEXT_PUBLIC_USE_CONVEX=.*/, `NEXT_PUBLIC_USE_CONVEX=${finalUseConvex}`);
      } else {
        updatedEnvContent += `\nNEXT_PUBLIC_USE_CONVEX=${finalUseConvex}`;
      }
      
      // Update NEXT_PUBLIC_CLERK_PUBLISHABLE_KEY
      if (finalClerkKey) {
        if (existingEnvContent.includes('NEXT_PUBLIC_CLERK_PUBLISHABLE_KEY=')) {
          updatedEnvContent = updatedEnvContent.replace(/NEXT_PUBLIC_CLERK_PUBLISHABLE_KEY=.*/, `NEXT_PUBLIC_CLERK_PUBLISHABLE_KEY=${finalClerkKey}`);
        } else {
          updatedEnvContent += `\nNEXT_PUBLIC_CLERK_PUBLISHABLE_KEY=${finalClerkKey}`;
        }
      }
      
      // Write the updated content to .env.local
      fs.writeFileSync(envFilePath, updatedEnvContent.trim());
      
      console.log('\nEnvironment variables updated successfully!');
      console.log(`Environment file saved to: ${envFilePath}`);
      
      rl.close();
    });
  });
});
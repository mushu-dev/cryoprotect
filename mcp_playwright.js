// mcp_playwright.js
// This script will set up the environment to allow MCP Playwright to work

const fs = require('fs');
const path = require('path');
const { execSync } = require('child_process');

// Path to symbolic links
const libLinksPath = path.join(__dirname, 'lib_links');

// Check if we need to create the links
if (!fs.existsSync(libLinksPath)) {
  console.log('Creating library symbolic links...');
  fs.mkdirSync(libLinksPath, { recursive: true });
  
  const links = [
    { target: '/usr/lib64/libicudata.so.76', link: 'libicudata.so.66' },
    { target: '/usr/lib64/libicui18n.so.76', link: 'libicui18n.so.66' },
    { target: '/usr/lib64/libicuuc.so.76', link: 'libicuuc.so.66' },
    { target: '/usr/lib64/libjpeg.so.62', link: 'libjpeg.so.8' },
    { target: '/usr/lib64/libwebp.so.7', link: 'libwebp.so.6' },
    { target: '/usr/lib64/libffi.so.8', link: 'libffi.so.7' }
  ];
  
  links.forEach(({ target, link }) => {
    const linkPath = path.join(libLinksPath, link);
    if (!fs.existsSync(linkPath)) {
      try {
        fs.symlinkSync(target, linkPath);
        console.log(`Created symlink ${link} -> ${target}`);
      } catch (error) {
        console.error(`Failed to create symlink ${link}:`, error);
      }
    }
  });
}

// Set the environment variable for child processes
process.env.LD_LIBRARY_PATH = libLinksPath;

// Display success message
console.log('Environment prepared for Playwright');
console.log('LD_LIBRARY_PATH set to:', process.env.LD_LIBRARY_PATH);
console.log('\nInstructions:');
console.log('1. Now start your application with this preloader:');
console.log('   node mcp_playwright.js node your_app.js');
console.log('2. Or export the variable in your shell:');
console.log(`   export LD_LIBRARY_PATH=${libLinksPath}`);
console.log('   Then run your application normally');

// If arguments are provided, execute the command with the new environment
if (process.argv.length > 2) {
  const command = process.argv.slice(2).join(' ');
  console.log(`\nExecuting: ${command}`);
  try {
    execSync(command, { 
      env: process.env, 
      stdio: 'inherit' 
    });
  } catch (error) {
    console.error('Command execution failed:', error);
    process.exit(1);
  }
}
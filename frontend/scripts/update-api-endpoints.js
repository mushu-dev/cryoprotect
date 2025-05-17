// Script to update API endpoints for Netlify deployment
const fs = require('fs');
const path = require('path');

// Define the Heroku backend URL
const HEROKU_BACKEND_URL = 'https://cryoprotect-8030e4025428.herokuapp.com';

// Function to recursively scan and update files
function scanAndUpdateFiles(dir, extensions = ['.js', '.jsx', '.ts', '.tsx']) {
  const files = fs.readdirSync(dir);
  
  for (const file of files) {
    const filePath = path.join(dir, file);
    const stat = fs.statSync(filePath);
    
    if (stat.isDirectory()) {
      // Skip node_modules and .next directories
      if (file !== 'node_modules' && file !== '.next' && file !== '.git') {
        scanAndUpdateFiles(filePath, extensions);
      }
    } else if (extensions.includes(path.extname(file).toLowerCase())) {
      updateApiEndpoints(filePath);
    }
  }
}

// Function to update API endpoints in a file
function updateApiEndpoints(filePath) {
  let content = fs.readFileSync(filePath, 'utf8');
  let modified = false;
  
  // Look for direct Vercel URL references
  const vercelUrlPattern = /https:\/\/[a-zA-Z0-9-]+\.vercel\.app/g;
  if (vercelUrlPattern.test(content)) {
    console.log(`Found Vercel URL in ${filePath}`);
    content = content.replace(vercelUrlPattern, '${process.env.NEXTAUTH_URL}');
    modified = true;
  }
  
  // Look for direct API URL references
  const apiUrlPattern = /https:\/\/api\.cryoprotect\.app/g;
  if (apiUrlPattern.test(content)) {
    console.log(`Found API URL in ${filePath}`);
    content = content.replace(apiUrlPattern, '${process.env.NEXT_PUBLIC_API_URL}');
    modified = true;
  }
  
  // Save the file if modified
  if (modified) {
    console.log(`Updating file: ${filePath}`);
    fs.writeFileSync(filePath, content, 'utf8');
  }
}

// Main execution
console.log('Updating API endpoints for Netlify deployment...');
const srcDir = path.join(__dirname, '..', 'src');
scanAndUpdateFiles(srcDir);
console.log('API endpoint update complete!');
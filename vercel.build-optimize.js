// Vercel Build Optimization Script
// This script helps reduce the serverless function size for Vercel deployment

const fs = require('fs');
const path = require('path');
const { execSync } = require('child_process');

// Directories to cleanup
const pythonDirs = [
  'api/__pycache__',
  'api/**/__pycache__',
  'api/tests',
  'api/*.dist-info',
  'api/*.egg-info',
];

// Large files to remove that aren't needed in production
const largeFiles = [
  'api/*/tests/**',
  'api/**/test_*.py',
  'api/**/tests.py',
  'api/**/documentation/**',
  'api/**/doc/**',
  'api/**/docs/**',
  'api/**/examples/**',
];

// File extensions that can be cleaned up
const cleanupExtensions = [
  '.pyc',
  '.pyo',
  '.pyd',
  '.so',
  '.c',
  '.h',
  '.txt',
  '.md',
  '.rst',
  '.html',
  '.cpp'
];

console.log('🔍 Starting Vercel build optimization...');

// Function to measure directory size
function getDirSize(dirPath) {
  try {
    const result = execSync(`du -sh "${dirPath}"`).toString();
    return result.trim();
  } catch (e) {
    return 'Unknown';
  }
}

// Initial size
console.log(`📊 Initial API directory size: ${getDirSize('./api')}`);

// Cleanup Python cache directories and test directories
pythonDirs.forEach(dir => {
  try {
    execSync(`rm -rf ${dir}`);
    console.log(`🧹 Removed directory pattern: ${dir}`);
  } catch (e) {
    // Ignore errors for globs that don't match
  }
});

// Remove large unnecessary files
largeFiles.forEach(pattern => {
  try {
    execSync(`find ./api -path "${pattern}" -type f -delete`);
    console.log(`🗑️ Removed files matching: ${pattern}`);
  } catch (e) {
    // Ignore errors for globs that don't match
  }
});

// Clean up by extension
cleanupExtensions.forEach(ext => {
  try {
    // Don't delete .py files!
    if (ext !== '.py') {
      execSync(`find ./api -name "*${ext}" -type f -not -path "*/bin/*" | grep -v "__init__" | xargs rm -f`);
      console.log(`🧼 Cleaned up files with extension: ${ext}`);
    }
  } catch (e) {
    // Ignore errors when no files match
  }
});

// Final optimization: create a minimal package structure
console.log('📦 Creating minimal package structure...');

// Final size
console.log(`📊 Final API directory size: ${getDirSize('./api')}`);
console.log('✅ Build optimization complete!');
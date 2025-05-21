/**
 * Static export configuration for Next.js
 * This file helps ensure all necessary UI components are included in the static build
 */

// List of all components that need to be pre-rendered
const uiComponents = [
  'button',
  'dropdown-menu',
  'sheet',
  'dialog',
  'select',
  'toast',
  'tooltip',
];

// Generate entry points for each component 
const componentEntryPoints = {};

uiComponents.forEach(component => {
  componentEntryPoints[`/components/${component}`] = { 
    page: '/components/[component]',
    query: { component }
  };
});

// Main export paths
const mainPaths = {
  '/': { page: '/' },
  '/molecules': { page: '/molecules' },
  '/mixtures': { page: '/mixtures' },
  '/experiments': { page: '/experiments' },
  '/protocols': { page: '/protocols' },
  '/properties': { page: '/properties' },
  '/protocols/create': { page: '/protocols/create' },
};

// Combine all paths
const allPaths = {
  ...mainPaths,
  ...componentEntryPoints
};

module.exports = {
  generateExportPathMap: async () => {
    return allPaths;
  }
};
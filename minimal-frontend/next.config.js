/** @type {import('next').NextConfig} */
module.exports = {
  reactStrictMode: true,
  output: 'export',
  images: {
    unoptimized: true
  },
  // Disable trailing slashes to ensure links work correctly
  trailingSlash: false,
  // Disable source maps in production to reduce bundle size
  productionBrowserSourceMaps: false
};
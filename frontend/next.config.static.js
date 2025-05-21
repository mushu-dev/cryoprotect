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

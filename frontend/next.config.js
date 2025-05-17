/** @type {import('next').NextConfig} */
const nextConfig = {
  reactStrictMode: true,
  swcMinify: true,
  images: {
    domains: ['localhost', 'api.cryoprotect.app'],
    unoptimized: process.env.NODE_ENV !== 'production',
  },
  // API routes are handled by netlify.toml redirects
  async rewrites() {
    return [
      {
        source: '/api/:path*',
        destination: process.env.NEXT_PUBLIC_API_URL + '/:path*',
      },
    ]
  },
  eslint: {
    // We use ESLint on our own workflow before deploying to production
    ignoreDuringBuilds: true
  },
  typescript: {
    // For deployment, we'll ignore TypeScript errors in the build
    ignoreBuildErrors: true
  },
  // Enable static export if you want to deploy as pure static site
  // output: 'export', // Uncomment for static export (removes API routes)
}

module.exports = nextConfig
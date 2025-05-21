/** @type {import('next').NextConfig} */
const nextConfig = {
  reactStrictMode: true,
  swcMinify: true,
  images: {
    domains: ['localhost', 'api.cryoprotect.app'],
    unoptimized: process.env.NODE_ENV !== 'production',
  },
  async rewrites() {
    return [
      {
        source: '/api/:path*',
        destination: process.env.NEXT_PUBLIC_API_URL + '/:path*',
      },
    ]
  },
  eslint: {
    ignoreDuringBuilds: true
  },
  typescript: {
    ignoreBuildErrors: true
  },
  // Force dynamic rendering for protected pages
  pageExtensions: ['js', 'jsx', 'ts', 'tsx'],
  // Add export directive to force static generation where possible
  output: process.env.NEXT_STATIC === 'true' ? 'export' : undefined
}

module.exports = nextConfig

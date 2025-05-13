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
    // We use ESLint on our own workflow before deploying to production
    ignoreDuringBuilds: true
  },
  typescript: {
    // For Vercel deployment, we'll ignore TypeScript errors in the build
    ignoreBuildErrors: true
  },
  // Server Actions are available by default in Next.js 14+
}

module.exports = nextConfig
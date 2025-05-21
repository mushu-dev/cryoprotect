/** @type {import('next').NextConfig} */
const nextConfig = {
  reactStrictMode: true,
  swcMinify: true,
  images: {
    unoptimized: true
  },
  eslint: {
    ignoreDuringBuilds: true
  },
  typescript: {
    ignoreBuildErrors: true
  },
  // This tells Next.js to use trailing slashes in URLs
  trailingSlash: true,
  // This will add environment variables to the browser
  env: {
    NEXT_PUBLIC_ENVIRONMENT: 'production',
    NEXT_PUBLIC_MINIMAL_TEST: 'true'
  }
};

module.exports = nextConfig;
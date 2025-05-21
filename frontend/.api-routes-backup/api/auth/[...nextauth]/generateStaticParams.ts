// This file is needed for static export
// It tells Next.js not to pre-render any auth API routes

export function generateStaticParams() {
  // Return an empty array for static export
  return []
}
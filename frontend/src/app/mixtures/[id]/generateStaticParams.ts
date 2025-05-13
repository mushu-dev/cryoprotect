// This file is used by Next.js for static site generation
// It tells Next.js which dynamic routes to pre-generate at build time

export function generateStaticParams() {
  // For static export, we return an empty array
  // This way, we don't pre-render any mixture pages
  // They will be generated on-demand at runtime client-side
  return []
}
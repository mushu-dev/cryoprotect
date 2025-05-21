// This file is used for static site generation with Next.js

// Generate static params for Netlify static export
export function generateStaticParams() {
  return [
    { id: 'placeholder' },
    { id: '1' },  // Example mixture 1
    { id: '2' },  // Example mixture 2
    { id: '3' }   // Example mixture 3
  ]
}
// This file is used for static site generation with Next.js

// Generate static params for Netlify static export
export function generateStaticParams() {
  return [
    { id: 'placeholder' },
    { id: '962' },   // Glycerol
    { id: '176' },   // DMSO
    { id: '6276' },  // Ethylene Glycol 
    { id: '8857' }   // Propylene Glycol
  ]
}
// This file is required for Next.js to correctly handle dynamic routes in static rendering
// For [...nextauth], we need to provide the possible route patterns

export async function generateStaticParams() {
  return [
    { nextauth: ['signin'] },
    { nextauth: ['signout'] },
    { nextauth: ['callback'] },
    { nextauth: ['error'] },
    { nextauth: ['session'] },
    { nextauth: ['csrf'] },
    { nextauth: ['providers'] },
  ];
}
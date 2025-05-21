// This file is a dummy route to tell Next.js to handle the API directory
// When using static export, API routes are not included

export const dynamic = 'force-static';

export function GET() {
  return new Response(JSON.stringify({ message: 'Using external API endpoints' }), {
    headers: { 'content-type': 'application/json' },
  });
}
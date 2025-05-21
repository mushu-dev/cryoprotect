import { NextResponse } from 'next/server';

export function middleware(request) {
  // Log the request path for debugging
  console.log(`Middleware handling path: ${request.nextUrl.pathname}`);
  
  // Continue to the next middleware or route handler
  return NextResponse.next();
}

// This configures which paths middleware will be invoked on
export const config = {
  matcher: [
    // Skip all internal paths
    '/((?!api|_next/static|_next/image|favicon.ico).*)',
  ],
};
import { NextResponse } from 'next/server'
import type { NextRequest } from 'next/server'
import { getToken } from 'next-auth/jwt'

// Paths that require authentication
const protectedPaths = [
  '/profile',
  '/settings',
  '/mixtures/new',
  '/mixtures/edit',
  '/molecules/import'
]

// Paths that are only accessible to non-authenticated users
const authPaths = [
  '/auth/signin',
  '/auth/signup',
]

export async function middleware(request: NextRequest) {
  const { pathname } = request.nextUrl
  
  // Check if the path is protected
  const isProtectedPath = protectedPaths.some(path => 
    pathname === path || pathname.startsWith(`${path}/`)
  )
  
  // Check if the path is auth-only
  const isAuthPath = authPaths.some(path => 
    pathname === path || pathname.startsWith(`${path}/`)
  )
  
  // Get the session token
  const token = await getToken({ req: request })
  
  // For protected routes, redirect to signin if not authenticated
  if (isProtectedPath && !token) {
    const url = new URL('/auth/signin', request.url)
    url.searchParams.set('callbackUrl', pathname)
    return NextResponse.redirect(url)
  }
  
  // For auth routes, redirect to home if already authenticated
  if (isAuthPath && token) {
    return NextResponse.redirect(new URL('/', request.url))
  }
  
  // Continue for all other routes
  return NextResponse.next()
}

// Configure the middleware to run only on the specified routes
export const config = {
  matcher: [
    // Protected routes
    '/profile/:path*',
    '/settings/:path*',
    '/mixtures/new/:path*',
    '/mixtures/edit/:path*',
    '/molecules/import/:path*',
    // Auth routes
    '/auth/:path*',
  ],
}
import { useState, useEffect } from 'react'
import Link from 'next/link'

export function AuthStatus() {
  const [mounted, setMounted] = useState(false)
  
  // Prevent hydration mismatch
  useEffect(() => {
    setMounted(true)
  }, [])
  
  if (!mounted) {
    return null
  }
  
  // Simplified version for Next.js 12 compatibility
  return (
    <div className="flex items-center gap-2">
      <Link href="/login">
        <a className="px-4 py-2 text-sm font-medium text-gray-700 bg-white border border-gray-300 rounded-md hover:bg-gray-50">
          Sign In
        </a>
      </Link>
      <Link href="/register">
        <a className="px-4 py-2 text-sm font-medium text-white bg-blue-600 rounded-md hover:bg-blue-700">
          Sign Up
        </a>
      </Link>
    </div>
  )
}
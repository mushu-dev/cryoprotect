import Link from 'next/link'
import { Beaker } from 'lucide-react'

export function Footer() {
  const year = new Date().getFullYear()
  
  return (
    <footer className="border-t py-6 md:py-0">
      <div className="container flex flex-col items-center justify-between gap-4 md:h-14 md:flex-row">
        <div className="flex items-center gap-2">
          <Beaker className="h-5 w-5" />
          <p className="text-sm text-muted-foreground">
            © {year} CryoProtect. All rights reserved.
          </p>
        </div>
        
        <nav className="flex gap-4 sm:gap-6">
          <Link href="/about" className="text-sm text-muted-foreground underline-offset-4 hover:underline">
            About
          </Link>
          <Link href="/privacy" className="text-sm text-muted-foreground underline-offset-4 hover:underline">
            Privacy
          </Link>
          <Link href="/terms" className="text-sm text-muted-foreground underline-offset-4 hover:underline">
            Terms
          </Link>
          <Link href="/contact" className="text-sm text-muted-foreground underline-offset-4 hover:underline">
            Contact
          </Link>
        </nav>
      </div>
    </footer>
  )
}
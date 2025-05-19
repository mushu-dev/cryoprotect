import Link from 'next/link'
import { useRouter } from 'next/router'
import { useState } from 'react'

// Simple utility function for class concatenation
const cn = (...classes: string[]) => classes.filter(Boolean).join(' ');

interface NavItem {
  title: string
  href: string
  exact?: boolean
}

export default function NavigationHeader() {
  const router = useRouter();
  const pathname = router.pathname;
  const [mobileNavOpen, setMobileNavOpen] = useState(false);
  
  const navItems: NavItem[] = [
    {
      title: 'Dashboard',
      href: '/',
      exact: true
    },
    {
      title: 'Molecules',
      href: '/molecules'
    },
    {
      title: 'Mixtures',
      href: '/mixtures'
    },
    {
      title: 'Experiments',
      href: '/experiments'
    },
    {
      title: 'Protocols',
      href: '/protocols'
    },
    {
      title: 'Properties',
      href: '/properties'
    }
  ]
  
  const toggleMobileNav = () => {
    setMobileNavOpen(!mobileNavOpen)
  }
  
  return (
    <header className="sticky top-0 z-50 border-b bg-white dark:bg-gray-800">
      <div className="container mx-auto flex h-14 items-center justify-between px-4">
        <div className="flex items-center">
          <Link href="/">
            <a className="mr-6 flex items-center space-x-2 font-bold">
              <span>CryoProtect</span>
            </a>
          </Link>
          
          <nav className="hidden md:flex items-center gap-6">
            {navItems.map((item) => (
              <Link key={item.href} href={item.href}>
                <a className={cn(
                  "text-sm font-medium transition-colors hover:text-gray-600 dark:hover:text-gray-300",
                  item.exact
                    ? pathname === item.href
                      ? "text-black dark:text-white"
                      : "text-gray-500 dark:text-gray-400"
                    : pathname?.startsWith(item.href)
                      ? "text-black dark:text-white"
                      : "text-gray-500 dark:text-gray-400"
                )}>
                  {item.title}
                </a>
              </Link>
            ))}
          </nav>
        </div>
        
        <div className="flex items-center gap-2">
          <button
            className="md:hidden p-2"
            onClick={toggleMobileNav}
            aria-label="Toggle menu"
          >
            <svg xmlns="http://www.w3.org/2000/svg" width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
              <line x1="3" y1="12" x2="21" y2="12"></line>
              <line x1="3" y1="6" x2="21" y2="6"></line>
              <line x1="3" y1="18" x2="21" y2="18"></line>
            </svg>
          </button>
        </div>
      </div>
      
      {/* Mobile navigation */}
      {mobileNavOpen && (
        <div className="fixed inset-0 z-50 bg-white dark:bg-gray-800 md:hidden">
          <div className="container mx-auto flex h-14 items-center justify-between px-4">
            <Link href="/">
              <a className="flex items-center space-x-2 font-bold">
                <span>CryoProtect</span>
              </a>
            </Link>
            
            <button
              className="p-2"
              onClick={toggleMobileNav}
              aria-label="Close menu"
            >
              <svg xmlns="http://www.w3.org/2000/svg" width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
                <line x1="18" y1="6" x2="6" y2="18"></line>
                <line x1="6" y1="6" x2="18" y2="18"></line>
              </svg>
            </button>
          </div>
          
          <nav className="container mx-auto grid gap-6 p-6">
            {navItems.map((item) => (
              <Link key={item.href} href={item.href}>
                <a onClick={toggleMobileNav} className={cn(
                  "flex items-center py-3 text-lg font-medium",
                  item.exact
                    ? pathname === item.href
                      ? "text-black dark:text-white"
                      : "text-gray-500 dark:text-gray-400"
                    : pathname?.startsWith(item.href)
                      ? "text-black dark:text-white"
                      : "text-gray-500 dark:text-gray-400"
                )}>
                  {item.title}
                </a>
              </Link>
            ))}
          </nav>
        </div>
      )}
    </header>
  )
}
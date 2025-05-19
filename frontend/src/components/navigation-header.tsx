'use client'

import Link from 'next/link'
import { usePathname } from 'next/navigation'
import { useState } from 'react'
import { cn } from '@/lib/utils'
import { Button } from '@/components/ui/button'
import { ThemeToggle } from '@/components/theme-toggle'
import { AuthStatus } from '@/features/auth/components/auth-status'
import { Beaker, Database, LineChart, LayoutDashboard, Menu, X, CloudCog, Flask, ScrollText } from 'lucide-react'
import { FlaskConical } from '@/components/ui/icons'

interface NavItem {
  title: string
  href: string
  icon: React.ReactNode
  description?: string
  exact?: boolean
}

export default function NavigationHeader() {
  const pathname = usePathname()
  const [mobileNavOpen, setMobileNavOpen] = useState(false)
  
  const isConvexEnabled = process.env.NEXT_PUBLIC_USE_CONVEX === 'true'
  
  const navItems: NavItem[] = [
    {
      title: 'Dashboard',
      href: '/',
      icon: <LayoutDashboard className="h-4 w-4 mr-2" />,
      description: 'Overview and key stats',
      exact: true
    },
    {
      title: 'Molecules',
      href: '/molecules',
      icon: <Database className="h-4 w-4 mr-2" />,
      description: 'Browse and search molecules'
    },
    {
      title: 'Mixtures',
      href: '/mixtures',
      icon: <Flask className="h-4 w-4 mr-2" />,
      description: 'Analyze cryoprotectant mixtures'
    },
    {
      title: 'Experiments',
      href: '/experiments',
      icon: <FlaskConical className="h-4 w-4 mr-2" />,
      description: 'Design and manage experiments'
    },
    {
      title: 'Protocols',
      href: '/protocols',
      icon: <ScrollText className="h-4 w-4 mr-2" />,
      description: 'Standardized cryopreservation protocols'
    },
    {
      title: 'Properties',
      href: '/properties',
      icon: <LineChart className="h-4 w-4 mr-2" />,
      description: 'Explore molecular properties'
    },
    ...(isConvexEnabled ? [
      {
        title: 'Convex',
        href: '/convex-test',
        icon: <CloudCog className="h-4 w-4 mr-2" />,
        description: 'Test Convex integration'
      }
    ] : [])
  ]
  
  const toggleMobileNav = () => {
    setMobileNavOpen(!mobileNavOpen)
  }
  
  return (
    <header className="sticky top-0 z-50 border-b bg-background/95 backdrop-blur supports-[backdrop-filter]:bg-background/60">
      <div className="container flex h-14 items-center justify-between">
        <div className="flex items-center">
          <Link href="/" className="mr-6 flex items-center space-x-2">
            <Beaker className="h-6 w-6" />
            <span className="font-bold">CryoProtect</span>
          </Link>
          
          <nav className="hidden md:flex items-center gap-6">
            {navItems.map((item) => (
              <Link
                key={item.href}
                href={item.href}
                className={cn(
                  "text-sm font-medium flex items-center transition-colors hover:text-foreground/80",
                  item.exact
                    ? pathname === item.href
                      ? "text-foreground"
                      : "text-foreground/60"
                    : pathname?.startsWith(item.href)
                      ? "text-foreground"
                      : "text-foreground/60"
                )}
              >
                {item.icon}
                {item.title}
              </Link>
            ))}
          </nav>
        </div>
        
        <div className="flex items-center gap-2">
          <ThemeToggle />
          <AuthStatus />
          
          <Button
            variant="ghost"
            size="icon"
            className="md:hidden"
            onClick={toggleMobileNav}
          >
            <Menu className="h-5 w-5" />
            <span className="sr-only">Toggle menu</span>
          </Button>
        </div>
      </div>
      
      {/* Mobile navigation */}
      {mobileNavOpen && (
        <div className="fixed inset-0 z-50 bg-background md:hidden">
          <div className="container flex h-14 items-center justify-between">
            <Link href="/" className="flex items-center space-x-2">
              <Beaker className="h-6 w-6" />
              <span className="font-bold">CryoProtect</span>
            </Link>
            
            <Button
              variant="ghost"
              size="icon"
              onClick={toggleMobileNav}
            >
              <X className="h-5 w-5" />
              <span className="sr-only">Close menu</span>
            </Button>
          </div>
          
          <nav className="container grid gap-6 p-6">
            {navItems.map((item) => (
              <Link
                key={item.href}
                href={item.href}
                onClick={toggleMobileNav}
                className={cn(
                  "flex items-center py-3 text-lg font-medium",
                  item.exact
                    ? pathname === item.href
                      ? "text-foreground"
                      : "text-foreground/60"
                    : pathname?.startsWith(item.href)
                      ? "text-foreground"
                      : "text-foreground/60"
                )}
              >
                {item.icon}
                {item.title}
                
                {item.description && (
                  <span className="ml-2 text-sm text-muted-foreground">
                    {item.description}
                  </span>
                )}
              </Link>
            ))}
          </nav>
        </div>
      )}
    </header>
  )
}
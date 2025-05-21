#!/bin/bash
set -e

echo "🔄 Updating navigation component for build..."

# Backup original navigation-header.tsx
if [ -f "src/components/navigation-header.tsx" ]; then
  cp src/components/navigation-header.tsx .typescript-backup/navigation-header.tsx
  
  # Create a modified version that uses the placeholder
  cat > src/components/navigation-header.tsx << EOL
import React from 'react';
import Link from 'next/link';
import { useState } from 'react';
import { 
  Database, 
  Beaker, 
  Flask, 
  TestTube, 
  Microscope,
  BarChart,
  Settings,
  Menu,
  X,
  ChevronDown
} from 'lucide-react';

import { cn } from '../lib/utils';
import { Button } from './ui/button.jsx';
import {
  DropdownMenu,
  DropdownMenuContent,
  DropdownMenuItem,
  DropdownMenuTrigger,
} from './ui/dropdown-menu.jsx';
import { Sheet, SheetContent, SheetTrigger } from './ui/sheet.jsx';
// Import the placeholder component
import CircuitBreakerIndicator from './placeholder/circuit-breaker-indicator';

export default function NavigationHeader() {
  const [isMobileNavOpen, setIsMobileNavOpen] = useState(false);
  
  // Use a simple string for pathname instead of router.pathname
  const pathname = typeof window !== 'undefined' ? window.location.pathname : '/';
  
  const navItems = [
    {
      title: 'Dashboard',
      href: '/dashboard',
      icon: <BarChart className="h-4 w-4 mr-2" />,
    },
    {
      title: 'Molecules',
      href: '/molecules',
      icon: <Flask className="h-4 w-4 mr-2" />,
    },
    {
      title: 'Mixtures',
      href: '/mixtures',
      icon: <Beaker className="h-4 w-4 mr-2" />,
    },
    {
      title: 'Experiments',
      href: '/experiments',
      icon: <TestTube className="h-4 w-4 mr-2" />,
    },
    {
      title: 'Protocols',
      href: '/protocols',
      icon: <Microscope className="h-4 w-4 mr-2" />,
    },
    {
      title: 'Properties',
      href: '/properties',
      icon: <Database className="h-4 w-4 mr-2" />,
    }
  ];
  
  return (
    <header className="sticky top-0 z-40 border-b bg-background">
      <div className="container mx-auto flex h-16 items-center justify-between px-4">
        <div className="flex items-center gap-6">
          <Link href="/">
            <span className="flex items-center space-x-2 font-bold text-xl text-primary">
              <span>CryoProtect</span>
            </span>
          </Link>
          
          <nav className="hidden md:flex items-center gap-6">
            {navItems.map((item) => {
              const isActive = 
                pathname === item.href || 
                (item.href !== '/' && pathname.startsWith(item.href));
                
              return (
                <Link key={item.href} href={item.href}>
                  <span
                    className={cn(
                      "flex items-center text-sm font-medium transition-colors hover:text-primary",
                      isActive
                        ? "text-primary"
                        : "text-muted-foreground"
                    )}
                  >
                    {item.icon}
                    {item.title}
                  </span>
                </Link>
              );
            })}
          </nav>
        </div>
        
        <div className="flex items-center gap-4">
          {/* Circuit Breaker Status */}
          <div className="hidden md:flex">
            <CircuitBreakerIndicator 
              circuitName="api" 
              showLabel={false}
            />
          </div>
          
          {/* User menu */}
          <DropdownMenu>
            <DropdownMenuTrigger asChild>
              <Button variant="ghost" className="flex items-center gap-1 h-8 w-8 rounded-full">
                <span className="sr-only">Open user menu</span>
                <div className="h-8 w-8 rounded-full bg-primary/10 flex items-center justify-center text-primary font-medium">
                  U
                </div>
              </Button>
            </DropdownMenuTrigger>
            <DropdownMenuContent align="end" className="w-56">
              <div className="flex items-center justify-start gap-2 p-2">
                <div className="flex flex-col space-y-0.5">
                  <p className="text-sm font-medium">User</p>
                  <p className="text-xs text-muted-foreground">user@example.com</p>
                </div>
              </div>
              <DropdownMenuItem>
                <Settings className="mr-2 h-4 w-4" />
                <span>Settings</span>
              </DropdownMenuItem>
              <DropdownMenuItem>
                <span>Sign out</span>
              </DropdownMenuItem>
              <DropdownMenuItem className="flex items-center justify-between">
                <span>System Status</span>
                <CircuitBreakerIndicator 
                  circuitName="api" 
                  showLabel={false} 
                  className="ml-2"
                />
              </DropdownMenuItem>
            </DropdownMenuContent>
          </DropdownMenu>
          
          {/* Mobile navigation */}
          <Sheet>
            <SheetTrigger asChild>
              <Button variant="ghost" size="icon" className="md:hidden">
                <Menu className="h-5 w-5" />
                <span className="sr-only">Toggle menu</span>
              </Button>
            </SheetTrigger>
            <SheetContent side="left" className="w-[300px] sm:max-w-sm">
              <Link href="/">
                <span className="flex items-center space-x-2 font-bold text-xl text-primary mb-8 mt-4">
                  <span>CryoProtect</span>
                </span>
              </Link>
              <nav className="flex flex-col gap-4">
                {navItems.map((item) => {
                  const isActive = 
                    pathname === item.href || 
                    (item.href !== '/' && pathname.startsWith(item.href));
                    
                  return (
                    <Link key={item.href} href={item.href}>
                      <span
                        className={cn(
                          "flex items-center py-2 text-sm font-medium transition-colors hover:text-primary",
                          isActive
                            ? "text-primary"
                            : "text-muted-foreground"
                        )}
                      >
                        {item.icon}
                        {item.title}
                      </span>
                    </Link>
                  );
                })}
                
                {/* Mobile System Status */}
                <Link href="/system-status">
                  <span className="flex items-center justify-between py-2 text-sm font-medium transition-colors hover:text-primary text-muted-foreground">
                    <div className="flex items-center">
                      <Settings className="h-4 w-4 mr-2" />
                      System Status
                    </div>
                    <CircuitBreakerIndicator 
                      circuitName="api" 
                      showLabel={false}
                    />
                  </span>
                </Link>
              </nav>
            </SheetContent>
          </Sheet>
        </div>
      </div>
    </header>
  );
}
EOL
fi

echo "✅ Navigation component updated for build"
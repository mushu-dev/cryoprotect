'use client';

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
  Moon,
  Sun
} from 'lucide-react';

// Simple utility function for class concatenation
const cn = (...classes) => classes.filter(Boolean).join(' ');

export default function NavigationHeader() {
  const [mobileNavOpen, setMobileNavOpen] = useState(false);
  const [currentTheme, setCurrentTheme] = useState('light');
  
  // Get pathname
  const pathname = typeof window !== 'undefined' ? window.location.pathname : '/';
  
  // Toggle theme
  const toggleTheme = () => {
    const newTheme = currentTheme === 'light' ? 'dark' : 'light';
    setCurrentTheme(newTheme);
    if (typeof document !== 'undefined') {
      document.documentElement.classList.toggle('dark');
    }
    if (typeof localStorage !== 'undefined') {
      localStorage.setItem('theme', newTheme);
    }
  };

  const navItems = [
    {
      title: 'Dashboard',
      href: '/',
      exact: true,
      icon: <BarChart className="w-4 h-4 mr-1" />
    },
    {
      title: 'Molecules',
      href: '/molecules',
      icon: <Flask className="w-4 h-4 mr-1" />
    },
    {
      title: 'Mixtures',
      href: '/mixtures',
      icon: <Beaker className="w-4 h-4 mr-1" />
    },
    {
      title: 'Experiments',
      href: '/experiments',
      icon: <TestTube className="w-4 h-4 mr-1" />
    },
    {
      title: 'Protocols',
      href: '/protocols',
      icon: <Microscope className="w-4 h-4 mr-1" />
    },
    {
      title: 'Properties',
      href: '/properties',
      icon: <Database className="w-4 h-4 mr-1" />
    }
  ];
  
  const toggleMobileNav = () => {
    setMobileNavOpen(!mobileNavOpen);
  };
  
  return (
    <header className={cn(
      "sticky top-0 z-50 border-b transition-all duration-200",
      "bg-white dark:bg-gray-800"
    )}>
      <div className="container mx-auto flex h-14 items-center justify-between px-4">
        <div className="flex items-center">
          <Link href="/" className="mr-6 flex items-center space-x-2 font-bold">
            <span>CryoProtect</span>
          </Link>
          
          <nav className="hidden md:flex items-center gap-4 lg:gap-6">
            {navItems.map((item) => (
              <Link 
                key={item.href} 
                href={item.href} 
                className={cn(
                  "text-sm font-medium transition-colors hover:text-primary flex items-center px-1 py-1.5 relative group",
                  item.exact
                    ? pathname === item.href
                      ? "text-primary"
                      : "text-gray-500 dark:text-gray-400"
                    : pathname?.startsWith(item.href)
                      ? "text-primary"
                      : "text-gray-500 dark:text-gray-400"
                )}
              >
                {item.icon}
                <span className="lg:ml-1">{item.title}</span>
                <span className={cn(
                  "absolute inset-x-0 -bottom-px h-0.5 bg-primary scale-x-0 group-hover:scale-x-100 transition-transform",
                  (item.exact ? pathname === item.href : pathname?.startsWith(item.href)) && "scale-x-100"
                )}></span>
              </Link>
            ))}
          </nav>
        </div>
        
        <div className="flex items-center gap-3">
          <button
            onClick={toggleTheme}
            className="p-2 rounded-full hover:bg-gray-100 dark:hover:bg-gray-700 focus:outline-none"
            aria-label="Toggle theme"
          >
            {currentTheme === 'dark' ? 
              <Sun className="h-5 w-5" /> : 
              <Moon className="h-5 w-5" />}
          </button>
          
          <button
            className="md:hidden p-2 rounded-full hover:bg-gray-100 dark:hover:bg-gray-700"
            onClick={toggleMobileNav}
            aria-label="Toggle menu"
          >
            <svg xmlns="http://www.w3.org/2000/svg" width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
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
            <Link href="/" className="flex items-center space-x-2 font-bold">
              <span>CryoProtect</span>
            </Link>
            
            <button
              className="p-2 rounded-full hover:bg-gray-100 dark:hover:bg-gray-700"
              onClick={toggleMobileNav}
              aria-label="Close menu"
            >
              <X className="h-5 w-5" />
            </button>
          </div>
          
          <nav className="container mx-auto mt-4 pb-20">
            <div className="space-y-1 px-3">
              {navItems.map((item) => (
                <Link 
                  key={item.href} 
                  href={item.href} 
                  onClick={toggleMobileNav} 
                  className={cn(
                    "flex items-center py-3 px-4 text-base font-medium rounded-md transition-colors",
                    item.exact
                      ? pathname === item.href
                        ? "bg-primary/10 text-primary"
                        : "text-gray-700 dark:text-gray-300 hover:bg-gray-100 dark:hover:bg-gray-700/30"
                      : pathname?.startsWith(item.href)
                        ? "bg-primary/10 text-primary"
                        : "text-gray-700 dark:text-gray-300 hover:bg-gray-100 dark:hover:bg-gray-700/30"
                  )}
                >
                  <span className="flex items-center justify-center w-8 h-8">
                    {item.icon}
                  </span>
                  <span className="ml-3">{item.title}</span>
                </Link>
              ))}
            </div>
            
            <div className="flex justify-center mt-8">
              <button
                onClick={() => {
                  toggleTheme();
                  // Don't close mobile nav
                }}
                className="flex items-center px-4 py-2 rounded-md bg-gray-100 dark:bg-gray-700"
              >
                {currentTheme === 'dark' ? 
                  <>
                    <Sun className="h-5 w-5 mr-2" />
                    <span>Light Mode</span>
                  </> : 
                  <>
                    <Moon className="h-5 w-5 mr-2" />
                    <span>Dark Mode</span>
                  </>}
              </button>
            </div>
          </nav>
        </div>
      )}
    </header>
  );
}
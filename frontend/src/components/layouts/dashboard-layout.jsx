import React, { useState } from 'react';
import Link from 'next/link';
import { useRouter } from 'next/router';
import { 
  Database, 
  Beaker, 
  Flask, 
  TestTube, 
  Microscope,
  BarChart,
  Settings,
  ChevronRight,
  Menu,
  X
} from 'lucide-react';
import { cn } from '../../lib/utils';

// Navigation items with icons
const navItems = [
  {
    title: 'Dashboard',
    href: '/dashboard',
    icon: <BarChart className="h-5 w-5" />,
  },
  {
    title: 'Molecules',
    href: '/molecules',
    icon: <Flask className="h-5 w-5" />,
  },
  {
    title: 'Mixtures',
    href: '/mixtures',
    icon: <Beaker className="h-5 w-5" />,
  },
  {
    title: 'Experiments',
    href: '/experiments',
    icon: <TestTube className="h-5 w-5" />,
  },
  {
    title: 'Protocols',
    href: '/protocols',
    icon: <Microscope className="h-5 w-5" />,
  },
  {
    title: 'Properties',
    href: '/properties',
    icon: <Database className="h-5 w-5" />,
  },
  {
    title: 'Settings',
    href: '/settings',
    icon: <Settings className="h-5 w-5" />,
  },
];

export default function DashboardLayout({ children }) {
  const router = useRouter();
  const [sidebarOpen, setSidebarOpen] = useState(false);

  return (
    <div className="flex h-screen bg-gray-100 dark:bg-gray-900">
      {/* Mobile sidebar backdrop */}
      {sidebarOpen && (
        <div 
          className="fixed inset-0 z-40 bg-gray-600 bg-opacity-75 lg:hidden"
          onClick={() => setSidebarOpen(false)}
        />
      )}

      {/* Sidebar */}
      <aside
        className={cn(
          "fixed inset-y-0 left-0 z-50 w-64 transform bg-white dark:bg-gray-800 shadow-lg transition-transform duration-300 ease-in-out lg:relative lg:translate-x-0",
          sidebarOpen ? "translate-x-0" : "-translate-x-full"
        )}
      >
        {/* Sidebar header */}
        <div className="flex h-16 items-center justify-between px-4">
          <Link href="/" className="flex items-center">
            <span className="text-lg font-bold text-primary">CryoProtect</span>
          </Link>
          <button
            className="lg:hidden rounded p-1 hover:bg-gray-100 dark:hover:bg-gray-700"
            onClick={() => setSidebarOpen(false)}
          >
            <X className="h-5 w-5" />
          </button>
        </div>

        {/* Sidebar navigation */}
        <nav className="space-y-1 px-2 py-4">
          {navItems.map((item) => {
            const isActive = 
              router.pathname === item.href || 
              (item.href !== '/' && router.pathname.startsWith(item.href));
              
            return (
              <Link 
                key={item.href} 
                href={item.href}
                className={cn(
                  "group flex items-center rounded-md px-3 py-2 text-sm font-medium",
                  isActive
                    ? "bg-primary/10 text-primary"
                    : "text-gray-700 hover:bg-gray-100 hover:text-gray-900 dark:text-gray-300 dark:hover:bg-gray-700 dark:hover:text-white"
                )}
              >
                {item.icon}
                <span className="ml-3">{item.title}</span>
                {isActive && (
                  <ChevronRight className="ml-auto h-4 w-4 text-primary" />
                )}
              </Link>
            );
          })}
        </nav>
      </aside>

      {/* Main content */}
      <div className="flex flex-1 flex-col overflow-hidden">
        {/* Top navigation */}
        <header className="bg-white shadow dark:bg-gray-800">
          <div className="px-4 py-3 flex items-center justify-between">
            <button
              className="lg:hidden rounded-md p-2 text-gray-600 hover:bg-gray-100 hover:text-gray-900 dark:text-gray-400 dark:hover:bg-gray-700 dark:hover:text-white"
              onClick={() => setSidebarOpen(true)}
            >
              <Menu className="h-6 w-6" />
            </button>
            <div className="flex-1 text-xl font-semibold px-4">
              {/* Page title would go here */}
            </div>
            <div className="flex items-center gap-x-4">
              {/* User profile, notifications, etc. would go here */}
            </div>
          </div>
        </header>

        {/* Page content */}
        <main className="flex-1 overflow-auto">
          <div className="container mx-auto px-4 py-8">
            {children}
          </div>
        </main>
      </div>
    </div>
  );
}